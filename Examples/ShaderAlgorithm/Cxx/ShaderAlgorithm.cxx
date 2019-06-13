#include <vtkSmartPointer.h>
#include <vtkImageViewer2.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkShaderProgram.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkImageMapToWindowLevelColors.h>
#include <vtkLookupTable.h>
#include <vtkTextureObject.h>
#include <vtkOpenGLShaderProperty.h>
#include <vtkUniforms.h>

#include <vtkImageToGPUImageFilter.h>
#include <vtkGPUSimpleImageFilter.h>
#include <vtkGPUImageToImageFilter.h>
#include <vtkGPUImageGaussianFilter.h>

#include <string>
#include <fstream>

#include <vtkNrrdReader.h>
#include <vtkMetaImageWriter.h>

int main(int argc, char* argv[])
{
  if (argc < 1)
  {
    return 0;
  }

  // Shader for generating a box
  std::string cFragShader = R"(
//VTK::System::Dec
varying vec2 tcoordVSOutput;
uniform float zPos;
//VTK::AlgTexUniforms::Dec
//VTK::CustomUniforms::Dec
//VTK::Output::Dec
void main(void) {
  float total = floor(tcoordVSOutput.x*float(outputSize.x/boxSize)) +
                floor(tcoordVSOutput.y*float(outputSize.y/boxSize)) + 
                floor(zPos*float(outputSize.z/boxSize));
  if (mod(total,2.0) == 0.0)
  {
    gl_FragData[0] = vec4(1.0);
  }
  else
  {
    gl_FragData[0] = vec4(vec3(0.0), 1.0);
  }
}
)";

  // Shader for generating a box
  std::string sFragShader = R"(
//VTK::System::Dec
varying vec2 tcoordVSOutput;
uniform float zPos;
//VTK::AlgTexUniforms::Dec
//VTK::CustomUniforms::Dec
//VTK::Output::Dec
void main(void) {
  vec4 value1 = texture3D(inputTex0, vec3(tcoordVSOutput.x, tcoordVSOutput.y, zPos));
  vec4 value2 = texture3D(inputTex1, vec3(tcoordVSOutput.x, tcoordVSOutput.y, zPos));
  gl_FragData[0] = vec4(vec3(value1.r * value2.r), 1.0);
}
)";

  vtkNew<vtkNrrdReader> reader;
  reader->SetFileName(argv[1]);
  reader->Update();
  vtkImageData* inputImage = reader->GetOutput();

  vtkNew<vtkImageToGPUImageFilter> inputConvert;
  inputConvert->SetInputDataObject(inputImage);
  inputConvert->Update();

  vtkNew<vtkGPUSimpleImageFilter> checkerPatternGenerator;
  checkerPatternGenerator->GetShaderProperty()->SetFragmentShaderCode(cFragShader.c_str());
  checkerPatternGenerator->SetOutputExtent(inputImage->GetExtent());
  checkerPatternGenerator->SetRenderWindow(inputConvert->GetRenderWindow());
  vtkOpenGLShaderProperty* property = checkerPatternGenerator->GetShaderProperty();
  property->GetFragmentCustomUniforms()->SetUniformi("boxSize", 10);

  vtkNew<vtkGPUImageGaussianFilter> gaussianAlgorithm;
  gaussianAlgorithm->SetOutputScalarTypeToShort();
  gaussianAlgorithm->AddInputConnection(inputConvert->GetOutputPort());

  vtkNew<vtkGPUSimpleImageFilter> shaderAlgorithm;
  shaderAlgorithm->GetShaderProperty()->SetFragmentShaderCode(sFragShader.c_str());
  shaderAlgorithm->AddInputConnection(gaussianAlgorithm->GetOutputPort());
  shaderAlgorithm->AddInputConnection(checkerPatternGenerator->GetOutputPort());
  shaderAlgorithm->SetOutputScalarTypeToShort();

  vtkNew<vtkGPUImageToImageFilter> outputConvert;
  outputConvert->SetInputConnection(shaderAlgorithm->GetOutputPort());
  outputConvert->Update();
  vtkSmartPointer<vtkImageData> outputImage = outputConvert->GetOutput();

  // Visualize
  vtkSmartPointer<vtkImageViewer2> imageViewer =
    vtkSmartPointer<vtkImageViewer2>::New();
  imageViewer->SetInputData(outputImage);
  imageViewer->GetRenderWindow()->SetSize(500, 500);
  imageViewer->GetWindowLevel()->SetWindow(imageViewer->GetInput()->GetScalarRange()[1] - imageViewer->GetInput()->GetScalarRange()[0]);
  imageViewer->GetWindowLevel()->SetLevel((imageViewer->GetInput()->GetScalarRange()[1] - imageViewer->GetInput()->GetScalarRange()[0]) / 2.0);
  imageViewer->GetRenderer()->ResetCamera();
  imageViewer->SetSlice(50);

  // Set up an interactor that does not respond to mouse events
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetInteractorStyle(0);
  imageViewer->GetRenderWindow()->SetInteractor(renderWindowInteractor);
  imageViewer->Render();

  // Start the event loop
  renderWindowInteractor->Initialize();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}