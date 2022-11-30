#include <vtkCamera.h>
#include <vtkColorTransferFunction.h>
#include <vtkContourValues.h>
#include <vtkGPUImageGaussianFilter.h>
#include <vtkGPUImageToImageFilter.h>
#include <vtkGPUSimpleImageFilter.h>
#include <vtkImageData.h>
#include <vtkImageToGPUImageFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkMetaImageReader.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOpenGLGPUVolumeRayCastMapper.h>
#include <vtkOpenGLShaderProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkUniforms.h>
#include <vtkVolume.h>
#include <vtkVolumeProperty.h>

int main(int argc, char* argv[])
{
  double iso1 = 500.0;
  double iso2 = 1150.0;

  // Shader for generating a box
  std::string cFragShader = R"(
//VTK::System::Dec
varying vec2 tcoordVSOutput;
uniform float zPos;
//VTK::AlgTexUniforms::Dec
//VTK::CustomUniforms::Dec
//VTK::Output::Dec
void main(void) {
  float total = floor(tcoordVSOutput.x*float(outputSize0.x/boxSize)) +
                floor(tcoordVSOutput.y*float(outputSize0.y/boxSize)) + 
                floor(zPos*float(outputSize0.z/boxSize));
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

  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " file.mnd [iso1=500] [iso2=1150]" << std::endl;
    std::cout << "e.g. FullHead.mhd 500 1150" << std::endl;
    return EXIT_FAILURE;
  }

  // Setup Pipeline
  vtkNew<vtkMetaImageReader> reader;
  reader->SetFileName(argv[1]);
  reader->Update();
  vtkImageData* inputImage = reader->GetOutput();

  vtkNew<vtkImageToGPUImageFilter> inputConvert;
  inputConvert->SetInputDataObject(inputImage);

  vtkNew<vtkGPUSimpleImageFilter> checkerPatternGenerator;
  checkerPatternGenerator->GetShaderProperty()->SetFragmentShaderCode(cFragShader.c_str());
  checkerPatternGenerator->SetOutputExtent(inputImage->GetExtent());
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

  vtkNew<vtkOpenGLGPUVolumeRayCastMapper> mapper;
  mapper->SetInputConnection(shaderAlgorithm->GetOutputPort());
  mapper->AutoAdjustSampleDistancesOff();
  mapper->SetSampleDistance(0.5);
  mapper->SetBlendModeToIsoSurface();

  // Setup Volume properties for isosurface colors
  if (argc > 3)
  {
    iso1 = atof(argv[2]);
    iso2 = atof(argv[3]);
  }

  vtkNew<vtkNamedColors> colors;
  vtkNew<vtkColorTransferFunction> colorTransferFunction;
  colorTransferFunction->RemoveAllPoints();
  colorTransferFunction->AddRGBPoint(iso2, colors->GetColor3d("ivory").GetData()[0],
    colors->GetColor3d("ivory").GetData()[1], colors->GetColor3d("ivory").GetData()[2]);
  colorTransferFunction->AddRGBPoint(iso1, colors->GetColor3d("flesh").GetData()[0],
    colors->GetColor3d("flesh").GetData()[1], colors->GetColor3d("flesh").GetData()[2]);

  vtkNew<vtkPiecewiseFunction> scalarOpacity;
  scalarOpacity->AddPoint(iso1, .3);
  scalarOpacity->AddPoint(iso2, 0.6);

  vtkNew<vtkVolumeProperty> volumeProperty;
  volumeProperty->ShadeOn();
  volumeProperty->SetInterpolationTypeToLinear();
  volumeProperty->SetColor(colorTransferFunction);
  volumeProperty->SetScalarOpacity(scalarOpacity);

  // Add some contour values to draw iso surfaces
  volumeProperty->GetIsoSurfaceValues()->SetValue(0, iso1);
  volumeProperty->GetIsoSurfaceValues()->SetValue(1, iso2);

  vtkNew<vtkVolume> volume;
  volume->SetMapper(mapper);
  volume->SetProperty(volumeProperty);

  // Setup Renderer
  vtkNew<vtkRenderer> renderer;
  renderer->AddVolume(volume);
  renderer->SetBackground(colors->GetColor3d("cornflower").GetData());

  // Setup RenderWindow
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->SetSize(800, 600);
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("RayCastIsosurface");

  // Set the context of our GPUImageData
  // Those will execute an offscreen render
  inputConvert->SetRenderWindow(renderWindow);
  gaussianAlgorithm->SetRenderWindow(renderWindow);
  checkerPatternGenerator->SetRenderWindow(renderWindow);
  shaderAlgorithm->SetRenderWindow(renderWindow);

  // Setup interactor
  vtkNew<vtkRenderWindowInteractor> interactor;
  interactor->SetRenderWindow(renderWindow);
  vtkNew<vtkInteractorStyleTrackballCamera> style;
  interactor->SetInteractorStyle(style);

  /* Start the renderWindow to init OpenGL before ResetCamera
  since we need it to create Textures, necessary to compute
  the volume's bounding box and properly execute ResetCamera()
  */
  renderWindow->Start();

  // Setup a camera
  vtkNew<vtkCamera> aCamera;
  renderer->SetActiveCamera(aCamera);

  // Generate a good view
  aCamera->SetViewUp(0, -1, 0);
  renderer->ResetCamera();
  renderer->ResetCameraClippingRange();

  // Start interactor
  interactor->Start();

  return EXIT_SUCCESS;
}
