/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestWarpScalarGenerateEnclosure.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkDataSet.h"
#include "vtkExplicitStructuredGrid.h"
#include "vtkFloatArray.h"
#include "vtkLookupTable.h"
#include "vtkMath.h"
#include "vtkPlaneSource.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRTAnalyticSource.h"
#include "vtkRegressionTestImage.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkTypeList.h"
#include "vtkWarpScalar.h"

#include <array>
#include <cstdlib>

namespace
{

void AddScalarAttributeToOutput(vtkAlgorithm* algo)
{
  algo->Update();
  vtkDataSet* output = vtkDataSet::SafeDownCast(algo->GetOutputDataObject(0));
  if (!output)
  {
    std::cout << "Fail: failure to cast data set when adding scalar warping attribute" << std::endl;
  }
  vtkNew<vtkFloatArray> warping;
  warping->SetName("Warp");
  warping->SetNumberOfComponents(1);
  warping->SetNumberOfTuples(output->GetNumberOfPoints());
  for (int iP = 0; iP < output->GetNumberOfPoints(); iP++)
  {
    double point[3] = { 0.0 };
    output->GetPoint(iP, point);
    warping->SetValue(iP, vtkMath::Norm(point) + 1.0);
  }
  output->GetPointData()->AddArray(warping);
  output->GetPointData()->SetActiveScalars("Warp");
}
}

int TestWarpScalarGenerateEnclosure(int argc, char* argv[])
{
  //----------------------------------------------------------------------
  // PolyData
  //----------------------------------------------------------------------
  vtkNew<vtkPlaneSource> planeSrc;
  planeSrc->SetResolution(7, 7);
  ::AddScalarAttributeToOutput(planeSrc);
  vtkNew<vtkWarpScalar> warper;
  warper->SetInputConnection(planeSrc->GetOutputPort());
  warper->GenerateEnclosureOn();
  warper->SetScaleFactor(0.5);
  warper->Update();
  vtkPolyData* output = vtkPolyData::SafeDownCast(warper->GetOutputDataObject(0));
  if (!output)
  {
    std::cout << "Did not output a poly data for plane transformation" << std::endl;
    return EXIT_FAILURE;
  }
  vtkNew<vtkLookupTable> surfaceLUT;
  surfaceLUT->SetRange(output->GetPointData()->GetScalars()->GetRange());
  surfaceLUT->Build();
  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(warper->GetOutputPort());
  mapper->ScalarVisibilityOn();
  mapper->SetScalarRange(output->GetPointData()->GetScalars()->GetRange());
  mapper->SetLookupTable(surfaceLUT);
  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  actor->GetProperty()->SetOpacity(1.0);
  vtkNew<vtkRenderer> renderer;
  renderer->AddActor(actor);
  vtkNew<vtkRenderWindow> renWin;
  renWin->AddRenderer(renderer);
  renWin->Render();

  vtkCamera* camera = renderer->GetActiveCamera();
  camera->SetPosition(9, 9, 9);
  renderer->ResetCamera();

  return (vtkRegressionTester::Test(argc, argv, renWin, 10) == vtkRegressionTester::PASSED)
    ? EXIT_SUCCESS
    : EXIT_FAILURE;
}
