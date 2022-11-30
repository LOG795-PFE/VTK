/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGPUImageData.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkGPUImageData.h"
#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkOpenGLRenderWindow.h"
#include "vtkTextureObject.h"
#include <vtkCellType.h>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkGPUImageData);

vtkInformationKeyMacro(vtkGPUImageData, CONTEXT_OBJECT, ObjectBase);

//----------------------------------------------------------------------------
vtkGPUImageData::vtkGPUImageData()
{
  this->DataDescription = VTK_XYZ_GRID;

  for (int idx = 0; idx < 3; ++idx)
  {
    this->Dimensions[idx] = 0;
    this->Origin[idx] = 0.0;
    this->Spacing[idx] = 1.0;
  }

  this->DirectionMatrix = vtkMatrix3x3::New();
  this->IndexToPhysicalMatrix = vtkMatrix4x4::New();
  this->PhysicalToIndexMatrix = vtkMatrix4x4::New();
  this->DirectionMatrix->Identity();
  this->ComputeTransforms();

  int extent[6] = { 0, -1, 0, -1, 0, -1 };
  memcpy(this->Extent, extent, 6 * sizeof(int));

  this->Information->Set(vtkDataObject::DATA_EXTENT_TYPE(), VTK_3D_EXTENT);
  this->Information->Set(vtkDataObject::DATA_EXTENT(), this->Extent, 6);

  this->TextureObject = nullptr;
}

//----------------------------------------------------------------------------
vtkGPUImageData::~vtkGPUImageData()
{
  if (this->DirectionMatrix)
  {
    this->DirectionMatrix->Delete();
  }
  if (this->IndexToPhysicalMatrix)
  {
    this->IndexToPhysicalMatrix->Delete();
  }
  if (this->PhysicalToIndexMatrix)
  {
    this->PhysicalToIndexMatrix->Delete();
  }
}

//----------------------------------------------------------------------------
// DataSet implementation
void vtkGPUImageData::CopyStructure(vtkDataSet* ds)
{
  vtkGPUImageData* sPts = static_cast<vtkGPUImageData*>(ds);
  this->SetDimensions(0, 0, 0); // this->Initialize();

  for (int idx = 0; idx < 3; ++idx)
  {
    this->Dimensions[idx] = sPts->Dimensions[idx];
    this->Spacing[idx] = sPts->Spacing[idx];
    this->Origin[idx] = sPts->Origin[idx];
  }
  this->DirectionMatrix->DeepCopy(sPts->GetDirectionMatrix());
  this->ComputeTransforms();
  this->SetExtent(sPts->GetExtent());
}

//------------------------------------------------------------------------------
// DataSet implementation
vtkIdType vtkGPUImageData::GetNumberOfCells()
{
  vtkIdType nCells = 1;
  const int* dims = this->GetDimensions();

  for (int i = 0; i < 3; i++)
  {
    if (dims[i] == 0)
    {
      return 0;
    }
    if (dims[i] > 1)
    {
      nCells *= (static_cast<long long>(dims[i]) - 1);
    }
  }

  return nCells;
}

//----------------------------------------------------------------------------
// DataSet implementation
void vtkGPUImageData::GetPoint(vtkIdType ptId, double x[3]) {
  const int* extent = this->Extent;

  vtkIdType dims[3];
  this->GetDimensions(dims);

  x[0] = x[1] = x[2] = 0.0;
  if (dims[0] == 0 || dims[1] == 0 || dims[2] == 0)
  {
    vtkErrorMacro("Requesting a point from an empty image.");
    return;
  }

  // "loc" holds the point x,y,z indices
  int loc[3];
  loc[0] = loc[1] = loc[2] = 0;

  switch (this->DataDescription)
  {
    case VTK_EMPTY:
      return;

    case VTK_SINGLE_POINT:
      break;

    case VTK_X_LINE:
      loc[0] = ptId;
      break;

    case VTK_Y_LINE:
      loc[1] = ptId;
      break;

    case VTK_Z_LINE:
      loc[2] = ptId;
      break;

    case VTK_XY_PLANE:
      loc[0] = ptId % dims[0];
      loc[1] = ptId / dims[0];
      break;

    case VTK_YZ_PLANE:
      loc[1] = ptId % dims[1];
      loc[2] = ptId / dims[1];
      break;

    case VTK_XZ_PLANE:
      loc[0] = ptId % dims[0];
      loc[2] = ptId / dims[0];
      break;

    case VTK_XYZ_GRID:
      loc[0] = ptId % dims[0];
      loc[1] = (ptId / dims[0]) % dims[1];
      loc[2] = ptId / (dims[0] * dims[1]);
      break;
  }

  int i, j, k;
  i = loc[0] + extent[0];
  j = loc[1] + extent[2];
  k = loc[2] + extent[4];
  this->TransformIndexToPhysicalPoint(i, j, k, x);
}

//----------------------------------------------------------------------------
// DataSet implementation
vtkCell* vtkGPUImageData::GetCell(vtkIdType cellId)
{
  vtkErrorMacro("TODO");
  return nullptr;
}

//----------------------------------------------------------------------------
// DataSet implementation
vtkCell* vtkGPUImageData::GetCell(int iMin, int jMin, int kMin)
{
  vtkErrorMacro("TODO");
  return nullptr;
}

//----------------------------------------------------------------------------
// DataSet implementation
void vtkGPUImageData::GetCell(vtkIdType cellId, vtkGenericCell* cell)
{
  vtkErrorMacro("TODO");
}

//------------------------------------------------------------------------------
vtkIdType vtkGPUImageData::FindPoint(double x[3])
{
  vtkErrorMacro("TODO");
  return -1;
}

//----------------------------------------------------------------------------
void vtkGPUImageData::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

}

//----------------------------------------------------------------------------
// Set dimensions of structured points dataset.
void vtkGPUImageData::SetDimensions(int i, int j, int k)
{
  this->Dimensions[0] = i;
  this->Dimensions[1] = j;
  this->Dimensions[2] = k;

  this->SetExtent(0, i - 1, 0, j - 1, 0, k - 1);
}

//----------------------------------------------------------------------------
// Set dimensions of structured points dataset.
void vtkGPUImageData::SetDimensions(const int dim[3])
{
  this->SetDimensions(dim[0], dim[1], dim[2]);
}

//------------------------------------------------------------------------------
// Set dimensions of structured points dataset.
int vtkGPUImageData::GetCellType(vtkIdType vtkNotUsed(cellId)) {
  vtkErrorMacro("TODO");
  return VTK_EMPTY_CELL;
}

//------------------------------------------------------------------------------
vtkIdType vtkGPUImageData::FindCell(double x[3], vtkCell* vtkNotUsed(cell),
  vtkIdType vtkNotUsed(cellId), double tol2, int& subId, double pcoords[3], double* weights)
{
  vtkErrorMacro("TODO");
  return -1;
}

vtkIdType vtkGPUImageData::FindCell(double x[3], vtkCell* vtkNotUsed(cell),
  vtkGenericCell* vtkNotUsed(gencell), vtkIdType vtkNotUsed(cellId), double tol2, int& subId,
  double pcoords[3], double* weights)
{
  vtkErrorMacro("TODO");
  return -1;
}

//------------------------------------------------------------------------------
void vtkGPUImageData::SetDirectionMatrix(vtkMatrix3x3* m)
{
  vtkMTimeType lastModified = this->GetMTime();
  vtkSetObjectBodyMacro(DirectionMatrix, vtkMatrix3x3, m);
  if (lastModified < this->GetMTime())
  {
    this->ComputeTransforms();
  }
}

//------------------------------------------------------------------------------
void vtkGPUImageData::ComputeTransforms()
{
  vtkMatrix4x4* m4 = vtkMatrix4x4::New();
  if (this->DirectionMatrix->IsIdentity())
  {
    m4->Zero();
    m4->SetElement(0, 0, this->Spacing[0]);
    m4->SetElement(1, 1, this->Spacing[1]);
    m4->SetElement(2, 2, this->Spacing[2]);
    m4->SetElement(3, 3, 1);
  }
  else
  {
    const double* m3 = this->DirectionMatrix->GetData();
    m4->SetElement(0, 0, m3[0] * this->Spacing[0]);
    m4->SetElement(0, 1, m3[1] * this->Spacing[1]);
    m4->SetElement(0, 2, m3[2] * this->Spacing[2]);
    m4->SetElement(1, 0, m3[3] * this->Spacing[0]);
    m4->SetElement(1, 1, m3[4] * this->Spacing[1]);
    m4->SetElement(1, 2, m3[5] * this->Spacing[2]);
    m4->SetElement(2, 0, m3[6] * this->Spacing[0]);
    m4->SetElement(2, 1, m3[7] * this->Spacing[1]);
    m4->SetElement(2, 2, m3[8] * this->Spacing[2]);
    m4->SetElement(3, 0, 0);
    m4->SetElement(3, 1, 0);
    m4->SetElement(3, 2, 0);
    m4->SetElement(3, 3, 1);
  }
  m4->SetElement(0, 3, this->Origin[0]);
  m4->SetElement(1, 3, this->Origin[1]);
  m4->SetElement(2, 3, this->Origin[2]);

  this->IndexToPhysicalMatrix->DeepCopy(m4);
  vtkMatrix4x4::Invert(m4, this->PhysicalToIndexMatrix);
  m4->Delete();
}

//------------------------------------------------------------------------------
void vtkGPUImageData::SetDirectionMatrix(const double elements[9])
{
  this->SetDirectionMatrix(elements[0], elements[1], elements[2], elements[3], elements[4],
    elements[5], elements[6], elements[7], elements[8]);
}

//------------------------------------------------------------------------------
void vtkGPUImageData::SetDirectionMatrix(double e00, double e01, double e02, double e10, double e11,
  double e12, double e20, double e21, double e22)
{
  vtkMatrix3x3* m3 = this->DirectionMatrix;
  vtkMTimeType lastModified = m3->GetMTime();

  m3->SetElement(0, 0, e00);
  m3->SetElement(0, 1, e01);
  m3->SetElement(0, 2, e02);
  m3->SetElement(1, 0, e10);
  m3->SetElement(1, 1, e11);
  m3->SetElement(1, 2, e12);
  m3->SetElement(2, 0, e20);
  m3->SetElement(2, 1, e21);
  m3->SetElement(2, 2, e22);

  if (lastModified < m3->GetMTime())
  {
    this->ComputeTransforms();
    this->Modified();
  }
}

//------------------------------------------------------------------------------
template <typename T1, typename T2>
inline static void TransformCoordinates(
  T1 input0, T1 input1, T1 input2, T2 output[3], vtkMatrix4x4* m4)
{
  double* mdata = m4->GetData();
  output[0] = mdata[0] * input0 + mdata[1] * input1 + mdata[2] * input2 + mdata[3];
  output[1] = mdata[4] * input0 + mdata[5] * input1 + mdata[6] * input2 + mdata[7];
  output[2] = mdata[8] * input0 + mdata[9] * input1 + mdata[10] * input2 + mdata[11];
}

//------------------------------------------------------------------------------
void vtkGPUImageData::TransformContinuousIndexToPhysicalPoint(double i, double j, double k,
  double const origin[3], double const spacing[3], double const direction[9], double xyz[3])
{
  for (int c = 0; c < 3; ++c)
  {
    xyz[c] = i * spacing[0] * direction[c * 3] + j * spacing[1] * direction[c * 3 + 1] +
      k * spacing[2] * direction[c * 3 + 2] + origin[c];
  }
}

//------------------------------------------------------------------------------
void vtkGPUImageData::TransformIndexToPhysicalPoint(int i, int j, int k, double xyz[3])
{
  TransformCoordinates<int, double>(i, j, k, xyz, this->IndexToPhysicalMatrix);
}

//------------------------------------------------------------------------------
void vtkGPUImageData::TransformIndexToPhysicalPoint(const int ijk[3], double xyz[3])
{
  TransformCoordinates<int, double>(ijk[0], ijk[1], ijk[2], xyz, this->IndexToPhysicalMatrix);
}

//----------------------------------------------------------------------------
int *vtkGPUImageData::GetDimensions()
{
  this->GetDimensions(this->Dimensions);
  return this->Dimensions;
}

//----------------------------------------------------------------------------
void vtkGPUImageData::GetDimensions(int *dOut)
{
  const int* extent = this->Extent;
  dOut[0] = extent[1] - extent[0] + 1;
  dOut[1] = extent[3] - extent[2] + 1;
  dOut[2] = extent[5] - extent[4] + 1;
}

//----------------------------------------------------------------------------
void vtkGPUImageData::GetDimensions(vtkIdType dims[3])
{
  // Use vtkIdType to avoid overflow on large images
  const int* extent = this->Extent;
  dims[0] = extent[1] - extent[0] + 1;
  dims[1] = extent[3] - extent[2] + 1;
  dims[2] = extent[5] - extent[4] + 1;
}

//----------------------------------------------------------------------------
void vtkGPUImageData::SetContext(vtkOpenGLRenderWindow* renWin)
{
  if (this->TextureObject)
  {
    this->TextureObject->SetContext(renWin);
  }
  this->Context = renWin;
}

//----------------------------------------------------------------------------
vtkOpenGLRenderWindow* vtkGPUImageData::GetContext()
{
  if (!this->TextureObject)
  {
    vtkErrorMacro("TODO");
    return nullptr;
  }

  return this->TextureObject->GetContext();
}

//----------------------------------------------------------------------------
bool vtkGPUImageData::AllocateScalars(int dataType, int numComponents)
{
  int dimensions[3] = { 0,0,0 };
  this->GetDimensions(dimensions);
  if (this->TextureObject && (
    this->TextureObject->GetWidth() != dimensions[0] ||
    this->TextureObject->GetHeight() != dimensions[1] ||
    this->TextureObject->GetDepth() != dimensions[2] ||
    this->TextureObject->GetVTKDataType() != dataType ||
    this->TextureObject->GetComponents() != numComponents ))
  {
    this->TextureObject = nullptr;
  }

  if (!this->TextureObject)
  {
    this->TextureObject = vtkSmartPointer<vtkTextureObject>::New();
    this->TextureObject->SetContext(this->Context);
    bool textureInt = false;
    if (!this->TextureObject->Create3D(dimensions[0], dimensions[1], dimensions[2], numComponents, dataType, false))
    {
      this->TextureObject = nullptr;
      return false;
    }
  }
  return true;
}

//----------------------------------------------------------------------------
bool vtkGPUImageData::AllocateScalarsFromPointer(int dataType, int numComponents, void *data)
{
  int dimensions[3] = { 0,0,0 };
  this->GetDimensions(dimensions);

  this->TextureObject = vtkSmartPointer<vtkTextureObject>::New();
  this->TextureObject->SetContext(this->Context);
  if (!this->TextureObject->Create3DFromRaw(dimensions[0], dimensions[1], dimensions[2], numComponents, dataType, data))
  {
    this->TextureObject = nullptr;
    return false;
  }
  return true;
}

//----------------------------------------------------------------------------
bool vtkGPUImageData::AllocateScalarsFromObject(int dataType, int numComponents, vtkDataObject* object)
{
  int dimensions[3] = { 0, 0, 0 };
  this->GetDimensions(dimensions);

  auto image = vtkImageData::SafeDownCast(object);

  this->TextureObject = vtkSmartPointer<vtkTextureObject>::New();
  this->TextureObject->SetContext(this->Context);
  if (!this->TextureObject->Create3DFromRaw(
        dimensions[0], dimensions[1], dimensions[2], numComponents, dataType, image->GetScalarPointer()))
  {
    this->TextureObject = nullptr;
    return false;
  }

  vtkIdType imageSize = dimensions[0] * dimensions[1] * dimensions[2];
  vtkDataArray* scalars = this->PointData->GetScalars();
  if (scalars && scalars->GetDataType() == dataType && scalars->GetReferenceCount() == 1)
  {
    scalars->SetNumberOfComponents(numComponents);
    scalars->SetNumberOfTuples(imageSize);
    // Since the execute method will be modifying the scalars
    // directly.
    scalars->Modified();
    return false;
  }

  // allocate the new scalars
  scalars = vtkDataArray::CreateDataArray(dataType);
  scalars->SetNumberOfComponents(numComponents);
  scalars->SetName("ImageScalars");

  // allocate enough memory
  scalars->SetNumberOfTuples(imageSize);

  this->PointData->SetScalars(scalars);
  scalars->Delete();

  return true;
}

//----------------------------------------------------------------------------
int vtkGPUImageData::GetScalarType()
{
  if (!this->TextureObject)
  {
    vtkErrorMacro("TODO");
    return -1;
  }
  return this->TextureObject->GetVTKDataType();
}

//----------------------------------------------------------------------------
void vtkGPUImageData::CopyInformationFromPipeline(vtkInformation* info)
{
  if (info->Has(vtkDataObject::SPACING()))
  {
    this->SetSpacing(info->Get(vtkDataObject::SPACING()));
  }
  if (info->Has(vtkDataObject::ORIGIN()))
  {
    this->SetOrigin(info->Get(vtkDataObject::ORIGIN()));
  }
  if (info->Has(vtkGPUImageData::CONTEXT_OBJECT()))
  {
    vtkOpenGLRenderWindow* contextObject = vtkOpenGLRenderWindow::SafeDownCast(info->Get(vtkGPUImageData::CONTEXT_OBJECT()));
    this->SetContext(contextObject);
  }
}

//----------------------------------------------------------------------------
void vtkGPUImageData::CopyInformationToPipeline(vtkInformation* info)
{
  info->Set(vtkDataObject::SPACING(), this->Spacing, 3);
  info->Set(vtkDataObject::ORIGIN(), this->Origin, 3);
  if (this->GetContext())
  {
    info->Set(vtkGPUImageData::CONTEXT_OBJECT(), this->GetContext());
  }
}