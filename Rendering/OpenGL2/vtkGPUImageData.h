/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGPUImageData.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkGPUImageData
 * @brief   abstracts an OpenGL texture object.
 *
 * vtkGPUImageData represents an OpenGL texture object. It provides API to
 * create textures using data already loaded into pixel buffer objects. It can
 * also be used to create textures without uploading any data.
*/

#ifndef vtkGPUImageData_h
#define vtkGPUImageData_h

#include "vtkRenderingOpenGL2Module.h" // For export macro
#include "vtkDataSet.h"
#include "vtkInformationObjectBaseKey.h"
#include "vtkSmartPointer.h"
#include "vtkTextureObject.h"

class vtkOpenGLRenderWindow;

// \ingroup InformationKeys

class VTKRENDERINGOPENGL2_EXPORT vtkGPUImageData : public vtkDataSet
{
public:

  static vtkInformationObjectBaseKey* CONTEXT_OBJECT();

  static vtkGPUImageData* New();
  vtkTypeMacro(vtkGPUImageData, vtkDataObject);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void CopyStructure(vtkDataSet* ds) override;

///@{
  /**
   * Standard vtkDataSet API methods. See vtkDataSet for more information.
   */
  vtkIdType GetNumberOfCells() override;
  vtkIdType GetNumberOfPoints() override;
  double* GetPoint(vtkIdType ptId) VTK_SIZEHINT(3) override;
  void GetPoint(vtkIdType id, double x[3]) override;
  vtkCell* GetCell(vtkIdType cellId) override;
  vtkCell* GetCell(int i, int j, int k) override;
  void GetCell(vtkIdType cellId, vtkGenericCell* cell) override;
  vtkIdType FindPoint(double x[3]) override;
  vtkIdType FindCell(double x[3], vtkCell* cell, vtkIdType cellId, double tol2, int& subId,
    double pcoords[3], double* weights) override;
  vtkIdType FindCell(double x[3], vtkCell* cell, vtkGenericCell* gencell, vtkIdType cellId,
    double tol2, int& subId, double pcoords[3], double* weights) override;
  int GetCellType(vtkIdType cellId) override;
  using vtkDataSet::GetCellPoints;
  void GetCellPoints(vtkIdType cellId, vtkIdList* ptIds) override
  {
    int dimensions[3];
    this->GetDimensions(dimensions);
    //vtkStructuredData::GetCellPoints(cellId, ptIds, this->DataDescription, dimensions);
  }
  void GetPointCells(vtkIdType ptId, vtkIdList* cellIds) override
  {
    int dimensions[3];
    this->GetDimensions(dimensions);
    //vtkStructuredData::GetPointCells(ptId, cellIds, dimensions);
  }
  int GetMaxCellSize() override { return 6; }
  ///@}

  /**
  * Same as SetExtent(0, i-1, 0, j-1, 0, k-1)
  */
  virtual void SetDimensions(int i, int j, int k);

  /**
  * Same as SetExtent(0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1)
  */
  virtual void SetDimensions(const int dims[3]);

  /**
  * Get dimensions of this structured points dataset.
  * It is the number of points on each axis.
  * Dimensions are computed from Extents during this call.
  * \warning Non thread-safe, use second signature if you want it to be.
  */
  virtual int *GetDimensions() VTK_SIZEHINT(3);

  /**
  * Get dimensions of this structured points dataset.
  * It is the number of points on each axis.
  * This method is thread-safe.
  * \warning The Dimensions member variable is not updated during this call.
  */
  virtual void GetDimensions(int dims[3]);
  //@}

  //@{
  /**
  * Set the spacing (width,height,length) of the cubical cells that
  * compose the data set.
  */
  vtkSetVector3Macro(Spacing, double);
  vtkGetVector3Macro(Spacing, double);
  //@}

  //@{
  /**
  * Set/Get the origin of the dataset. The origin is the position in world
  * coordinates of the point of extent (0,0,0). This point does not have to be
  * part of the dataset, in other words, the dataset extent does not have to
  * start at (0,0,0) and the origin can be outside of the dataset bounding
  * box.
  * The origin plus spacing determine the position in space of the points.
  */
  vtkSetVector3Macro(Origin, double);
  vtkGetVector3Macro(Origin, double);
  //@}

  //@{
  /**
  * Set/Get the extent. On each axis, the extent is defined by the index
  * of the first point and the index of the last point.  The extent should
  * be set before the "Scalars" are set or allocated.  The Extent is
  * stored in the order (X, Y, Z).
  * The dataset extent does not have to start at (0,0,0). (0,0,0) is just the
  * extent of the origin.
  * The first point (the one with Id=0) is at extent
  * (Extent[0],Extent[2],Extent[4]). As for any dataset, a data array on point
  * data starts at Id=0.
  */
  vtkSetVector6Macro(Extent, int);
  vtkGetVector6Macro(Extent, int);
  //@}

  void SetContext(vtkOpenGLRenderWindow*);
  vtkOpenGLRenderWindow* GetContext();

  vtkTextureObject* GetTextureObject()
  {
    return this->TextureObject;
  };

  bool AllocateScalars(int dataType, int numComponents);
  bool AllocateScalarsFromPointer(int dataType, int numComponents, void* data);
  bool AllocateScalarsFromObject(int dataType, int numComponents, vtkDataObject* data);

  int GetScalarType();

  virtual void CopyInformationFromPipeline(vtkInformation* info) override;
  virtual void CopyInformationToPipeline(vtkInformation* info) override;

protected:
  vtkGPUImageData();
  ~vtkGPUImageData() override;

  vtkSmartPointer<vtkTextureObject> TextureObject;

  double Origin[3];
  double Spacing[3];
  int Extent[6];
  int Dimensions[3];



  vtkSmartPointer<vtkOpenGLRenderWindow> Context;

private:
  vtkGPUImageData(const vtkGPUImageData&) = delete;
  void operator=(const vtkGPUImageData&) = delete;
};

//----------------------------------------------------------------------------
inline double* vtkGPUImageData::GetPoint(vtkIdType id)
{
  vtkErrorMacro("todo");
  return 0;
}

//----------------------------------------------------------------------------
inline vtkIdType vtkGPUImageData::GetNumberOfPoints()
{
  const int* extent = this->Extent;
  const int* dims = this->GetDimensions();

  return dims[0] * dims[1] * dims[2];
}
#endif
