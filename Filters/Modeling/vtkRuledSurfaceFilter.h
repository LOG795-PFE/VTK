/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRuledSurfaceFilter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkRuledSurfaceFilter
 * @brief   generates a surface from a set of lines
 *
 * vtkRuledSurfaceFilter is a filter that generates a surface from a set of
 * lines. The lines are assumed to be "parallel" in the sense that they do
 * not intersect and remain somewhat close to one another. A surface is
 * generated by connecting the points defining each pair of lines with
 * straight lines. This creates a strip for each pair of lines (i.e., a
 * triangulation is created from two generating lines). The filter can handle
 * an arbitrary number of lines, with lines i and i+1 assumed connected.
 * Note that there are several different approaches for creating the ruled
 * surface, the method for creating the surface can either use the input
 * points or resample from the polylines (using a user-specified resolution).
 *
 * This filter offers some other important features. A DistanceFactor ivar is
 * used to decide when two lines are too far apart to connect. (The factor is
 * a multiple of the distance between the first two points of the two lines
 * defining the strip.) If the distance between the two generating lines
 * becomes too great, then the surface is not generated in that
 * region. (Note: if the lines separate and then merge, then a hole can be
 * generated in the surface.) In addition, the Offset and OnRation ivars can
 * be used to create nifty striped surfaces. Closed surfaces (e.g., tubes) can
 * be created by setting the CloseSurface ivar. (The surface can be closed
 * in the other direction by repeating the first and last point in the
 * polylines defining the surface.)
 *
 * An important use of this filter is to combine it with vtkStreamTracer to
 * generate stream surfaces. It can also be used to create surfaces from
 * contours.
 *
 * @warning
 * The number of lines must be greater than two if a surface is to be
 * generated.  sides (i.e., a ribbon), use vtkRibbonFilter.
 *
 * @sa
 * vtkRibbonFilter vtkStreamTracer
 */

#ifndef vtkRuledSurfaceFilter_h
#define vtkRuledSurfaceFilter_h

#include "vtkFiltersModelingModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

VTK_ABI_NAMESPACE_BEGIN
class vtkIdList;
class vtkPoints;
class vtkPolyData;

#define VTK_RULED_MODE_RESAMPLE 0
#define VTK_RULED_MODE_POINT_WALK 1

class VTKFILTERSMODELING_EXPORT vtkRuledSurfaceFilter : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkRuledSurfaceFilter, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /**
   * Construct object with OnRatio=1, Offset=0. DistanceFactor=3.0,
   * CloseSurface off, and PassLines off.
   */
  static vtkRuledSurfaceFilter* New();

  ///@{
  /**
   * Set/Get the factor that controls tearing of the surface.
   */
  vtkSetClampMacro(DistanceFactor, double, 1.0, VTK_DOUBLE_MAX);
  vtkGetMacro(DistanceFactor, double);
  ///@}

  ///@{
  /**
   * Control the striping of the ruled surface. If OnRatio is greater
   * than 1, then every nth strip is turned on, beginning with the Offset
   * strip.
   */
  vtkSetClampMacro(OnRatio, int, 1, VTK_INT_MAX);
  vtkGetMacro(OnRatio, int);
  ///@}

  ///@{
  /**
   * Control the striping of the ruled surface. The offset sets the
   * first stripe that is visible. Offset is generally used with
   * OnRatio to create nifty striping effects.
   */
  vtkSetClampMacro(Offset, int, 0, VTK_INT_MAX);
  vtkGetMacro(Offset, int);
  ///@}

  ///@{
  /**
   * Indicate whether the surface is to be closed. If this boolean is
   * on, then the first and last polyline are used to generate a stripe
   * that closes the surface. (Note: to close the surface in the other
   * direction, repeat the first point in the polyline as the last
   * point in the polyline.)
   */
  vtkSetMacro(CloseSurface, vtkTypeBool);
  vtkGetMacro(CloseSurface, vtkTypeBool);
  vtkBooleanMacro(CloseSurface, vtkTypeBool);
  ///@}

  ///@{
  /**
   * Set the mode by which to create the ruled surface. (Dramatically
   * different results are possible depending on the chosen mode.) The
   * resample mode evenly resamples the polylines (based on length) and
   * generates triangle strips. The point walk mode uses the existing
   * points and walks around the polyline using existing points.
   */
  vtkSetClampMacro(RuledMode, int, VTK_RULED_MODE_RESAMPLE, VTK_RULED_MODE_POINT_WALK);
  vtkGetMacro(RuledMode, int);
  void SetRuledModeToResample() { this->SetRuledMode(VTK_RULED_MODE_RESAMPLE); }
  void SetRuledModeToPointWalk() { this->SetRuledMode(VTK_RULED_MODE_POINT_WALK); }
  const char* GetRuledModeAsString();
  ///@}

  ///@{
  /**
   * If the ruled surface generation mode is RESAMPLE, then these parameters
   * are used to determine the resample rate. Resolution[0] defines the
   * resolution in the direction of the polylines; Resolution[1] defines
   * the resolution across the polylines (i.e., direction orthogonal to
   * Resolution[0]).
   */
  vtkSetVector2Macro(Resolution, int);
  vtkGetVectorMacro(Resolution, int, 2);
  ///@}

  ///@{
  /**
   * Indicate whether the generating lines are to be passed to the output.
   * By default lines are not passed to the output.
   */
  vtkSetMacro(PassLines, vtkTypeBool);
  vtkGetMacro(PassLines, vtkTypeBool);
  vtkBooleanMacro(PassLines, vtkTypeBool);
  ///@}

  ///@{
  /**
   * Indicate whether the starting points of the loops need to be determined.
   * If set to 0, then its assumes that the 0th point of each loop should be
   * always connected
   * By default the loops are not oriented.
   */
  vtkSetMacro(OrientLoops, vtkTypeBool);
  vtkGetMacro(OrientLoops, vtkTypeBool);
  vtkBooleanMacro(OrientLoops, vtkTypeBool);
  ///@}

protected:
  vtkRuledSurfaceFilter();
  ~vtkRuledSurfaceFilter() override;

  // Usual data generation method
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  double DistanceFactor;
  int OnRatio;
  int Offset;
  vtkTypeBool CloseSurface;
  int RuledMode;
  int Resolution[2];
  vtkTypeBool PassLines;
  vtkTypeBool OrientLoops;

private:
  vtkIdList* Ids;
  double Weights[4];

  void Resample(vtkPolyData* output, vtkPolyData* input, vtkPoints* inPts, vtkPoints* newPts,
    int npts, const vtkIdType* pts, int npts2, const vtkIdType* pts2);
  void PointWalk(vtkPolyData* output, vtkPoints* inPts, int npts, const vtkIdType* pts, int npts2,
    const vtkIdType* pts2);

private:
  vtkRuledSurfaceFilter(const vtkRuledSurfaceFilter&) = delete;
  void operator=(const vtkRuledSurfaceFilter&) = delete;
};

VTK_ABI_NAMESPACE_END
#endif
