/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkParticleTracerBase.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkParticleTracerBase
 * @brief   A parallel particle tracer for vector fields
 *
 * vtkPParticleTracerBase is the base class for parallel filters that advect particles
 * in a vector field. Note that the input vtkPointData structure must
 * be identical on all datasets.
 * @sa
 * vtkRibbonFilter vtkRuledSurfaceFilter vtkInitialValueProblemSolver
 * vtkRungeKutta2 vtkRungeKutta4 vtkRungeKutta45 vtkStreamTracer
 */

#ifndef vtkPParticleTracerBase_h
#define vtkPParticleTracerBase_h

#include "vtkParticleTracerBase.h"
#include "vtkSmartPointer.h" // For protected ivars.

#include <vector> // STL Header

#include "vtkFiltersParallelFlowPathsModule.h" // For export macro

VTK_ABI_NAMESPACE_BEGIN
class VTKFILTERSPARALLELFLOWPATHS_EXPORT vtkPParticleTracerBase : public vtkParticleTracerBase
{
public:
  vtkTypeMacro(vtkPParticleTracerBase, vtkParticleTracerBase);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  ///@{
  /**
   * Set/Get the controller used when sending particles between processes
   * The controller must be an instance of vtkMPIController.
   */
  virtual void SetController(vtkMultiProcessController* controller);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  ///@}

protected:
  struct RemoteParticleInfo
  {
    vtkParticleTracerBaseNamespace::ParticleInformation Current;
    vtkParticleTracerBaseNamespace::ParticleInformation Previous;
    vtkSmartPointer<vtkPointData> PreviousPD;
  };

  typedef std::vector<RemoteParticleInfo> RemoteParticleVector;

  vtkPParticleTracerBase();
  ~vtkPParticleTracerBase() override;

  virtual int RequestUpdateExtent(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

  //

  vtkPolyData* Execute(vtkInformationVector** inputVector) override;
  bool SendParticleToAnotherProcess(vtkParticleTracerBaseNamespace::ParticleInformation& info,
    vtkParticleTracerBaseNamespace::ParticleInformation& previous, vtkPointData*) override;

  /**
   * Before starting the particle trace, classify
   * all the injection/seed points according to which processor
   * they belong to. This saves us retesting at every injection time
   * providing 1) The volumes are static, 2) the seed points are static
   * If either are non static, then this step is skipped.
   */
  void AssignSeedsToProcessors(double time, vtkDataSet* source, int sourceID, int ptId,
    vtkParticleTracerBaseNamespace::ParticleVector& localSeedPoints,
    int& localAssignedCount) override;

  /**
   * give each one a unique ID. We need to use MPI to find out
   * who is using which numbers.
   */
  void AssignUniqueIds(vtkParticleTracerBaseNamespace::ParticleVector& localSeedPoints) override;

  /**
   * this is used during classification of seed points and also between iterations
   * of the main loop as particles leave each processor domain. Returns
   * true if particles were migrated to any new process.
   */
  bool SendReceiveParticles(RemoteParticleVector& outofdomain, RemoteParticleVector& received);

  bool UpdateParticleListFromOtherProcesses() override;

  /**
   * Method that checks that the input arrays are ordered the
   * same on all data sets. This needs to be true for all
   * blocks in a composite data set as well as across all processes.
   */
  bool IsPointDataValid(vtkDataObject* input) override;

  //

  //

  // MPI controller needed when running in parallel
  vtkMultiProcessController* Controller;

  // List used for transmitting between processors during parallel operation
  RemoteParticleVector MPISendList;

  RemoteParticleVector Tail; // this is to receive the "tails" of traces from other processes
private:
  vtkPParticleTracerBase(const vtkPParticleTracerBase&) = delete;
  void operator=(const vtkPParticleTracerBase&) = delete;
};
VTK_ABI_NAMESPACE_END
#endif
