/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOSPRayViewNodeFactory.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOSPRayViewNodeFactory.h"
#include "vtkObjectFactory.h"

#include "vtkOSPRayAMRVolumeMapperNode.h"
#include "vtkOSPRayActorNode.h"
#include "vtkOSPRayCameraNode.h"
#include "vtkOSPRayCompositePolyDataMapper2Node.h"
#include "vtkOSPRayLightNode.h"
#include "vtkOSPRayMoleculeMapperNode.h"
#include "vtkOSPRayPointGaussianMapperNode.h"
#include "vtkOSPRayPolyDataMapperNode.h"
#include "vtkOSPRayRendererNode.h"
#include "vtkOSPRayUnstructuredVolumeMapperNode.h"
#include "vtkOSPRayVolumeMapperNode.h"
#include "vtkOSPRayVolumeNode.h"

VTK_ABI_NAMESPACE_BEGIN
vtkViewNode* ren_maker()
{
  vtkOSPRayRendererNode* vn = vtkOSPRayRendererNode::New();
  return vn;
}

vtkViewNode* amrm_maker()
{
  vtkOSPRayAMRVolumeMapperNode* vn = vtkOSPRayAMRVolumeMapperNode::New();
  return vn;
}

vtkViewNode* act_maker()
{
  vtkOSPRayActorNode* vn = vtkOSPRayActorNode::New();
  return vn;
}

vtkViewNode* vol_maker()
{
  return vtkOSPRayVolumeNode::New();
}

vtkViewNode* cam_maker()
{
  vtkOSPRayCameraNode* vn = vtkOSPRayCameraNode::New();
  return vn;
}

vtkViewNode* light_maker()
{
  vtkOSPRayLightNode* vn = vtkOSPRayLightNode::New();
  return vn;
}

vtkViewNode* pd_maker()
{
  vtkOSPRayPolyDataMapperNode* vn = vtkOSPRayPolyDataMapperNode::New();
  return vn;
}

vtkViewNode* molecule_maker()
{
  vtkOSPRayMoleculeMapperNode* vn = vtkOSPRayMoleculeMapperNode::New();
  return vn;
}

vtkViewNode* vm_maker()
{
  vtkOSPRayVolumeMapperNode* vn = vtkOSPRayVolumeMapperNode::New();
  return vn;
}

vtkViewNode* cpd_maker()
{
  vtkOSPRayCompositePolyDataMapper2Node* vn = vtkOSPRayCompositePolyDataMapper2Node::New();
  return vn;
}

vtkViewNode* tetm_maker()
{
  vtkOSPRayUnstructuredVolumeMapperNode* vn = vtkOSPRayUnstructuredVolumeMapperNode::New();
  return vn;
}

vtkViewNode* particle_maker()
{
  vtkOSPRayPointGaussianMapperNode* vn = vtkOSPRayPointGaussianMapperNode::New();
  return vn;
}

//============================================================================
vtkStandardNewMacro(vtkOSPRayViewNodeFactory);

//------------------------------------------------------------------------------
vtkOSPRayViewNodeFactory::vtkOSPRayViewNodeFactory()
{
  // see vtkRenderWindow::GetRenderLibrary
  this->RegisterOverride("vtkOpenGLRenderer", ren_maker);
  this->RegisterOverride("vtkOpenGLActor", act_maker);
  this->RegisterOverride("vtkPVLODActor", act_maker);
  this->RegisterOverride("vtkPVLODVolume", vol_maker);
  this->RegisterOverride("vtkVolume", vol_maker);
  this->RegisterOverride("vtkOpenGLCamera", cam_maker);
  this->RegisterOverride("vtkPVCamera", cam_maker);
  this->RegisterOverride("vtkOpenGLLight", light_maker);
  this->RegisterOverride("vtkPVLight", light_maker);
  this->RegisterOverride("vtkPainterPolyDataMapper", pd_maker);
  this->RegisterOverride("vtkOpenGLPolyDataMapper", pd_maker);
  this->RegisterOverride("vtkSmartVolumeMapper", vm_maker);
  this->RegisterOverride("vtkOSPRayVolumeMapper", vm_maker);
  this->RegisterOverride("vtkOpenGLGPUVolumeRayCastMapper", vm_maker);
  this->RegisterOverride("vtkMultiBlockVolumeMapper", vm_maker);
  this->RegisterOverride("vtkCompositePolyDataMapper2", cpd_maker);
  this->RegisterOverride("vtkOpenGLProjectedTetrahedraMapper", tetm_maker);
  this->RegisterOverride("vtkUnstructuredGridVolumeZSweepMapper", tetm_maker);
  this->RegisterOverride("vtkUnstructuredGridVolumeRayCastMapper", tetm_maker);
  this->RegisterOverride("vtkAMRVolumeMapper", amrm_maker);
  this->RegisterOverride("vtkMoleculeMapper", molecule_maker);
  this->RegisterOverride("vtkOpenGLPointGaussianMapper", particle_maker);
}

//------------------------------------------------------------------------------
vtkOSPRayViewNodeFactory::~vtkOSPRayViewNodeFactory() = default;

//------------------------------------------------------------------------------
void vtkOSPRayViewNodeFactory::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
VTK_ABI_NAMESPACE_END
