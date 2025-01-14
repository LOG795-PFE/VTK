/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTkWidgetsInit.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkTcl.h"
#include "vtkTk.h"

#include "vtkImageData.h"
#include "vtkTkImageViewerWidget.h"
#include "vtkTkRenderWidget.h"
#include "vtkVersionMacros.h"

//------------------------------------------------------------------------------
// Vtkrenderingtk_Init
// Called upon system startup to create the widget commands.
VTK_ABI_NAMESPACE_BEGIN
extern "C"
{
  VTK_EXPORT int Vtkrenderingtk_Init(Tcl_Interp* interp);
}

extern "C"
{
  VTK_EXPORT int Vtktkrenderwidget_Init(Tcl_Interp* interp);
}
extern "C"
{
  VTK_EXPORT int Vtktkimageviewerwidget_Init(Tcl_Interp* interp);
}

#define VTKTK_TO_STRING(x) VTKTK_TO_STRING0(x)
#define VTKTK_TO_STRING0(x) VTKTK_TO_STRING1(x)
#define VTKTK_TO_STRING1(x) #x
#define VTKTK_VERSION VTKTK_TO_STRING(VTK_MAJOR_VERSION) "." VTKTK_TO_STRING(VTK_MINOR_VERSION)

int Vtkrenderingtk_Init(Tcl_Interp* interp)
{
  // Forward the call to the real init functions.
  if (Vtktkrenderwidget_Init(interp) == TCL_OK && Vtktkimageviewerwidget_Init(interp) == TCL_OK)
  {
    // Report that the package is provided.
    return Tcl_PkgProvide(interp, (char*)"Vtkrenderingtk", (char*)VTKTK_VERSION);
  }
  else
  {
    // One of the widgets is not provided.
    return TCL_ERROR;
  }
}
VTK_ABI_NAMESPACE_END
