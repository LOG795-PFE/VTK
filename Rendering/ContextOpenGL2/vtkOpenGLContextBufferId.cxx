/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOpenGLContextBufferId.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOpenGLContextBufferId.h"

#include "vtkIntArray.h"
#include "vtkObjectFactory.h"
#include "vtkOpenGLError.h"
#include "vtkOpenGLRenderWindow.h"
#include "vtkOpenGLState.h"
#include "vtkOpenGLTexture.h"
#include "vtkTextureObject.h"
#include <cassert>

VTK_ABI_NAMESPACE_BEGIN
vtkStandardNewMacro(vtkOpenGLContextBufferId);

//------------------------------------------------------------------------------
vtkOpenGLContextBufferId::vtkOpenGLContextBufferId()
{
  this->Texture = nullptr;
  this->Context = nullptr;
}

//------------------------------------------------------------------------------
vtkOpenGLContextBufferId::~vtkOpenGLContextBufferId()
{
  if (this->Texture != nullptr)
  {
    vtkErrorMacro("texture should have been released.");
  }
}

//------------------------------------------------------------------------------
void vtkOpenGLContextBufferId::ReleaseGraphicsResources()
{
  if (this->Texture != nullptr)
  {
    this->Texture->Delete();
    this->Texture = nullptr;
  }
}

//------------------------------------------------------------------------------
void vtkOpenGLContextBufferId::SetContext(vtkRenderWindow* context)
{
  vtkOpenGLRenderWindow* c = vtkOpenGLRenderWindow::SafeDownCast(context);
  if (this->Context != c)
  {
    this->ReleaseGraphicsResources();
    this->Context = c;
    this->Modified();
  }
}

//------------------------------------------------------------------------------
vtkRenderWindow* vtkOpenGLContextBufferId::GetContext()
{
  return this->Context;
}

//------------------------------------------------------------------------------
bool vtkOpenGLContextBufferId::IsSupported()
{
  assert("pre: context_is_set" && this->GetContext() != nullptr);
  return vtkTextureObject::IsSupported(this->Context);
}

//------------------------------------------------------------------------------
void vtkOpenGLContextBufferId::Allocate()
{
  assert("pre: positive_width" && this->GetWidth() > 0);
  assert("pre: positive_height" && this->GetHeight() > 0);
  assert("pre: context_is_set" && this->GetContext() != nullptr);

  if (this->Texture == nullptr)
  {
    this->Texture = vtkTextureObject::New();
    this->Texture->SetContext(this->Context);
  }
  this->Context->MakeCurrent();
  // 3: RGB
  this->Texture->Allocate2D(static_cast<unsigned int>(this->GetWidth()),
    static_cast<unsigned int>(this->GetHeight()), 3, VTK_UNSIGNED_CHAR);
}

//------------------------------------------------------------------------------
bool vtkOpenGLContextBufferId::IsAllocated() const
{
  return this->Texture != nullptr &&
    this->Texture->GetWidth() == static_cast<unsigned int>(this->Width) &&
    this->Texture->GetHeight() == static_cast<unsigned int>(this->Height);
}

//------------------------------------------------------------------------------
void vtkOpenGLContextBufferId::SetValues(int srcXmin, int srcYmin)
{
  assert("pre: is_allocated" && this->IsAllocated());

  // copy the current read buffer to the texture.
  this->Texture->CopyFromFrameBuffer(srcXmin, srcYmin, 0, 0, this->Width, this->Height);
}

//------------------------------------------------------------------------------
vtkIdType vtkOpenGLContextBufferId::GetPickedItem(int x, int y)
{
  assert("pre: is_allocated" && this->IsAllocated());

  vtkOpenGLClearErrorMacro();

  vtkIdType result = -1;
  if (x < 0 || x >= this->Width)
  {
    vtkDebugMacro(<< "x mouse position out of range: x=" << x << " (width=" << this->Width << ")");
  }
  else
  {
    if (y < 0 || y >= this->Height)
    {
      vtkDebugMacro(<< "y mouse position out of range: y=" << y << " (height=" << this->Height
                    << ")");
    }
    else
    {
      this->Context->MakeCurrent();
      vtkOpenGLState* ostate = this->Context->GetState();

      // Render texture to current write buffer. Texel x,y is rendered at
      // pixel x,y (instead of pixel 0,0 to work around pixel ownership test).
      GLint savedDrawBuffer;
      glGetIntegerv(GL_DRAW_BUFFER, &savedDrawBuffer);

      vtkOpenGLState::ScopedglEnableDisable dsaver(ostate, GL_DEPTH_TEST);
      vtkOpenGLState::ScopedglEnableDisable ssaver(ostate, GL_STENCIL_TEST);
      vtkOpenGLState::ScopedglEnableDisable bsaver(ostate, GL_BLEND);

      if (savedDrawBuffer != GL_BACK_LEFT)
      {
        ostate->vtkglDrawBuffer(GL_BACK_LEFT);
      }
      ostate->vtkglDisable(GL_DEPTH_TEST);
      ostate->vtkglDisable(GL_STENCIL_TEST);
      ostate->vtkglDisable(GL_BLEND);

      this->Texture->CopyToFrameBuffer(x, y, x, y, x, y, this->Context->GetSize()[0],
        this->Context->GetSize()[1], nullptr, nullptr);

      GLint savedReadBuffer;
      glGetIntegerv(GL_READ_BUFFER, &savedReadBuffer);
      ostate->vtkglReadBuffer(GL_BACK_LEFT);

      // To workaround pixel ownership test,
      // get value from current read buffer at pixel (x,y) instead of just
      // (0,0).
      ostate->vtkglPixelStorei(GL_PACK_ALIGNMENT, 1);
      unsigned char rgb[3];
      rgb[0] = 5;
      rgb[1] = 1;
      rgb[2] = 8;
      glReadPixels(x, y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, rgb);

      if (savedReadBuffer != GL_BACK_LEFT)
      {
        ostate->vtkglReadBuffer(static_cast<GLenum>(savedReadBuffer));
      }
      if (savedDrawBuffer != GL_BACK_LEFT)
      {
        ostate->vtkglDrawBuffer(static_cast<GLenum>(savedDrawBuffer));
      }

      int value = (static_cast<int>(rgb[0]) << 16) | (static_cast<int>(rgb[1]) << 8) |
        static_cast<int>(rgb[2]);

      result = static_cast<vtkIdType>(value - 1);
    }
  }

  assert("post: valid_result" && result >= -1);

  vtkOpenGLCheckErrorMacro("failed after GetPickedItem");

  return result;
}

//------------------------------------------------------------------------------
void vtkOpenGLContextBufferId::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
VTK_ABI_NAMESPACE_END
