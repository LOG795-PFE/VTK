/*=========================================================================

Program:   Visualization Toolkit
Module:    vtkMemoryResourceStream.cxx

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkMemoryResourceStream.h"

#include "vtkObjectFactory.h"

#include <cstring>

VTK_ABI_NAMESPACE_BEGIN

vtkStandardNewMacro(vtkMemoryResourceStream);

//------------------------------------------------------------------------------
vtkMemoryResourceStream::vtkMemoryResourceStream()
  : vtkResourceStream{ true }
{
}

//------------------------------------------------------------------------------
void vtkMemoryResourceStream::SetBuffer(const void* buffer, std::size_t size)
{
  if (buffer == nullptr && size != 0)
  {
    vtkErrorMacro("buffer must not be nullptr if size > 0");
    return;
  }

  this->Buffer = static_cast<const unsigned char*>(buffer);
  this->Size = size;
  this->Pos = 0;
  this->Eos = (this->Size == 0);

  this->Modified();
}

//------------------------------------------------------------------------------
std::size_t vtkMemoryResourceStream::Read(void* buffer, std::size_t bytes)
{
  if (bytes == 0)
  {
    return 0;
  }

  const auto sbytes = static_cast<vtkTypeInt64>(bytes);
  const auto ssize = static_cast<vtkTypeInt64>(this->Size);
  const auto read = std::min(sbytes, ssize - this->Pos);

  if (read <= 0)
  {
    this->Eos = true;
    return 0;
  }

  std::memcpy(buffer, this->Buffer + this->Pos, static_cast<std::size_t>(read));
  this->Pos += read;
  this->Eos = read != sbytes;

  return read;
}

//------------------------------------------------------------------------------
bool vtkMemoryResourceStream::EndOfStream()
{
  return this->Eos;
}

//------------------------------------------------------------------------------
vtkTypeInt64 vtkMemoryResourceStream::Seek(vtkTypeInt64 pos, SeekDirection dir)
{
  if (dir == SeekDirection::Begin)
  {
    this->Pos = pos;
  }
  else if (dir == SeekDirection::Current)
  {
    this->Pos += pos;
  }
  else
  {
    this->Pos = static_cast<vtkTypeInt64>(this->Size) + pos;
  }

  this->Eos = false;
  return this->Pos;
}

//------------------------------------------------------------------------------
vtkTypeInt64 vtkMemoryResourceStream::Tell()
{
  return this->Pos;
}

//------------------------------------------------------------------------------
void vtkMemoryResourceStream::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Buffer: " << this->Buffer << "\n";
  os << indent << "Size: " << this->Size << "o\n";
  os << indent << "Position: " << this->Pos << "\n";
}

VTK_ABI_NAMESPACE_END
