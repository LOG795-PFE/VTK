cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)


PROJECT (ShaderAlgorithm)

if(NOT VTK_SOURCE_DIR)
  message(ERROR " Cannot build ShaderAlgorithm examples without VTK_SOURCE_DIR")
endif()
find_package(VTK REQUIRED)
vtk_module_config(VTK
  vtkCommonCore
  vtkFiltersSources
  vtkInteractionStyle
  vtkRenderingOpenGL2
  vtkInteractionImage
  vtkIOImage
  vtkImagingColor
)
include(${VTK_USE_FILE})

set(RENDERING_EXAMPLES_SRCS
  ShaderAlgorithm
  )

foreach(name ${RENDERING_EXAMPLES_SRCS})
  add_executable(${name} MACOSX_BUNDLE ${name}.cxx)
  target_link_libraries(${name} ${VTK_LIBRARIES} )
endforeach()
