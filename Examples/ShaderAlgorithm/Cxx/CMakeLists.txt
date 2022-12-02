cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)
PROJECT (ShaderAlgorithm)

find_package(VTK 
    COMPONENTS
        CommonCore
        FiltersSources
        InteractionStyle
        RenderingCore
        RenderingOpenGL2
        InteractionImage
        IOImage
        ImagingColor
        vtksys
)
if (NOT VTK_FOUND)
  message("Skipping example: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif ()

add_executable(ShaderAlgorithm MACOSX_BUNDLE
    ShaderAlgorithm.cxx)
target_link_libraries(ShaderAlgorithm
    PRIVATE
        ${VTK_LIBRARIES})

set(RENDERING_EXAMPLES_SRCS
  ShaderAlgorithm
  )

vtk_module_autoinit(
  TARGETS ShaderAlgorithm
  MODULES ${VTK_LIBRARIES})