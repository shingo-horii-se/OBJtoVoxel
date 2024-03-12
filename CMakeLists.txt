
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(OBJtoVoxel)

find_package(VTK COMPONENTS 
  # CommonColor
  CommonCore
  CommonDataModel
  FiltersCore
  FiltersGeneral
  FiltersModeling
  FiltersSources
  IOGeometry
  IOXML
  IOImage
  ImagingCore
  ImagingHybrid
  ImagingStencil
  # InteractionStyle
  # RenderingContextOpenGL2
  # RenderingCore
  # RenderingFreeType
  # RenderingGL2PSOpenGL2
  # RenderingOpenGL2
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "OBJtoVoxel: Unable to find the VTK build folder.")
endif()

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")
add_executable(OBJtoVoxel MACOSX_BUNDLE OBJtoVoxel.cxx )
  target_link_libraries(OBJtoVoxel PRIVATE ${VTK_LIBRARIES}
)
# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS OBJtoVoxel
  MODULES ${VTK_LIBRARIES}
)
