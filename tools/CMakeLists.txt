cmake_minimum_required(VERSION 3.11)

option(PNDL_SAMPLERS "Build PapillonNDL and OpenMC samplers" OFF)
option(PNDL_PANGLOS "Build Panglos TSL processing code" OFF)

if(PNDL_SAMPLERS)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/samplers)
endif()

if(PNDL_PANGLOS)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/panglos)
endif()
