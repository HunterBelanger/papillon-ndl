cmake_minimum_required(VERSION 3.9)
project(PapillonNDL
  #VERSION 0.0.1
  DESCRIPTION "Library access ACE nuclear data"
  LANGUAGES CXX
)

# Add options
option(PNDL_GO_FASTER "Enable -Ofast -ffast-math compiler optimizations" OFF)
option(PNDL_PYTHON "Enable Python interface to PapillonNDL" OFF)
option(PNDL_TESTS "Build tests" OFF)

# Create PapillonNDL library target
add_library(PapillonNDL SHARED
  src/frame.cpp
  src/interpolation.cpp
  src/integration.cpp
  src/region_1d.cpp
  src/multi_region_1d.cpp
  src/polynomial_1d.cpp
  src/ace.cpp
  src/isotropic.cpp
  src/equiprobable_angle_bins.cpp
  src/angle_table.cpp
  src/equiprobable_energy_bins.cpp
  src/level_inelastic_scatter.cpp
  src/general_evaporation.cpp
  src/evaporation.cpp
  src/maxwellian.cpp
  src/watt.cpp
  src/tabular_energy.cpp
  src/pctable.cpp
  src/angle_distribution.cpp
  src/uncorrelated.cpp
  src/nbody.cpp
  src/kalbach_table.cpp
  src/kalbach.cpp
  src/energy_angle_table.cpp
  src/tabular_energy_angle.cpp
  src/energy_grid.cpp
  src/cross_section.cpp
  src/reaction.cpp
  src/delayed_group.cpp
  src/fission_data.cpp
  src/nuclide.cpp
)

# Add public interface includes directory
target_include_directories(PapillonNDL PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Require C++17 standard
target_compile_features(PapillonNDL PRIVATE cxx_std_17)

# Get compiler info
# TODO, for now, only support Clang or GNU

# Set compile options for library
target_compile_options(PapillonNDL PRIVATE -W -Wall -Wextra -Wpedantic -Weffc++)
target_compile_options(PapillonNDL PRIVATE $<$<CONFIG:DEBUG>:-g>)
target_compile_options(PapillonNDL PRIVATE $<$<CONFIG:RELEASE>:-O2>)
target_compile_options(PapillonNDL PRIVATE $<$<BOOL:PNDL_GO_FASTER>:-Ofast -ffast-math>)

if(PNDL_TESTS)
  add_subdirectory(tests)
endif()

if(PNDL_PYTHON)
  find_package(pybind11 REQUIRED)

  pybind11_add_module(pyPapillonNDL src/python/pyPapillonNDL.cpp
                                    src/python/interpolation.cpp
                                    src/python/ace.cpp
                                    src/python/function_1d.cpp
                                    src/python/frame.cpp
                                    src/python/angle_law.cpp
                                    src/python/angle_distribution.cpp
                                    src/python/energy_law.cpp
                                    src/python/angle_energy.cpp
                                    src/python/pctable.cpp
                                    src/python/energy_grid.cpp
                                    src/python/cross_section.cpp
                                    src/python/reaction.cpp
                                    src/python/delayed_group.cpp
                                    src/python/fission_data.cpp
                                    src/python/nuclide.cpp
                                    src/python/prng.cpp)

  target_link_libraries(pyPapillonNDL PRIVATE PapillonNDL)

  # TODO install pyPapillonNDL

endif()

# Install
install(TARGETS PapillonNDL EXPORT PapillonNDL DESTINATION lib)
install(EXPORT PapillonNDL FILE PapillonNDL-config.cmake DESTINATION share/cmake/PapillonNDL)
install(DIRECTORY include/PapillonNDL DESTINATION include)
