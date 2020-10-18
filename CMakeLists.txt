cmake_minimum_required(VERSION 3.9)
project(PapillonNDL
  #VERSION 0.0.1
  DESCRIPTION "Library access ACE nuclear data"
  LANGUAGES CXX
)

# Add options
option(PNDL_GO_FASTER "Enable -Ofast -ffast-math compiler optimizations" OFF)
option(PNDL_TESTS "Build tests" OFF)

# Create PapillonNDL library target
add_library(PapillonNDL SHARED
  src/interpolation.cpp
  src/integration.cpp
  src/region_1d.cpp
  src/multi_region_1d.cpp
  src/polynomial_1d.cpp
  src/ace.cpp
  src/equiprobable_energy_bins.cpp
  src/level_inelastic_scatter.cpp
  src/general_evaporation.cpp
  src/evaporation.cpp
  src/maxwellian.cpp
  src/watt.cpp
  src/tabular_energy.cpp
  src/pctable.cpp
  src/angle_table.cpp
  src/equiprobable_angle_bins.cpp
  src/isotropic.cpp
  src/angle_distribution.cpp
  src/uncorrelated.cpp
  src/nbody.cpp
  src/kalbach_table.cpp
  src/kalbach.cpp
  src/energy_angle_table.cpp
  src/tabular_energy_angle.cpp
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
target_compile_options(PapillonNDL PRIVATE $<$<CONFIG:DEBUG>:-g -fsanitizer=address,undefined>)
target_compile_options(PapillonNDL PRIVATE $<$<CONFIG:RELEASE>:-O2>)
target_compile_options(PapillonNDL PRIVATE $<$<BOOL:PNDL_GO_FASTER>:-Ofast -ffast-math>)

if(PNDL_TESTS)
  add_subdirectory(tests)
endif()

# Install
install(TARGETS PapillonNDL EXPORT PapillonNDL DESTINATION lib)
install(EXPORT PapillonNDL FILE PapillonNDL-config.cmake DESTINATION share/cmake/PapillonNDL)
install(DIRECTORY include/PapillonNDL DESTINATION include)
