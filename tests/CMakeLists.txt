cmake_minimum_required(VERSION 3.9)

project(PapillonNDL-Tests
  DESCRIPTION "Tests for PapillonNDL library"
  LANGUAGES CXX
)

find_package(GTest REQUIRED)
# TODO Use FetchContent if GTest isn't found

# Create PapillonNDL library target
add_executable(PapillonNDL-Tests
  main.cpp
  interpolation_tests.cpp
  integration_tests.cpp
  polynomial_1d_tests.cpp
  region_1d_tests.cpp
  multi_region_1d_tests.cpp
)

target_compile_features(PapillonNDL-Tests PRIVATE cxx_std_17)

# Need to "link" agains cereal (this includes it's header files)
target_link_libraries(PapillonNDL-Tests PUBLIC PapillonNDL gtest)

if(WIN32) # Comile options for Windows
  target_compile_options(PapillonNDL-Tests PRIVATE /W4)
  target_compile_options(PapillonNDL-Tests PRIVATE $<$<CONFIG:DEBUG>:/DEBUG>)
  target_compile_options(PapillonNDL-Tests PRIVATE $<$<CONFIG:RELEASE>:/O2>)
  target_compile_options(PapillonNDL-Tests PRIVATE $<$<BOOL:PNDL_GO_FASTER>:/O2>)
else() # Compile options for GCC and Clang on Unix-like systems
  target_compile_options(PapillonNDL-Tests PRIVATE -W -Wall -Wextra -Wpedantic -Weffc++)
  target_compile_options(PapillonNDL-Tests PRIVATE $<$<CONFIG:DEBUG>:-g>)
  target_compile_options(PapillonNDL-Tests PRIVATE $<$<CONFIG:RELEASE>:-O2>)
  target_compile_options(PapillonNDL-Tests PRIVATE $<$<BOOL:PNDL_GO_FASTER>:-Ofast -ffast-math>)
endif()