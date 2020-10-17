cmake_minimum_required(VERSION 3.9)
project(PapillonNDL-Tests
  DESCRIPTION "Tests for PapillonNDL library"
  LANGUAGES CXX
)

find_package(GTest REQUIRED)

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

# Set compile options for library
target_compile_options(PapillonNDL-Tests PRIVATE -W -Wall -Wextra -Wpedantic -Weffc++)
target_compile_options(PapillonNDL-Tests PRIVATE $<$<CONFIG:DEBUG>:-g -fsanitizer=address,undefined>)
target_compile_options(PapillonNDL-Tests PRIVATE $<$<CONFIG:RELEASE>:-O2>)
