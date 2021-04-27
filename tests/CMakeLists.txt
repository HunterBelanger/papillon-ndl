cmake_minimum_required(VERSION 3.11)

project(PapillonNDL-Tests
  DESCRIPTION "Tests for PapillonNDL library"
  LANGUAGES CXX
)

find_package(GTest QUIET)
if(NOT GTest_FOUND)
  message(STATUS "Could not find a local install of GTest")
  message(STATUS "Will download GTest instead")
  
  # We don't want to force install GTest, so we turn off
  # this option.
  set(INSTALL_GTEST OFF CACHE BOOL "Enable installation of googletest.")
  # We also don't need googlemock so we turn that off as well
  set(BUILD_GMOCK OFF CACHE BOOL "Builds the googlemock subproject")

  include(FetchContent)

  FetchContent_Declare(GTest
    GIT_REPOSITORY https://github.com/google/googletest
    GIT_TAG        release-1.10.0
  )

  FetchContent_MakeAvailable(GTest)
else()
  message(STATUS "Using local install of GTest")
endif()

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
else() # Compile options for GCC and Clang on Unix-like systems
  target_compile_options(PapillonNDL-Tests PRIVATE -W -Wall -Wextra -Wpedantic)
endif()
