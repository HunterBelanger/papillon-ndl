cmake_minimum_required(VERSION 3.11)

project(PapillonNDLTests
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

  if(WIN32)
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  endif()

  include(FetchContent)

  FetchContent_Declare(GTest
    GIT_REPOSITORY https://github.com/google/googletest
    GIT_TAG        release-1.10.0
  )

  FetchContent_MakeAvailable(GTest)
else()
  message(STATUS "Using local install of GTest")
endif()

# Interpolation Tests
add_executable(InterpolationTests interpolation.cpp)
target_compile_features(InterpolationTests PRIVATE cxx_std_17)
target_link_libraries(InterpolationTests PUBLIC PapillonNDL gtest_main)
add_test(InterpolationTests InterpolationTests)

# Function1D Tests
add_executable(Function1DTests function_1d.cpp)
target_compile_features(Function1DTests PRIVATE cxx_std_17)
target_link_libraries(Function1DTests PUBLIC PapillonNDL gtest_main)
add_test(Function1DTests Function1DTests)

# PCTable Tests
add_executable(PCTableTests pctable.cpp)
target_compile_features(PCTableTests PRIVATE cxx_std_17)
target_link_libraries(PCTableTests PUBLIC PapillonNDL gtest_main)
add_test(PCTableTests PCTableTests)

# AngleLaw Tests
add_executable(AngleLawTests angle_law.cpp)
target_compile_features(AngleLawTests PRIVATE cxx_std_17)
target_link_libraries(AngleLawTests PUBLIC PapillonNDL gtest_main)
add_test(AngleLawTests AngleLawTests)
