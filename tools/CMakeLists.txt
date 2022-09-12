cmake_minimum_required(VERSION 3.9)

# Get FetchContent for downloading dependencies
include(FetchContent)

project(spndl LANGUAGES CXX)

#===============================================================================
# Get docopt
message(STATUS "Downloading docopt.cpp v0.6.3")
FetchContent_Declare(docopt.cpp
  GIT_REPOSITORY https://github.com/docopt/docopt.cpp
  GIT_TAG        v0.6.3
)
if(NOT docopt.cpp_POPULATED)
  FetchContent_Populate(docopt.cpp)
  add_subdirectory(${docopt.cpp_SOURCE_DIR} ${docopt.cpp_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

#===============================================================================
# Get NDArray
message(STATUS "Downloading NDArray")
FetchContent_Declare(NDArray
  GIT_REPOSITORY https://github.com/HunterBelanger/ndarray
  GIT_TAG        master
)
set(NDARRAY_INSTALL OFF CACHE BOOL "Install NDArray")
FetchContent_MakeAvailable(NDArray)

#===============================================================================
# Get OpenMC
message(STATUS "Downloading OpenMC v0.13.1")
include(FetchContent)
FetchContent_Declare(OpenMC
  GIT_REPOSITORY https://github.com/openmc-dev/openmc
  GIT_TAG        v0.13.1
)
FetchContent_MakeAvailable(OpenMC)

# Build PapillonNDL sampler
add_executable(spndl src/spndl.cpp)
target_compile_features(spndl PRIVATE cxx_std_17)
target_compile_options(spndl PRIVATE -W -Wall -Wextra -Wpedantic)
target_compile_options(spndl PRIVATE $<$<CONFIG:DEBUG>:-g>)
target_link_libraries(spndl PUBLIC PapillonNDL::PapillonNDL docopt NDArray)

# Build OpenMC sampler
add_executable(sopenmc src/sopenmc.cpp)
target_compile_features(sopenmc PRIVATE cxx_std_17)
target_compile_options(sopenmc PRIVATE -W -Wall -Wextra -Wpedantic)
target_compile_options(sopenmc PRIVATE $<$<CONFIG:DEBUG>:-g>)
target_link_libraries(sopenmc PUBLIC libopenmc docopt NDArray)

# Copy python comparison file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compare.py ${CMAKE_CURRENT_BINARY_DIR}/compare.py COPYONLY)
