cmake_minimum_required(VERSION 3.11)
project(panglos LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(interpolation
    GIT_REPOSITORY https://github.com/njoy/interpolation
    GIT_TAG        develop
)

FetchContent_Declare(ENDFtk
    GIT_REPOSITORY https://github.com/njoy/ENDFtk
    GIT_TAG        v0.3.0
)

FetchContent_Declare(Log
    GIT_REPOSITORY  https://github.com/njoy/Log
    GIT_TAG         52962b7796afe37ef1d8f7edb4bf9ecb1b868d15
)

FetchContent_MakeAvailable(
    interpolation
    ENDFtk
    Log
)

FetchContent_Declare(NDArray
  GIT_REPOSITORY https://github.com/HunterBelanger/ndarray.git
  GIT_TAG        master
)
set(NDARRAY_INSTALL OFF CACHE BOOL "Install NDArray")
FetchContent_MakeAvailable(NDArray)

add_executable(panglos src/main.cpp
                       src/short_collision_time_sab.cpp
                       src/free_gas_sab.cpp
                       src/tabulated_sab.cpp
                       src/incoherent_inelastic.cpp
                       src/incoherent_elastic.cpp
                       src/coherent_elastic.cpp
                       src/gauss_kronrod.cpp
)
target_compile_features(panglos PRIVATE cxx_std_17)
target_compile_options(panglos PRIVATE -W -Wall -Wextra -Wpedantic)
target_link_libraries(panglos PRIVATE ENDFtk Log interpolation NDArray)

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  target_link_libraries(panglos PRIVATE OpenMP::OpenMP_CXX)
endif()
