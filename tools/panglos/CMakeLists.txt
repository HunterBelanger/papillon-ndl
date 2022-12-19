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

FetchContent_Declare(disco
  GIT_REPOSITORY    https://github.com/njoy/disco
  GIT_TAG           2606933
)

FetchContent_Declare(docopt.cpp
  GIT_REPOSITORY https://github.com/docopt/docopt.cpp.git
  GIT_TAG        v0.6.3
)

FetchContent_MakeAvailable(
    interpolation
    ENDFtk
    Log
    disco
    docopt.cpp
)

add_executable(panglos src/main.cpp
                       src/short_collision_time_sab.cpp
                       src/free_gas_sab.cpp
                       src/tabulated_sab.cpp
                       src/incoherent_inelastic.cpp
                       src/incoherent_elastic.cpp
                       src/coherent_elastic.cpp
                       src/gauss_kronrod.cpp
                       src/ace.cpp
)
target_compile_features(panglos PRIVATE cxx_std_17)
target_compile_options(panglos PRIVATE -W -Wall -Wextra -Wpedantic)
target_link_libraries(panglos PRIVATE ENDFtk Log interpolation disco docopt_s)

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  target_link_libraries(panglos PRIVATE OpenMP::OpenMP_CXX)
endif()
