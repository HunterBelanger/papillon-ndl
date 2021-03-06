cmake_minimum_required(VERSION 3.11)

project(PapillonNDL
  VERSION "0.2.1"
  DESCRIPTION "Library to interact with ACE nuclear data files."
  LANGUAGES CXX
)

include(CMakePackageConfigHelpers)

# Add options
option(PNDL_PYTHON "Enable Python interface to PapillonNDL" ON)
option(PNDL_TESTS "Build tests" OFF)

# List of source files for PapillonNDL
set(PNDL_SOURCE_LIST src/region_1d.cpp
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
                     src/multiple_distribution.cpp
                     src/cm_distribution.cpp
                     src/energy_grid.cpp
                     src/cross_section.cpp
                     src/reaction.cpp
                     src/delayed_group.cpp
                     src/ce_neutron.cpp
                     src/continuous_energy_discrete_cosines.cpp
                     src/discrete_cosines_energies.cpp
                     src/st_incoherent_inelastic.cpp
                     src/st_coherent_elastic.cpp
                     src/st_incoherent_elastic.cpp
                     src/st_thermal_scattering_law.cpp
                     src/version.cpp
)

# Create PapillonNDL library target
if(WIN32)
  # On windows, make a static lib to avoid errors using cmake
  add_library(PapillonNDL STATIC ${PNDL_SOURCE_LIST})
else()
  # On Unix-like systems, use a shared lib
  add_library(PapillonNDL SHARED ${PNDL_SOURCE_LIST})
endif()

# Add alias to make more friendly with FetchConent
add_library(PapillonNDL::PapillonNDL ALIAS PapillonNDL)

# Add public interface includes directory
target_include_directories(PapillonNDL PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Require C++17 standard
target_compile_features(PapillonNDL PRIVATE cxx_std_17)

if(WIN32) # Comile options for Windows
  target_compile_options(PapillonNDL PRIVATE /W4)
else() # Compile options for GCC and Clang on Unix-like systems
  target_compile_options(PapillonNDL PRIVATE -W -Wall -Wextra -Wpedantic)
endif()

# If building tests, add the tests subdirectory
if(PNDL_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# If building the Python bindings
if(PNDL_PYTHON)
  find_package(pybind11 QUIET)
  if(NOT pybind11_FOUND)
    message(STATUS "Could not find a local install of Pybind11")
    message(STATUS "Will download Pybind11 instead")

    include(FetchContent)

    FetchContent_Declare(pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11
      GIT_TAG        v2.6.2
    )

    FetchContent_MakeAvailable(pybind11)
  else()
    message(STATUS "Using local install of Pybind11")
  endif()

  find_package(Python REQUIRED)

  pybind11_add_module(pyPapillonNDL src/python/pyPapillonNDL.cpp
                                    src/python/interpolation.cpp
                                    src/python/ace.cpp
                                    src/python/function_1d.cpp
                                    src/python/frame.cpp
                                    src/python/angle_law.cpp
                                    src/python/angle_distribution.cpp
                                    src/python/energy_law.cpp
                                    src/python/angle_energy.cpp
                                    src/python/thermal_scattering.cpp
                                    src/python/pctable.cpp
                                    src/python/energy_grid.cpp
                                    src/python/cross_section.cpp
                                    src/python/reaction.cpp
                                    src/python/delayed_group.cpp
                                    src/python/ce_neutron.cpp
                                    src/python/prng.cpp
  )
  
  # Require C++17 standard
  target_compile_features(pyPapillonNDL PRIVATE cxx_std_17)
  
  # Must link Python bindings with full library
  target_link_libraries(pyPapillonNDL PRIVATE PapillonNDL)

  # Get location for python lib
  execute_process(
    COMMAND ${Python_EXECUTABLE}
        -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib(True, False, ''))"
    OUTPUT_VARIABLE PNDL_PYTHON_DEFAULT_RELATIVE_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Inform user where python lib will be going
  message(NOTICE "\nInstalling pyPapillonNDL to " ${CMAKE_INSTALL_PREFIX}/${PNDL_PYTHON_DEFAULT_RELATIVE_PATH}
          "\nEnsure this directory is included in your PYTHONPATH environment variable !\n")

  # Install pyPapillonNDL
  configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/pyPapillonNDLConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/pyPapillonNDLConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pyPapillonNDL
  )

  install(TARGETS pyPapillonNDL
          EXPORT pyPapillonNDLTargets
          LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${PNDL_PYTHON_DEFAULT_RELATIVE_PATH}
  )

  install(EXPORT pyPapillonNDLTargets
          FILE pyPapillonNDLTargets.cmake
          NAMESPACE pyPapillonNDL::
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pyPapillonNDL
  )

  install(FILES
            "${CMAKE_CURRENT_BINARY_DIR}/pyPapillonNDLConfig.cmake"
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pyPapillonNDL
  )
endif()

# Install PapillonNDL
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/PapillonNDLConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/PapillonNDLConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PapillonNDL
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/PapillonNDLConfigVersion.cmake"
  VERSION "${PapillonNDL_VERSION}"
  COMPATIBILITY SameMinorVersion
)

install(TARGETS PapillonNDL
        EXPORT PapillonNDLTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/PapillonNDL DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/)

install(EXPORT PapillonNDLTargets
        FILE PapillonNDLTargets.cmake
        NAMESPACE PapillonNDL::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PapillonNDL
)

install(FILES
          "${CMAKE_CURRENT_BINARY_DIR}/PapillonNDLConfig.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/PapillonNDLConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PapillonNDL
)
