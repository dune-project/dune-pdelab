# general stuff
cmake_minimum_required(VERSION 3.13)

# set up project
project(dune-pdelab C CXX)

# guess dune-common build dir
if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# deactivate global include-directories
if (dune-common_VERSION VERSION_GREATER_EQUAL 2.11)
  dune_policy(SET DP_DEFAULT_INCLUDE_DIRS NEW)
endif()

# start a dune project with information from dune.module
dune_project()

# Create a target for dune-pdelab with a Dune::PDELab alias
dune_add_library(dunepdelab EXPORT_NAME PDELab)
target_link_libraries(dunepdelab PUBLIC ${DUNE_LIBS})

# set include directories to target
if (dune-common_VERSION VERSION_GREATER_EQUAL 2.11)
  dune_default_include_directories(dunepdelab PUBLIC)
endif()

if(NOT MINGW)
  set(clock_cc "dune/pdelab/common/clock.cc")
endif()

target_sources(dunepdelab PRIVATE
  ${clock_cc}
  dune/pdelab/common/hostname.cc
  dune/pdelab/common/logtag.cc
)

# Configuration options
option(COVERAGE_REPORT "Enable coverage compiler flags")

if(COVERAGE_REPORT)
  target_compile_options(dunepdelab PUBLIC --coverage)
  target_link_libraries(dunepdelab PUBLIC gcov)
endif()

add_dune_metis_flags(dunepdelab)

if (METIS_FOUND)
  list(APPEND DUNE_PDELAB_PACKAGE_DEPENDENCIES "find_dependency(METIS)")
endif()


add_subdirectory(doc)
add_subdirectory(dune)
add_subdirectory(cmake/modules)
add_subdirectory(lib)

# write contents into DUNE_CUSTOM_PKG_CONFIG_SECTION, which will be injected into dune-common-config.cmake
string(JOIN "\n" DUNE_CUSTOM_PKG_CONFIG_SECTION
  # make sure that Find<module>.cmake provided by dune-common can be found by cmake
  [[list(APPEND CMAKE_MODULE_PATH "${dune-common_MODULE_PATH}")]]
  ${DUNE_PDELAB_PACKAGE_DEPENDENCIES}
)


# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project()
