# set up project
project("dune-pdelab" C CXX)

# general stuff
cmake_minimum_required(VERSION 3.6)

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

# start a dune project with information from dune.module
dune_project()

# Disable the upstream testing magic, that builds tests during "make test"
set(DUNE_TEST_MAGIC OFF)

# Use the enable all packages feature
dune_enable_all_packages(MODULE_LIBRARIES dunepdelab)

if(NOT MINGW)
  set(clock_cc "dune/pdelab/common/clock.cc")
endif()

dune_library_add_sources(dunepdelab SOURCES
  ${clock_cc}
  dune/pdelab/common/hostname.cc
  dune/pdelab/common/logtag.cc
)

add_subdirectory(doc)
add_subdirectory(dune/pdelab)
add_subdirectory(cmake/modules)
add_subdirectory(lib)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
