# set up project
project("dune-pdelab" C CXX)

# general stuff
cmake_minimum_required(VERSION 2.8.12)

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

# To have enable all packages work with cmake-2.8, the module
# libraries need to be added before the call to dune_enable_all_packages()
# With cmake 3.1, one would use the MODULE_LIBRARIES argument of
# dune_enable_all_packages() instead.
dune_add_library(dunepdelab
  dune/pdelab/common/clock.cc
  dune/pdelab/common/hostname.cc
  dune/pdelab/common/logtag.cc
  )
dune_target_enable_all_packages(dunepdelab)

dune_register_package_flags(LIBRARIES dunepdelab)

# Use the enable all packages feature
dune_enable_all_packages()

add_subdirectory(doc)
add_subdirectory(dune)
add_subdirectory(cmake/modules)
add_subdirectory(lib)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
