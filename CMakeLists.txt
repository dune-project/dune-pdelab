# set up project
project("dune-pdelab" C CXX)

# general stuff
cmake_minimum_required(VERSION 3.13)

# configuration options
option(DUNE_PDELAB_ENABLE_TRACING "Enable tracing with Perfetto")

# guess dune-common build dir
if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

# find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${dune-common_MODULE_PATH}")

# include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module, and find dune dependencies
dune_project()

# Disable the upstream testing magic, that builds tests during "make test"
set(DUNE_TEST_MAGIC OFF)

# create main library target and set up its properties
dune_add_library(dunepdelab EXPORT_NAME PDELab)

target_link_libraries(dunepdelab
  PUBLIC
    ${DUNE_LIBS}
    $<$<BOOL:${TBB_FOUND}>:TBB::tbb>
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
    $<$<BOOL:${DUNE_PDELAB_ENABLE_TRACING}>:Perfetto::SDK>
)

target_sources(dunepdelab
  PRIVATE
    $<$<NOT:$<BOOL:${MINGW}>>:dune/pdelab/common/clock.cc>
    $<$<BOOL:${DUNE_PDELAB_ENABLE_TRACING}>:dune/pdelab/common/trace.cc>
    dune/pdelab/common/hostname.cc
    dune/pdelab/common/logtag.cc
    dune/pdelab/common/dictionary.cc
)

add_subdirectory(doc)
add_subdirectory(dune)
add_subdirectory(cmake/modules)
add_subdirectory(lib)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
