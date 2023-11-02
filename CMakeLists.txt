# general stuff
cmake_minimum_required(VERSION 3.13)

# set up project
project("dune-pdelab" C CXX)

# configuration options
option(DUNE_PDELAB_ENABLE_TRACING "Enable tracing with Perfetto")
option(DUNE_PDELAB_USE_BOOST_ATOMIC_REF "Use boost atomic_ref" OFF)

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
    $<TARGET_NAME_IF_EXISTS:TBB::tbb>
    $<TARGET_NAME_IF_EXISTS:OpenMP::OpenMP_CXX>
    $<TARGET_NAME_IF_EXISTS:Dune::Perfetto::SDK>
)

target_sources(dunepdelab
  PRIVATE
    $<$<NOT:$<BOOL:${MINGW}>>:dune/pdelab/common/clock.cc>
    $<$<TARGET_EXISTS:Dune::Perfetto::SDK>:dune/pdelab/common/trace.cc>
    dune/pdelab/common/hostname.cc
    dune/pdelab/common/logtag.cc
    dune/pdelab/common/property.cc
    dune/pdelab/common/property_tree.cc
)

target_compile_features(dunepdelab
  PUBLIC
    cxx_std_20
)

if(DUNE_PDELAB_USE_BOOST_ATOMIC_REF)
    target_compile_definitions(dunepdelab PUBLIC DUNE_PDELAB_USE_BOOST_ATOMIC_REF)
endif()

add_subdirectory(doc)
add_subdirectory(dune)
add_subdirectory(cmake/modules)
add_subdirectory(lib)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
