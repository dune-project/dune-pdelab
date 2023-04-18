# set up project
project("dune-pdelab" C CXX)

# general stuff
cmake_minimum_required(VERSION 3.6)

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
find_package(OpenMP)

set(DUNE_CUSTOM_PKG_CONFIG_SECTION "
include(CMakeFindDependencyMacro)
find_dependency(OpenMP)
")

if (DUNE_PDELAB_ENABLE_TRACING)
  # declare and populate perfetto
  include(FetchContent)
  message("-- Declaring Perfetto")
  FetchContent_Declare(
    perfetto
    GIT_REPOSITORY https://android.googlesource.com/platform/external/perfetto
    GIT_TAG        v32.1
  )

  # configure perfetto targets
  if(NOT perfetto_POPULATED)
    message("-- Populating Perfetto")
    FetchContent_Populate(perfetto)
    # add_library(perfetto-sdk STATIC "${perfetto_SOURCE_DIR}/sdk/perfetto.cc")
    dune_add_library(perfetto-sdk STATIC
      SOURCES "${perfetto_SOURCE_DIR}/sdk/perfetto.cc"
      EXPORT_NAME Perfetto::SDK)
    target_link_libraries(perfetto-sdk PUBLIC $<$<BOOL:${Threads_FOUND}>:Threads::Threads>)
    target_compile_definitions(perfetto-sdk INTERFACE HAVE_PERFETTO)
    install(DIRECTORY ${perfetto_SOURCE_DIR}/sdk DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/external/perfetto-sdk)
    target_include_directories(perfetto-sdk
      PUBLIC
        $<BUILD_INTERFACE:${perfetto_SOURCE_DIR}/sdk>
        $<INSTALL_INTERFACE:dune/external/perfetto-sdk>
      )
  endif()
endif()

# include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module, and find dune dependencies
dune_project()

# Disable the upstream testing magic, that builds tests during "make test"
set(DUNE_TEST_MAGIC OFF)

# create main library target and set up its properties
dune_add_library(dunepdelab EXPORT_NAME PDELab)

target_compile_features(dunepdelab
  PUBLIC
    cxx_std_20
)

target_link_libraries(dunepdelab
  PUBLIC
    ${DUNE_LIBS}
    $<$<BOOL:${TBB_FOUND}>:TBB::tbb>
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
    $<$<BOOL:${DUNE_PDELAB_ENABLE_TRACING}>:Dune::Perfetto::SDK>
)

target_sources(dunepdelab
  PRIVATE
    $<$<NOT:$<BOOL:${MINGW}>>:dune/pdelab/common/clock.cc>
    $<$<BOOL:${DUNE_PDELAB_ENABLE_TRACING}>:dune/pdelab/common/trace.cc>
    dune/pdelab/common/hostname.cc
    dune/pdelab/common/logtag.cc
)

add_subdirectory(doc)
add_subdirectory(dune)
add_subdirectory(cmake/modules)
add_subdirectory(lib)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
