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

if(NOT MINGW)
  set(clock_cc "dune/pdelab/common/clock.cc")
endif()


# Pull in vendored version of {fmt}
# The original upstream repository can be found in vendor/fmt, so make sure it actually exists,
# as people might have checked out PDELab without the --recursive switch
# We then copy the required headers to a nested location under dune/pdelab to avoid polluting
# the paths which belong to the upstream library and add the two compiled source files to the
# dunepdelab library
if(NOT TARGET fmt::fmt)

  message(STATUS "Could not find installed {fmt} library, vendoring included version")

  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/vendor/fmt/src/format.cc)
    message(FATAL_ERROR "Could not find vendored {fmt} library, did you check out the submodule?")
  endif()

  set(
    FMT_SOURCES
      vendor/fmt/src/format.cc
      vendor/fmt/src/posix.cc
    )

  file(
    GLOB FMT_HEADERS
    vendor/fmt/include/fmt/*.h
    )

  file(
    COPY ${FMT_HEADERS}
    DESTINATION dune/pdelab/vendor/fmt/fmt
    )

  install(
    FILES ${FMT_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/vender/fmt/fmt
    )

endif()

dune_add_library(dunepdelab
  ${clock_cc}
  dune/pdelab/common/hostname.cc
  dune/pdelab/common/logtag.cc
  ${FMT_SOURCES}
  )
dune_target_enable_all_packages(dunepdelab)

if(TARGET fmt::fmt)
  # We are using an upstream version of {fmt}, so link against that
  target_link_libraries(dunepdelab PUBLIC fmt::fmt)
else()

  # We need to set some preprocessor definions
  target_compile_definitions(
    dunepdelab
    PRIVATE FMT_EXPORT
    INTERFACE FMT_SHARED
    PUBLIC DUNE_PDELAB_VENDORED_FMT=1
    )

  # Tell the compiler about the custom path of the fmt headers
  target_include_directories(
    dunepdelab
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/dune/pdelab/vendor/fmt>
      $<INSTALL_INTERFACE:dune/pdelab/vendor/fmt>
    )

endif()

# turn off warning about string literal operator templates caused by _fmt
target_compile_options(
  dunepdelab
  PUBLIC -Wno-gnu-string-literal-operator-template
  )

dune_register_package_flags(LIBRARIES dunepdelab)

# Use the enable all packages feature
dune_enable_all_packages()

add_subdirectory(doc)
add_subdirectory(dune/pdelab)
add_subdirectory(cmake/modules)
add_subdirectory(lib)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
