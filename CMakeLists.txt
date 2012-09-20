# set up project
project("dune-pdelab" C CXX)

# general stuff
cmake_minimum_required(VERSION 2.8)

# make sure our own modules are found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# $Id: duneproject 5370 2008-11-15 17:43:24Z christi $

# we need the module file to be able to build via dunecontrol
set(EXTRA_DIST                     
        dune.module             
        GPL-2                   
        GPL-3                   
        LGPL-3                  
        RELEASE_NOTES)

# don't follow the full GNU-standard
# we need automake 1.5
set(AUTOMAKE_OPTIONS  foreign 1.5)

# All subdirectories that must (or can) be handled before the library is built
set(SUBDIRS                                        
        doc                                     
        dune                                    
        m4)
# The library subdirectory itself
list(APPEND SUBDIRS lib)
# All subdirectories that must (or can) be handled after the library is built
# none yet...

# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune=${DUNEROOT} CXX="${CXX}" CXXFLAGS="${CXXFLAGS}" CC="${CC}")

# include not needed for CMake
# include $(top_srcdir)/am/top-rules
# include not needed for CMake
# include $(top_srcdir)/am/global-rules

foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})

add_subdirectory(cmake/modules)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
