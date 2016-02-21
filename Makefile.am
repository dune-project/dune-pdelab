# $Id: duneproject 5370 2008-11-15 17:43:24Z christi $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST = 					\
	CHANGELOG.md				\
	CMakeLists.txt				\
	COPYING					\
	dune.module                             \
	GPL-2                                   \
	GPL-3                                   \
	LGPL-3                                  \
	LICENSE.md				\
	REAMDE.md

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# All subdirectories that must (or can) be handled before the library is built
SUBDIRS =					\
        cmake                                   \
	doc					\
	dune					\
	m4
# The library subdirectory itself
SUBDIRS += lib
# All subdirectories that must (or can) be handled after the library is built
# none yet...

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

# Distribute and install config.h.cmake
configdir = $(datadir)/dune-pdelab
dist_config_DATA = config.h.cmake

# Generate package configuration files for finding
# installed modules with CMake
include $(top_srcdir)/am/cmake-pkg-config
