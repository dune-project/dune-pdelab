# $Id: duneproject 5370 2008-11-15 17:43:24Z christi $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = doc dune

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

