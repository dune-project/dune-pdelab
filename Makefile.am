# $Id: duneproject 5370 2008-11-15 17:43:24Z christi $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)

noinst_PROGRAMS = dune_pdelab

dune_pdelab_SOURCES = dune_pdelab.cc

dune_pdelab_CXXFLAGS = $(MPI_CPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS)
dune_pdelab_LDADD = $(MPI_LDFLAGS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

