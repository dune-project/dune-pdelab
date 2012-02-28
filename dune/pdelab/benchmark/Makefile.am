# $Id: Makefile.am 586 2006-04-28 13:41:31Z mblatt $

#defined empty so we can add to it later
NORMALTESTS = 
MOSTLYCLEANFILES =
check_SCRIPTS = 

# list of tests to run
TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

# UG flags must be in front of ALUGrid flags, else when
# dune/grid/uggrid/ugincludes.hh does #include "parallel.h",
# ALUGrid-1.14-inst/include/serial/parallel.h will be included...
AM_CPPFLAGS = @AM_CPPFLAGS@				\
	$(UG_CPPFLAGS)					\
	$(ALUGRID_CPPFLAGS)				\
	$(AMIRAMESH_CPPFLAGS)				\
	$(DUNE_CPPFLAGS)				\
	$(DUNEMPICPPFLAGS)				\
	$(PARMETIS_CPPFLAGS)
AM_LDFLAGS = @AM_LDFLAGS@			\
	$(UG_LDFLAGS)				\
	$(ALUGRID_LDFLAGS)			\
	$(AMIRAMESH_LDFLAGS)			\
	$(DUNE_LDFLAGS)				\
	$(DUNEMPILDFLAGS)			\
	$(PARMETIS_LDFLAGS)
LDADD =						\
	$(PARMETIS_LIBS)			\
	$(ALUGRID_LIBS)				\
	$(AMIRAMESH_LIBS)			\
	$(DUNE_LIBS)				\
	$(DUNEMPILIBS)				\
	$(UG_LIBS)				\
	@LDADD@

#noinst_HEADERS =

# THIS IS A SEPARATION MARKER TO MINIMIZE SVN MERGE CONFLICTS
# TRUNK TARGETS FOLLOW

# THIS IS A SEPARATION MARKER TO MINIMIZE SVN MERGE CONFLICTS
# BRANCH TARGETS FOLLOW

NORMALTESTS += poisson
poisson_SOURCES = poisson.cc
MOSTLYCLEANFILES += poisson_*.vtu poisson_*timings.txt

include $(top_srcdir)/am/global-rules

DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"
