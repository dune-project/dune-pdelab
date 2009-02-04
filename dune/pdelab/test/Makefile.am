# $Id: Makefile.am 586 2006-04-28 13:41:31Z mblatt $

NORMALTESTS = testcountingptr testmultitypetree testfunction testanalytic testfiniteelementmap testgridfunctionspace

# list of tests to run
TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

testcountingptr_SOURCES = testcountingptr.cc
testcountingptr_CXXFLAGS = $(MPI_CPPFLAGS)
testcountingptr_LDADD = $(MPI_LDFLAGS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

testmultitypetree_SOURCES = testmultitypetree.cc
testmultitypetree_CXXFLAGS = $(MPI_CPPFLAGS)
testmultitypetree_LDADD = $(MPI_LDFLAGS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

testfunction_SOURCES = testfunction.cc
testfunction_CXXFLAGS = $(MPI_CPPFLAGS)
testfunction_LDADD = $(MPI_LDFLAGS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

testanalytic_SOURCES = testanalytic.cc
testanalytic_CXXFLAGS = $(MPI_CPPFLAGS)
testanalytic_LDADD = $(MPI_LDFLAGS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

testfiniteelementmap_SOURCES = testfiniteelementmap.cc
testfiniteelementmap_CXXFLAGS = $(MPI_CPPFLAGS)
testfiniteelementmap_LDADD = $(MPI_LDFLAGS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

testgridfunctionspace_SOURCES = testgridfunctionspace.cc
testgridfunctionspace_CXXFLAGS = $(MPI_CPPFLAGS)
testgridfunctionspace_LDADD = $(MPI_LDFLAGS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

include $(top_srcdir)/am/global-rules

DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"
