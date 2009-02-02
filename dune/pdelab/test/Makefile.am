# $Id: Makefile.am 586 2006-04-28 13:41:31Z mblatt $

NORMALTESTS = testcountingptr

# list of tests to run
TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

testcountingptr_SOURCES = testcountingptr.cc
testcountingptr_CXXFLAGS = $(MPI_CPPFLAGS)
testcountingptr_LDADD = $(MPI_LDFLAGS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

include $(top_srcdir)/am/global-rules

DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"
