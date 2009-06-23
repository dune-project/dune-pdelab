# $Id: Makefile.am 586 2006-04-28 13:41:31Z mblatt $

#defined empty so we can add to it later
NORMALTESTS = 
MOSTLYCLEANFILES =

# list of tests to run
TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

AM_CPPFLAGS = @AM_CPPFLAGS@			\
	$(ALUGRID_CPPFLAGS)			\
	$(AMIRAMESH_CPPFLAGS)			\
	$(DUNE_CPPFLAGS)			\
	$(MPI_CPPFLAGS)				\
	$(UG_CPPFLAGS)
LIBS = @LIBS@					\
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS)	\
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS)	\
	$(DUNE_LDFLAGS) $(DUNE_LIBS)		\
	$(MPI_LDFLAGS) $(MPILIBS)		\
	$(UG_LDFLAGS) $(UG_LIBS)

noinst_HEADERS = gridexamples.hh stddomains.hh


NORMALTESTS += testanalytic
testanalytic_SOURCES = testanalytic.cc
testanalytic_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testanalytic_LDADD =				\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)
MOSTLYCLEANFILES += channel.vtu

NORMALTESTS += testconstraints
testconstraints_SOURCES = testconstraints.cc
testconstraints_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testconstraints_LDADD = 			\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)
MOSTLYCLEANFILES += testconstraintsp1.vtu testconstraintspowerp1.vtu

NORMALTESTS += testcountingptr
testcountingptr_SOURCES = testcountingptr.cc
testcountingptr_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testcountingptr_LDADD = 			\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) 

NORMALTESTS += testedges02d
testedges02d_SOURCES = testedges02d.cc
testedges02d_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testedges02d_LDADD =				\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)
MOSTLYCLEANFILES += testedges02d.vtu

NORMALTESTS += testedges02dgridfunctionsspace
testedges02dgridfunctionsspace_SOURCES = testedges02dgridfunctionsspace.cc
testedges02dgridfunctionsspace_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=2 -DENABLE_ALBERTA
testedges02dgridfunctionsspace_LDADD =				\
	$(ALBERTA_LDFLAGS) -lalberta_2d $(ALBERTA_BASE_LIBS)
MOSTLYCLEANFILES += edges02dgridfunctionspace-*.vtu

NORMALTESTS += testedges02dinterpolation
testedges02dinterpolation_SOURCES = testedges02dinterpolation.cc
testedges02dinterpolation_CPPFLAGS = $(AM_CPPFLAGS)			\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=2 -DENABLE_ALBERTA 
testedges02dinterpolation_LDADD =				\
	$(ALBERTA_LDFLAGS) -lalberta_2d $(ALBERTA_BASE_LIBS)
MOSTLYCLEANFILES += edges02dinterpolation-*.vtu

NORMALTESTS += testedges02dinterpolationglobal
testedges02dinterpolationglobal_SOURCES = testedges02dinterpolationglobal.cc
testedges02dinterpolationglobal_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=2 -DENABLE_ALBERTA
testedges02dinterpolationglobal_LDADD = 			\
	$(ALBERTA_LDFLAGS) -lalberta_2d $(ALBERTA_BASE_LIBS)
MOSTLYCLEANFILES += edges02dinterpolationglobal-*.vtu

NORMALTESTS += testedges03dgridfunctionsspace
testedges03dgridfunctionsspace_SOURCES = testedges03dgridfunctionsspace.cc
testedges03dgridfunctionsspace_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA
testedges03dgridfunctionsspace_LDADD =				\
	$(ALBERTA_LDFLAGS) -lalberta_3d $(ALBERTA_BASE_LIBS)
MOSTLYCLEANFILES += edges03dgridfunctionspace-*.vtu

NORMALTESTS += testedges03dinterpolationglobal
testedges03dinterpolationglobal_SOURCES = testedges03dinterpolationglobal.cc
testedges03dinterpolationglobal_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA
testedges03dinterpolationglobal_LDADD = 			\
	$(ALBERTA_LDFLAGS) -lalberta_3d $(ALBERTA_BASE_LIBS)
MOSTLYCLEANFILES += edges03dinterpolationglobal-*.vtu

NORMALTESTS += testelectrostatic
testelectrostatic_SOURCES = testelectrostatic.cc
testelectrostatic_CPPFLAGS = $(AM_CPPFLAGS)				\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA
testelectrostatic_LDADD =					\
	$(ALBERTA_LDFLAGS) -lalberta_3d $(ALBERTA_BASE_LIBS)
MOSTLYCLEANFILES += electrostatic_*.vtu

NORMALTESTS += testfiniteelementmap
testfiniteelementmap_SOURCES = testfiniteelementmap.cc
testfiniteelementmap_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testfiniteelementmap_LDADD = 			\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)

NORMALTESTS += testfunction
testfunction_SOURCES = testfunction.cc
testfunction_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testfunction_LDADD = 				\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)
MOSTLYCLEANFILES += multi.vtu single.vtu

NORMALTESTS += testgridfunctionspace
testgridfunctionspace_SOURCES = testgridfunctionspace.cc
testgridfunctionspace_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testgridfunctionspace_LDADD = 			\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)

NORMALTESTS += testlaplacedirichletccfv
testlaplacedirichletccfv_SOURCES = testlaplacedirichletccfv.cc
testlaplacedirichletccfv_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testlaplacedirichletccfv_LDADD =		\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)
MOSTLYCLEANFILES += testlaplacedirichletccfv.vtu

NORMALTESTS += testlaplacedirichletp12d
testlaplacedirichletp12d_SOURCES = testlaplacedirichletp12d.cc
testlaplacedirichletp12d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testlaplacedirichletp12d_LDADD =		\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)
MOSTLYCLEANFILES += testlaplacedirichletp12d.vtu

NORMALTESTS += testlocalfunctionspace
testlocalfunctionspace_SOURCES = testlocalfunctionspace.cc
testlocalfunctionspace_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testlocalfunctionspace_LDADD =			\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)

NORMALTESTS += testmultitypetree
testmultitypetree_SOURCES = testmultitypetree.cc
testmultitypetree_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testmultitypetree_LDADD = 			\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)

NORMALTESTS += testp12dinterpolation
testp12dinterpolation_SOURCES = testp12dinterpolation.cc
testp12dinterpolation_CPPFLAGS = $(AM_CPPFLAGS)				\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=2 -DENABLE_ALBERTA
testp12dinterpolation_LDADD =					\
	$(ALBERTA_LDFLAGS) -lalberta_2d $(ALBERTA_BASE_LIBS)

NORMALTESTS += testpk
testpk_SOURCES = testpk.cc
testpk_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testpk_LDADD =					\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)
MOSTLYCLEANFILES += testpk.vtu

NORMALTESTS += testpoisson
testpoisson_SOURCES = testpoisson.cc
testpoisson_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)			\
	$(SUPERLU_CPPFLAGS)
testpoisson_LDADD =				\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)	\
	$(UG_LDFLAGS) $(UG_LIBS)
MOSTLYCLEANFILES += poisson_*.vtu

NORMALTESTS += testrt0
testrt0_SOURCES = testrt0.cc
testrt0_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testrt0_LDADD = 				\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)
MOSTLYCLEANFILES += testrt0.vtu

NORMALTESTS += testrt02dgridfunctionsspace
testrt02dgridfunctionsspace_SOURCES = testrt02dgridfunctionsspace.cc
testrt02dgridfunctionsspace_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testrt02dgridfunctionsspace_LDADD = 		\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)
MOSTLYCLEANFILES += rt02dgridfunctionspace-*.vtu

NORMALTESTS += testutilities
testutilities_SOURCES = testutilities.cc
testutilities_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testutilities_LDADD = 				\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)
MOSTLYCLEANFILES += edger.vtu interpolated.vtu q1.vtu taylorhood.vtu

include $(top_srcdir)/am/global-rules

DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"
