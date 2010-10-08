# $Id: Makefile.am 586 2006-04-28 13:41:31Z mblatt $

#defined empty so we can add to it later
NORMALTESTS = 
MOSTLYCLEANFILES =

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

noinst_HEADERS =				\
	fmt.hh					\
	gnuplotgraph.hh				\
	gridexamples.hh				\
	l2difference.hh				\
	l2norm.hh

noinst_SCRIPTS =				\
	make_pvd.sh

GRIDDIM = 2

# scripts are not distributed automatically
# but since this is need for make all:
EXTRA_DIST = make_pvd.sh

NORMALTESTS += testanalytic
testanalytic_SOURCES = testanalytic.cc
testanalytic_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testanalytic_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)
testanalytic_LDADD =				\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += channel.vtu

NORMALTESTS += testconstraints
testconstraints_SOURCES = testconstraints.cc
testconstraints_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testconstraints_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)
testconstraints_LDADD =				\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += testconstraintsp1.vtu testconstraintspowerp1.vtu

NORMALTESTS += testcountingptr
testcountingptr_SOURCES = testcountingptr.cc
testcountingptr_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testcountingptr_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)
testcountingptr_LDADD =				\
	$(ALBERTA_LIBS)				\
	$(LDADD)

NORMALTESTS += testedges02d
testedges02d_SOURCES = testedges02d.cc
testedges02d_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testedges02d_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)
testedges02d_LDADD =				\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += testedges02d.vtu

NORMALTESTS += testedges02dgridfunctionspace
testedges02dgridfunctionspace_SOURCES = testedges02dgridfunctionspace.cc
testedges02dgridfunctionspace_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA2D_CPPFLAGS)
testedges02dgridfunctionspace_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA2D_LDFLAGS)
testedges02dgridfunctionspace_LDADD =				\
	$(ALBERTA2D_LIBS)			\
	$(LDADD)
MOSTLYCLEANFILES += edges02dgridfunctionspace-*.vtu

NORMALTESTS += testedges02dinterpolation
testedges02dinterpolation_SOURCES = testedges02dinterpolation.cc
testedges02dinterpolation_CPPFLAGS = $(AM_CPPFLAGS)			\
	$(ALBERTA2D_CPPFLAGS)
testedges02dinterpolation_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA2D_LDFLAGS)
testedges02dinterpolation_LDADD =				\
	$(ALBERTA2D_LIBS)			\
	$(LDADD)
MOSTLYCLEANFILES += edges02dinterpolation-*.vtu

NORMALTESTS += testedges02dinterpolationglobal
testedges02dinterpolationglobal_SOURCES = testedges02dinterpolationglobal.cc
testedges02dinterpolationglobal_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA2D_CPPFLAGS)
testedges02dinterpolationglobal_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA2D_LDFLAGS)
testedges02dinterpolationglobal_LDADD = 			\
	$(ALBERTA2D_LIBS)			\
	$(LDADD)
MOSTLYCLEANFILES += edges02dinterpolationglobal-*.vtu

NORMALTESTS += testedges03dgridfunctionspaceglobal
testedges03dgridfunctionspaceglobal_SOURCES = testedges03dgridfunctionspaceglobal.cc
testedges03dgridfunctionspaceglobal_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA3D_CPPFLAGS)
testedges03dgridfunctionspaceglobal_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA3D_LDFLAGS)
testedges03dgridfunctionspaceglobal_LDADD =				\
	$(ALBERTA3D_LIBS)			\
	$(LDADD)
MOSTLYCLEANFILES += edges03dgridfunctionspaceglobal-*.vtu

NORMALTESTS += testedges03dinterpolationglobal
testedges03dinterpolationglobal_SOURCES = testedges03dinterpolationglobal.cc
testedges03dinterpolationglobal_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA3D_CPPFLAGS)
testedges03dinterpolationglobal_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA3D_LDFLAGS)
testedges03dinterpolationglobal_LDADD = 			\
	$(ALBERTA3D_LIBS)			\
	$(LDADD)
MOSTLYCLEANFILES += edges03dinterpolationglobal-*.vtu	\
	edges03dinterpolationglobal-*.dat		\
	edges03dinterpolationglobal.eps			\
	testedges03dinterpolationglobal.gnuplot

NORMALTESTS += testedges03dfem
testedges03dfem_SOURCES = testedges03dfem.cc
testedges03dfem_CPPFLAGS = $(AM_CPPFLAGS)				\
	$(ALBERTA3D_CPPFLAGS)
testedges03dfem_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA3D_LDFLAGS)
testedges03dfem_LDADD =						\
	$(ALBERTA3D_LIBS)			\
	$(LDADD)

NORMALTESTS += testfiniteelementmap
testfiniteelementmap_SOURCES = testfiniteelementmap.cc
testfiniteelementmap_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testfiniteelementmap_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA_LDFLAGS)
testfiniteelementmap_LDADD =			\
	$(ALBERTA_LIBS)				\
	$(LDADD)

NORMALTESTS += testfunction
testfunction_SOURCES = testfunction.cc
testfunction_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testfunction_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)
testfunction_LDADD =				\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += multi.vtu single.vtu

NORMALTESTS += testgridfunctionspace
testgridfunctionspace_SOURCES = testgridfunctionspace.cc

NORMALTESTS += testlaplacedirichletccfv
testlaplacedirichletccfv_SOURCES = testlaplacedirichletccfv.cc
testlaplacedirichletccfv_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)				\
	$(SUPERLU_CPPFLAGS)
testlaplacedirichletccfv_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA_LDFLAGS)				\
	$(SUPERLU_LDFLAGS)
testlaplacedirichletccfv_LDADD =		\
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)	\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += testlaplacedirichletccfv.vtu

NORMALTESTS += testlaplacedirichletp12d
testlaplacedirichletp12d_SOURCES = testlaplacedirichletp12d.cc
testlaplacedirichletp12d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)				\
	$(SUPERLU_CPPFLAGS)
testlaplacedirichletp12d_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA_LDFLAGS)				\
	$(SUPERLU_LDFLAGS)
testlaplacedirichletp12d_LDADD =		\
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)	\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += testlaplacedirichletp12d.vtu

NORMALTESTS += testlocalfunctionspace
testlocalfunctionspace_SOURCES = testlocalfunctionspace.cc
testlocalfunctionspace_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testlocalfunctionspace_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA_LDFLAGS)
testlocalfunctionspace_LDADD =			\
	$(ALBERTA_LIBS)				\
	$(LDADD)

NORMALTESTS += testmultitypetree
testmultitypetree_SOURCES = testmultitypetree.cc
testmultitypetree_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testmultitypetree_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA_LDFLAGS)
testmultitypetree_LDADD = 			\
	$(ALBERTA_LIBS)				\
	$(LDADD)

NORMALTESTS += testmultistep
testmultistep_SOURCES = testmultistep.cc
testmultistep_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(DUNEMPICPPFLAGS)
testmultistep_LDFLAGS = $(AM_LDFLAGS)		\
	$(DUNEMPILDFLAGS)
testmultistep_LDADD =				\
	$(DUNEMPILIBS)				\
	$(LDADD)

MOSTLYCLEANFILES +=					\
	s*:p*:testmultistep_yasp_P1_1d-*.vtp		\
	s*:testmultistep_yasp_P1_1d-*.pvtp		\
	testmultistep_yasp_P1_1d.pvd

NORMALTESTS += testp12dinterpolation
testp12dinterpolation_SOURCES = testp12dinterpolation.cc
testp12dinterpolation_CPPFLAGS = $(AM_CPPFLAGS)				\
	$(ALBERTA2D_CPPFLAGS)
testp12dinterpolation_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA2D_LDFLAGS)
testp12dinterpolation_LDADD =					\
	$(ALBERTA2D_LIBS)			\
	$(LDADD)

NORMALTESTS += testpk
testpk_SOURCES = testpk.cc
testpk_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testpk_LDFLAGS = $(AM_LDFLAGS)			\
	$(ALBERTA_LDFLAGS)
testpk_LDADD =					\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += testpk.vtu

NORMALTESTS += testpoisson
testpoisson_SOURCES = testpoisson.cc
testpoisson_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)			\
	$(SUPERLU_CPPFLAGS)
testpoisson_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)			\
	$(SUPERLU_LDFLAGS)
testpoisson_LDADD =				\
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)	\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += poisson_*.vtu

NORMALTESTS += testpoisson-globalfe
testpoisson_globalfe_SOURCES = testpoisson-globalfe.cc
testpoisson_globalfe_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)			\
	$(SUPERLU_CPPFLAGS)
testpoisson_globalfe_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA_LDFLAGS)			\
	$(SUPERLU_LDFLAGS)
testpoisson_globalfe_LDADD =			\
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)	\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += poisson_globalfe_*.vtu

NORMALTESTS += testrt0
testrt0_SOURCES = testrt0.cc
testrt0_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testrt0_LDFLAGS = $(AM_LDFLAGS)			\
	$(ALBERTA_LDFLAGS)
testrt0_LDADD = 				\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += testrt0.vtu

NORMALTESTS += testrt02dgridfunctionspace
testrt02dgridfunctionspace_SOURCES = testrt02dgridfunctionspace.cc
testrt02dgridfunctionspace_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testrt02dgridfunctionspace_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA_LDFLAGS)
testrt02dgridfunctionspace_LDADD = 		\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += rt02dgridfunctionspace-*.vtu

NORMALTESTS += testutilities
testutilities_SOURCES = testutilities.cc
testutilities_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testutilities_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)
testutilities_LDADD = 				\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += edger.vtu interpolated.vtu q1.vtu taylorhood.vtu

NORMALTESTS += testvectorwave
testvectorwave_SOURCES = testvectorwave.cc
MOSTLYCLEANFILES +=				\
	s*:p*:vectorwave_*.vtu			\
	s*:vectorwave_*.vtup			\
	vectorwave_*.pvd

include $(top_srcdir)/am/global-rules

DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"
