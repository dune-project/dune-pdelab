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
	$(MPI_CPPFLAGS)
LIBS = @LIBS@					\
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS)	\
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS)	\
	$(DUNE_LDFLAGS) $(DUNE_LIBS)		\
	$(MPI_LDFLAGS) $(MPILIBS)		\
	$(UG_LDFLAGS) $(UG_LIBS)

noinst_HEADERS =				\
	divergence-probe.hh			\
	electricenergy-probe.hh			\
	gnuplotgraph.hh				\
	gridexamples.hh				\
	l2difference.hh				\
	l2norm.hh				\
	physicalconstants.hh			\
	probe.hh				\
	resonatorsolution.hh

noinst_SCRIPTS =				\
	make_pvd.sh

# scripts are not distributed automatically
# but since this is need for make all:
EXTRA_DIST = make_pvd.sh

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

NORMALTESTS += testedges02dgridfunctionspace
testedges02dgridfunctionspace_SOURCES = testedges02dgridfunctionspace.cc
testedges02dgridfunctionspace_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=2 -DENABLE_ALBERTA
testedges02dgridfunctionspace_LDADD =				\
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

NORMALTESTS += testedges03dgridfunctionspaceglobal
testedges03dgridfunctionspaceglobal_SOURCES = testedges03dgridfunctionspaceglobal.cc
testedges03dgridfunctionspaceglobal_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA
testedges03dgridfunctionspaceglobal_LDADD =				\
	$(ALBERTA_LDFLAGS) -lalberta_3d $(ALBERTA_BASE_LIBS)
MOSTLYCLEANFILES += edges03dgridfunctionspaceglobal-*.vtu

NORMALTESTS += testedges03dinterpolationglobal
testedges03dinterpolationglobal_SOURCES = testedges03dinterpolationglobal.cc
testedges03dinterpolationglobal_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA
testedges03dinterpolationglobal_LDADD = 			\
	$(ALBERTA_LDFLAGS) -lalberta_3d $(ALBERTA_BASE_LIBS)
MOSTLYCLEANFILES += edges03dinterpolationglobal-*.vtu	\
	edges03dinterpolationglobal-*.dat		\
	edges03dinterpolationglobal.eps			\
	testedges03dinterpolationglobal.gnuplot

NORMALTESTS += testedges03dfem
testedges03dfem_SOURCES = testedges03dfem.cc
testedges03dfem_CPPFLAGS = $(AM_CPPFLAGS)				\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA
testedges03dfem_LDADD =						\
	$(ALBERTA_LDFLAGS) -lalberta_3d $(ALBERTA_BASE_LIBS)

NORMALTESTS += testelectrodynamic
testelectrodynamic_SOURCES = testelectrodynamic.cc
testelectrodynamic_CPPFLAGS = $(AM_CPPFLAGS)				\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA	\
	$(SUPERLU_CPPFLAGS)
testelectrodynamic_LDADD =					\
	$(ALBERTA_LDFLAGS) -lalberta_3d $(ALBERTA_BASE_LIBS)	\
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)
MOSTLYCLEANFILES +=							\
	electrodynamic-alberta-triangulated-cube-6.level*.pvd		\
	electrodynamic-alu-triangulated-cube-6.level*.pvd		\
	electrodynamic-electricenergy-evolution.dat			\
	electrodynamic-electricenergy-evolution.eps			\
	electrodynamic-electricenergy-evolution.gnuplot			\
	electrodynamic-globalerror.dat					\
	electrodynamic-globalerror.eps					\
	electrodynamic-globalerror.gnuplot				\
	electrodynamic-globalevolution.dat				\
	electrodynamic-globalevolution.eps				\
	electrodynamic-globalevolution.gnuplot				\
	electrodynamic-l2error.dat					\
	electrodynamic-l2error.eps					\
	electrodynamic-l2error.gnuplot					\
	electrodynamic-l2evolution.dat					\
	electrodynamic-l2evolution.eps					\
	electrodynamic-l2evolution.gnuplot				\
	electrodynamic-probe-*.eps					\
	electrodynamic-probe.dat					\
	electrodynamic-probe.gnuplot					\
	electrodynamic-ug-*.level*.pvd					\
	s0001:electrodynamic-alberta-triangulated-cube-6.level*-*.pvtu	\
	s0001:electrodynamic-alu-triangulated-cube-6.level*-*.pvtu	\
	s0001:electrodynamic-ug-triangulated-cube-6.level*-*.pvtu	\
	s0001:electrodynamic-ug-triangulated-square.level*-*.pvtu	\
	s0001:p0000:electrodynamic-alberta-triangulated-cube-6.level*-*.vtu \
	s0001:p0000:electrodynamic-alu-triangulated-cube-6.level*-*.vtu	\
	s0001:p0000:electrodynamic-ug-triangulated-cube-6.level*-*.vtu	\
	s0001:p0000:electrodynamic-ug-triangulated-square.level*-*.vtu	\
	testelectrodynamic

NORMALTESTS += testelectrodynamic-newmark
testelectrodynamic_newmark_SOURCES = testelectrodynamic-newmark.cc
testelectrodynamic_newmark_CPPFLAGS = $(AM_CPPFLAGS)			\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA	\
	$(SUPERLU_CPPFLAGS)
testelectrodynamic_newmark_LDADD =				\
	$(ALBERTA_LDFLAGS) -lalberta_3d $(ALBERTA_BASE_LIBS)	\
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)
MOSTLYCLEANFILES +=								\
	electrodynamic-newmark-electricenergy-evolution.eps			\
	electrodynamic-newmark-electricenergy-evolution.dat			\
	electrodynamic-newmark-electricenergy-evolution.gnuplot			\
	electrodynamic-newmark-globalerror.dat					\
	electrodynamic-newmark-globalerror.eps					\
	electrodynamic-newmark-globalerror.gnuplot				\
	electrodynamic-newmark-globalevolution.dat				\
	electrodynamic-newmark-globalevolution.eps				\
	electrodynamic-newmark-globalevolution.gnuplot				\
	electrodynamic-newmark-l2error.dat					\
	electrodynamic-newmark-l2error.eps					\
	electrodynamic-newmark-l2error.gnuplot					\
	electrodynamic-newmark-l2evolution.dat					\
	electrodynamic-newmark-l2evolution.eps					\
	electrodynamic-newmark-l2evolution.gnuplot				\
	electrodynamic-newmark-probe.dat					\
	electrodynamic-newmark-probe-*.eps					\
	electrodynamic-newmark-probe.gnuplot					\
	electrodynamic-newmark-ug-triangulated-cube-6.level*.pvd		\
	s0001:electrodynamic-newmark-ug-triangulated-cube-6.level*-*.pvtu	\
	s0001:p0000:electrodynamic-newmark-ug-triangulated-cube-6.level*-*.vtu	\
	testelectrodynamic-newmark

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
	$(ALBERTA_CPPFLAGS)				\
	$(SUPERLU_CPPFLAGS)
testlaplacedirichletccfv_LDADD =		\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)	\
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)
MOSTLYCLEANFILES += testlaplacedirichletccfv.vtu

NORMALTESTS += testlaplacedirichletp12d
testlaplacedirichletp12d_SOURCES = testlaplacedirichletp12d.cc
testlaplacedirichletp12d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)				\
	$(SUPERLU_CPPFLAGS)
testlaplacedirichletp12d_LDADD =		\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)	\
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)
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

NORMALTESTS += testrt02dgridfunctionspace
testrt02dgridfunctionspace_SOURCES = testrt02dgridfunctionspace.cc
testrt02dgridfunctionspace_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testrt02dgridfunctionspace_LDADD = 		\
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
