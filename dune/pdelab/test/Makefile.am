# $Id: Makefile.am 586 2006-04-28 13:41:31Z mblatt $

SUBDIRS =					\
	fem

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
	-DGRIDSDIR="\"$(srcdir)/grids\""			\
	$(UG_CPPFLAGS)					\
	$(DUNE_ALUGRID_CPPFLAGS)				\
	$(AMIRAMESH_CPPFLAGS)				\
	$(DUNE_CPPFLAGS)				\
	$(DUNEMPICPPFLAGS)				\
	$(PARMETIS_CPPFLAGS)
AM_LDFLAGS = @AM_LDFLAGS@			\
	$(UG_LDFLAGS)				\
	$(DUNE_ALUGRID_LDFLAGS)			\
	$(AMIRAMESH_LDFLAGS)			\
	$(DUNE_LDFLAGS)				\
	$(DUNEMPILDFLAGS)			\
	$(PARMETIS_LDFLAGS)
LDADD =						\
	@LDADD@                                 \
	$(PARMETIS_LIBS)			\
	$(DUNE_ALUGRID_LIBS)				\
	$(AMIRAMESH_LIBS)			\
	$(DUNEMPILIBS)				\
	$(UG_LIBS)

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
EXTRA_DIST = CMakeLists.txt make_pvd.sh \
             grids grids/*

NORMALTESTS += testanalytic
testanalytic_SOURCES = testanalytic.cc
testanalytic_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testanalytic_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)
testanalytic_LDADD =				\
	$(LDADD)                                \
	$(ALBERTA_LIBS)

MOSTLYCLEANFILES += channel.vtu

NORMALTESTS += testclock
testclock_SOURCES = testclock.cc
# don't include all the grid stuff
testclock_CPPFLAGS = @AM_CPPFLAGS@
testclock_LDFLAGS = @AM_LDFLAGS@
testclock_LDADD = @LDADD@

# disable IntersectionIndexSet test for now, we know it's broken
# if DUNE_ALUGRID
# NORMALTESTS += test-composed-iis-gfs
# test_composed_iis_gfs_SOURCES = test-composed-iis-gfs.cc
# endif DUNE_ALUGRID

NORMALTESTS += testconstraints
testconstraints_SOURCES = testconstraints.cc
testconstraints_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testconstraints_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)
testconstraints_LDADD =				\
        $(LDADD)                                \
	$(ALBERTA_LIBS)
MOSTLYCLEANFILES += testconstraintsp1.vtu testconstraintspowerp1.vtu

NORMALTESTS += testfunction
testfunction_SOURCES = testfunction.cc
testfunction_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testfunction_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)
testfunction_LDADD =				\
        $(LDADD)                                \
	$(ALBERTA_LIBS)
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
	$(LDADD)                                \
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)	\
	$(ALBERTA_LIBS)
MOSTLYCLEANFILES += testlaplacedirichletccfv-2d.vtu \
	testlaplacedirichletccfv-3d.vtu


# This test is broken right now
# NORMALTESTS += testlaplacedirichletccfv-simplified
# testlaplacedirichletccfv_simplified_SOURCES = $(testlaplacedirichletccfv_SOURCES)
# testlaplacedirichletccfv_simplified_CPPFLAGS = $(testlaplacedirichletccfv_CPPFLAGS)	\
	-DTEST_SIMPLIFIED_INFRASTRUCTURE
# testlaplacedirichletccfv_simplified_LDFLAGS = $(testlaplacedirichletccfv_LDFLAGS)
# testlaplacedirichletccfv_simplified_LDADD = $(testlaplacedirichletccfv_LDADD)
# MOSTLYCLEANFILES += testlaplacedirichletccfv-simplified-2d.vtu \
	testlaplacedirichletccfv-simplified-3d.vtu

# This test still uses GridOperatorSpace
# NORMALTESTS += testlaplacedirichletp12d
# testlaplacedirichletp12d_SOURCES = testlaplacedirichletp12d.cc
# testlaplacedirichletp12d_CPPFLAGS = $(AM_CPPFLAGS)	\
# 	$(ALBERTA_CPPFLAGS)				\
# 	$(SUPERLU_CPPFLAGS)
# testlaplacedirichletp12d_LDFLAGS = $(AM_LDFLAGS)	\
# 	$(ALBERTA_LDFLAGS)				\
# 	$(SUPERLU_LDFLAGS)
# testlaplacedirichletp12d_LDADD =		\
# 	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)	\
# 	$(ALBERTA_LIBS)				\
# 	$(LDADD)
# MOSTLYCLEANFILES += testlaplacedirichletp12d.vtu

NORMALTESTS += testelasticity
testelasticity_SOURCES = testelasticity.cc
testelasticity_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(SUPERLU_CPPFLAGS)
testelasticity_LDFLAGS = $(AM_LDFLAGS)
testelasticity_LDADD =		          \
	$(LDADD)                          \
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)
MOSTLYCLEANFILES += testelasticity.vtu

NORMALTESTS += testlocalfunctionspace
testlocalfunctionspace_SOURCES = testlocalfunctionspace.cc
testlocalfunctionspace_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testlocalfunctionspace_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA_LDFLAGS)
testlocalfunctionspace_LDADD =			\
        $(LDADD)                                \
	$(ALBERTA_LIBS)

NORMALTESTS += testpk2dinterpolation
testpk2dinterpolation_SOURCES = testpk2dinterpolation.cc
testpk2dinterpolation_CPPFLAGS = $(AM_CPPFLAGS)				\
	$(ALBERTA2D_CPPFLAGS)
testpk2dinterpolation_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA2D_LDFLAGS)
testpk2dinterpolation_LDADD =					\
	$(LDADD)                                \
	$(ALBERTA2D_LIBS)

NORMALTESTS += testpk
testpk_SOURCES = testpk.cc
testpk_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testpk_LDFLAGS = $(AM_LDFLAGS)			\
	$(ALBERTA_LDFLAGS)
testpk_LDADD =					\
	$(LDADD)                                \
	$(ALBERTA_LIBS)
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
	$(LDADD)                                \
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)	\
	$(ALBERTA_LIBS)
MOSTLYCLEANFILES += poisson_*.vtu

NORMALTESTS += testpoisson-old
testpoisson_old_SOURCES = testpoisson.cc
testpoisson_old_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)			\
	$(SUPERLU_CPPFLAGS)			\
	-DOLD_BACKEND
testpoisson_old_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)			\
	$(SUPERLU_LDFLAGS)
testpoisson_old_LDADD =				\
	$(LDADD)                                \
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)	\
	$(ALBERTA_LIBS)
MOSTLYCLEANFILES += poisson_*.vtu

NORMALTESTS += testrt0
testrt0_SOURCES = testrt0.cc
testrt0_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testrt0_LDFLAGS = $(AM_LDFLAGS)			\
	$(ALBERTA_LDFLAGS)
testrt0_LDADD = 				\
	$(LDADD)                                \
	$(ALBERTA_LIBS)
MOSTLYCLEANFILES += testrt0.vtu

NORMALTESTS += testrt02dgridfunctionspace
testrt02dgridfunctionspace_SOURCES = testrt02dgridfunctionspace.cc
testrt02dgridfunctionspace_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testrt02dgridfunctionspace_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA_LDFLAGS)
testrt02dgridfunctionspace_LDADD = 		\
	$(LDADD)                                \
	$(ALBERTA_LIBS)
MOSTLYCLEANFILES += rt02dgridfunctionspace-*.vtu

NORMALTESTS += testutilities
testutilities_SOURCES = testutilities.cc
testutilities_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)
testutilities_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)
testutilities_LDADD = 				\
	$(LDADD)                                \
	$(ALBERTA_LIBS)
MOSTLYCLEANFILES += edger.vtu interpolated.vtu q1.vtu taylorhood.vtu

check_PROGRAMS += testinterpolate
testinterpolate_SOURCES = testinterpolate.cc

if PETSC
check_PROGRAMS += testpoisson-petsc
testpoisson_petsc_SOURCES = testpoisson-petsc.cc
testpoisson_petsc_CPPFLAGS = @AM_CPPFLAGS@ $(PETSC_CPPFLAGS)
testpoisson_petsc_LDFLAGS = $(AM_LDFLAGS) $(PETSC_LDFLAGS)
testpoisson_petsc_LDADD = $(LDADD) $(PETSC_LIBS)
endif

if GMP
check_PROGRAMS += testopbfem
testopbfem_SOURCES = testopbfem.cc
testopbfem_CPPFLAGS = @AM_CPPFLAGS@ $(GMP_CPPFLAGS)
testopbfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopbfem_LDADD = $(LDADD) $(GMP_LIBS)
endif

if HAVE_DUNE_ALUGRID
check_PROGRAMS += testordering
testordering_SOURCES = testordering.cc
endif HAVE_DUNE_ALUGRID

check_PROGRAMS += testpermutedordering
testpermutedordering_SOURCES = testpermutedordering.cc

check_PROGRAMS += testdatahandle
testdatahandle_SOURCES = testdatahandle.cc

check_PROGRAMS += testrtfem
testrtfem_SOURCES = testrtfem.cc

check_PROGRAMS += testbdmfem
testbdmfem_SOURCES = testbdmfem.cc

MOSTLYCLEANFILES += poisson_petsc_*.vtu

if UG
check_PROGRAMS += testnonoverlapping
testnonoverlapping_SOURCES = testnonoverlapping.cc
endif

check_PROGRAMS += testvectoriterator
testvectoriterator_SOURCES = testvectoriterator.cc

check_PROGRAMS += testchunkedblockordering
testchunkedblockordering_SOURCES = testchunkedblockordering.cc

NORMALTESTS += testsimplebackend
testsimplebackend_SOURCES = testsimplebackend.cc
MOSTLYCLEANFILES += simplebackend_*.vtu

if EIGEN

NORMALTESTS += testeigenbackend
testeigenbackend_SOURCES = testeigenbackend.cc
testeigenbackend_CXXFLAGS = @AM_CPPFLAGS@ $(EIGEN_CPPFLAGS) $(EIGEN_CXXFLAGS)
testeigenbackend_LDFLAGS = $(AM_LDFLAGS) $(EIGEN_LDFLAGS)
testeigenbackend_LDADD = $(LDADD) $(EIGEN_LIBS)
MOSTLYCLEANFILES += eigenbackend_*.vtu

endif


if GMP

NORMALTESTS += test-dg-amg
test_dg_amg_SOURCES = test-dg-amg.cc
test_dg_amg_CXXFLAGS = @AM_CPPFLAGS@ $(GMP_CPPFLAGS)
test_dg_amg_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
test_dg_amg_LDADD = $(GMP_LIBS) $(LDADD)

endif

include $(top_srcdir)/am/global-rules

DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"
