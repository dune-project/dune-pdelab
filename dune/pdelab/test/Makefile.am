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
	l2norm.hh                               \
        typetreetargetnodes.hh                  \
        typetreetestswitch.hh                   \
        typetreetestutility.hh


noinst_SCRIPTS =				\
	make_pvd.sh

headercheck_IGNORE =				\
        typetreetargetnodes.hh			\
        typetreetestswitch.hh			\
        typetreetestutility.hh

GRIDDIM = 2

# scripts are not distributed automatically
# but since this is need for make all:
EXTRA_DIST = CMakeLists.txt make_pvd.sh

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

NORMALTESTS += testclock
testclock_SOURCES = testclock.cc
# don't include all the grid stuff
testclock_CPPFLAGS = @AM_CPPFLAGS@
testclock_LDFLAGS = @AM_LDFLAGS@
testclock_LDADD = @LDADD@

if ALUGRID
NORMALTESTS += test-composed-iis-gfs
test_composed_iis_gfs_SOURCES = test-composed-iis-gfs.cc
endif ALUGRID

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
MOSTLYCLEANFILES += testlaplacedirichletccfv-2d.vtu \
	testlaplacedirichletccfv-3d.vtu

NORMALTESTS += testlaplacedirichletccfv-simplified
testlaplacedirichletccfv_simplified_SOURCES = $(testlaplacedirichletccfv_SOURCES)
testlaplacedirichletccfv_simplified_CPPFLAGS = $(testlaplacedirichletccfv_CPPFLAGS)	\
	-DTEST_SIMPLIFIED_INFRASTRUCTURE
testlaplacedirichletccfv_simplified_LDFLAGS = $(testlaplacedirichletccfv_LDFLAGS)
testlaplacedirichletccfv_simplified_LDADD = $(testlaplacedirichletccfv_LDADD)
MOSTLYCLEANFILES += testlaplacedirichletccfv-simplified-2d.vtu \
	testlaplacedirichletccfv-simplified-3d.vtu

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

NORMALTESTS += testelasticity
testelasticity_SOURCES = testelasticity.cc
testelasticity_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(SUPERLU_CPPFLAGS)
testelasticity_LDFLAGS = $(AM_LDFLAGS)	\
	$(SUPERLU_LDFLAGS)
testelasticity_LDADD =		\
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)	\
	$(LDADD)
MOSTLYCLEANFILES += testelasticity.vtu

NORMALTESTS += testlocalfunctionspace
testlocalfunctionspace_SOURCES = testlocalfunctionspace.cc
testlocalfunctionspace_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testlocalfunctionspace_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA_LDFLAGS)
testlocalfunctionspace_LDADD =			\
	$(ALBERTA_LIBS)				\
	$(LDADD)

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

check_PROGRAMS += testtypetree testtypetree-no-sfinae testtypetree-no-variadic testtypetree-no-rvalue-refs testtypetree-legacy
testtypetree_SOURCES = testtypetree.cc
testtypetree_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_TYPETREE
testtypetree_no_sfinae_SOURCES = testtypetree.cc
testtypetree_no_sfinae_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_TYPETREE_NO_SFINAE
testtypetree_no_variadic_SOURCES = testtypetree.cc
testtypetree_no_variadic_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_TYPETREE_NO_VARIADIC
testtypetree_no_rvalue_refs_SOURCES = testtypetree.cc
testtypetree_no_rvalue_refs_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_TYPETREE_NO_RVALUE_REFERENCES
testtypetree_legacy_SOURCES = testtypetree.cc
testtypetree_legacy_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_TYPETREE_LEGACY

check_PROGRAMS += testtypetreetransformation
testtypetreetransformation_SOURCES = testtypetreetransformation.cc
testtypetreetransformation_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_TYPETREE

check_PROGRAMS += testinterpolate
testinterpolate_SOURCES = testinterpolate.cc

check_PROGRAMS += testfilteredcompositenode
testfilteredcompositenode_SOURCES = testfilteredcompositenode.cc
testfilteredcompositenode_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_TYPETREE

check_PROGRAMS += testproxynode
testproxynode_SOURCES = testproxynode.cc
testproxynode_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_TYPETREE

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

check_PROGRAMS += testordering
testordering_SOURCES = testordering.cc

check_PROGRAMS += testdatahandle
testdatahandle_SOURCES = testdatahandle.cc

check_PROGRAMS += testrtfem
testrtfem_SOURCES = testrtfem.cc

check_PROGRAMS += testbdmfem
testbdmfem_SOURCES = testbdmfem.cc

MOSTLYCLEANFILES += poisson_petsc_*.vtu

check_PROGRAMS += testnonoverlapping
testnonoverlapping_SOURCES = testnonoverlapping.cc

check_PROGRAMS += testvectoriterator
testvectoriterator_SOURCES = testvectoriterator.cc

check_PROGRAMS += testvectoriterator-no-variadic-templates
testvectoriterator_no_variadic_templates_SOURCES = testvectoriterator.cc
testvectoriterator_no_variadic_templates_CPPFLAGS = $(AM_CPPFLAGS) -DDISABLE_VARIADIC_TEMPLATES

NORMALTESTS += testdensebackend
testdensebackend_SOURCES = testdensebackend.cc
testdensebackend_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_CPPFLAGS)			\
	$(SUPERLU_CPPFLAGS)
testdensebackend_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)			\
	$(SUPERLU_LDFLAGS)
testdensebackend_LDADD =				\
	$(SUPERLU_LDFLAGS) $(SUPERLU_LIBS)	\
	$(ALBERTA_LIBS)				\
	$(LDADD)
MOSTLYCLEANFILES += densebackend_*.vtu


include $(top_srcdir)/am/global-rules

DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"
