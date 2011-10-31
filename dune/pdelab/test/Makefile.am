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

NORMALTESTS += testcountingptr
testcountingptr_SOURCES = testcountingptr.cc
testcountingptr_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA_CPPFLAGS)
testcountingptr_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA_LDFLAGS)
testcountingptr_LDADD =				\
	$(ALBERTA_LIBS)				\
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

NORMALTESTS += testmultistepcached
testmultistepcached_SOURCES = testmultistepcached.cc
testmultistepcached_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(DUNEMPICPPFLAGS)
testmultistepcached_LDFLAGS = $(AM_LDFLAGS)	\
	$(DUNEMPILDFLAGS)
testmultistepcached_LDADD =			\
	$(DUNEMPILIBS)				\
	$(LDADD)
MOSTLYCLEANFILES +=					\
	s*:p*:testmultistepcached_yasp_P1_1d-*.vtp	\
	s*:testmultistepcached_yasp_P1_1d-*.pvtp	\
	testmultistepcached_yasp_P1_1d.pvd

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

TESTS += testvectorwave.conf
check_SCRIPTS += testvectorwave.conf
check_PROGRAMS += testvectorwave
testvectorwave_SOURCES = testvectorwave.cc
MOSTLYCLEANFILES +=				\
	s*:p*:vectorwave_*.vtu			\
	s*:vectorwave_*.pvtu			\
	vectorwave_*.pvd

TESTS += testvectorwavecached.conf
check_SCRIPTS += testvectorwavecached.conf
check_PROGRAMS += testvectorwavecached
testvectorwavecached_SOURCES = testvectorwavecached.cc
MOSTLYCLEANFILES +=				\
	s*:p*:vectorwavecached_*.vtu		\
	s*:vectorwavecached_*.pvtu		\
	vectorwavecached_*.pvd

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

MOSTLYCLEANFILES += poisson_petsc_*.vtu

include $(top_srcdir)/am/global-rules

DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"
