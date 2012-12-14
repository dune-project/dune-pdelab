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
	$(DUNE_CPPFLAGS)				\
	$(DUNEMPICPPFLAGS)

AM_LDFLAGS = @AM_LDFLAGS@			\
	$(UG_LDFLAGS)				\
	$(ALUGRID_LDFLAGS)			\
	$(DUNE_LDFLAGS)				\
	$(DUNEMPILDFLAGS)

LDADD =						\
	$(ALUGRID_LIBS)				\
	$(DUNE_LIBS)				\
	$(DUNEMPILIBS)				\
	$(UG_LIBS)				\
	@LDADD@


check_PROGRAMS += testrt0simplex2dfem
testrt0simplex2dfem_SOURCES = testfemwithgfs.cc
testrt0simplex2dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=2 -DFEM_SIMPLEX -DUSE_RT_BDM_FEM_FACTORY -DFEM_FACTORY=RTBDMFEMFactory -DFEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap -DFEM_FACTORY_ORDER=0
TESTS += testrt0simplex2dfem

check_PROGRAMS += testrt1simplex2dfem
testrt1simplex2dfem_SOURCES = testfemwithgfs.cc
testrt1simplex2dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=2 -DFEM_SIMPLEX -DUSE_RT_BDM_FEM_FACTORY -DFEM_FACTORY=RTBDMFEMFactory -DFEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap -DFEM_FACTORY_ORDER=1
TESTS += testrt1simplex2dfem

check_PROGRAMS += testrt0cube2dfem
testrt0cube2dfem_SOURCES = testfemwithgfs.cc
testrt0cube2dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=2 -DFEM_CUBE -DUSE_RT_BDM_FEM_FACTORY -DFEM_FACTORY=RTBDMFEMFactory -DFEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap -DFEM_FACTORY_ORDER=0
TESTS += testrt0cube2dfem

check_PROGRAMS += testrt1cube2dfem
testrt1cube2dfem_SOURCES = testfemwithgfs.cc
testrt1cube2dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=2 -DFEM_CUBE -DUSE_RT_BDM_FEM_FACTORY -DFEM_FACTORY=RTBDMFEMFactory -DFEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap -DFEM_FACTORY_ORDER=1
TESTS += testrt1cube2dfem

check_PROGRAMS += testrt2cube2dfem
testrt2cube2dfem_SOURCES = testfemwithgfs.cc
testrt2cube2dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=2 -DFEM_CUBE -DUSE_RT_BDM_FEM_FACTORY -DFEM_FACTORY=RTBDMFEMFactory -DFEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap -DFEM_FACTORY_ORDER=2
TESTS += testrt2cube2dfem

check_PROGRAMS += testrt0cube3dfem
testrt0cube3dfem_SOURCES = testfemwithgfs.cc
testrt0cube3dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=3 -DFEM_CUBE -DUSE_RT_BDM_FEM_FACTORY -DFEM_FACTORY=RTBDMFEMFactory -DFEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap -DFEM_FACTORY_ORDER=0
TESTS += testrt0cube3dfem

check_PROGRAMS += testrt1cube3dfem
testrt1cube3dfem_SOURCES = testfemwithgfs.cc
testrt1cube3dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=3 -DFEM_CUBE -DUSE_RT_BDM_FEM_FACTORY -DFEM_FACTORY=RTBDMFEMFactory -DFEM_FACTORY_FEM_CLASS=Dune::PDELab::RaviartThomasLocalFiniteElementMap -DFEM_FACTORY_ORDER=1
TESTS += testrt1cube3dfem


check_PROGRAMS += testbdm1simplex2dfem
testbdm1simplex2dfem_SOURCES = testfemwithgfs.cc
testbdm1simplex2dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=2 -DFEM_SIMPLEX -DUSE_RT_BDM_FEM_FACTORY -DFEM_FACTORY=RTBDMFEMFactory -DFEM_FACTORY_FEM_CLASS=Dune::PDELab::BrezziDouglasMariniLocalFiniteElementMap -DFEM_FACTORY_ORDER=1
TESTS += testbdm1simplex2dfem


check_PROGRAMS += testbdm1cube2dfem
testbdm1cube2dfem_SOURCES = testfemwithgfs.cc
testbdm1cube2dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=2 -DFEM_CUBE -DUSE_RT_BDM_FEM_FACTORY -DFEM_FACTORY=RTBDMFEMFactory -DFEM_FACTORY_FEM_CLASS=Dune::PDELab::BrezziDouglasMariniLocalFiniteElementMap -DFEM_FACTORY_ORDER=1
TESTS += testbdm1cube2dfem


check_PROGRAMS += testopb0cube2dfem
testopb0cube2dfem_SOURCES = testfemwithgfs.cc
testopb0cube2dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=2 -DFEM_CUBE -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=0
testopb0cube2dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb0cube2dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb0cube2dfem

check_PROGRAMS += testopb1cube2dfem
testopb1cube2dfem_SOURCES = testfemwithgfs.cc
testopb1cube2dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=2 -DFEM_CUBE -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=1
testopb1cube2dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb1cube2dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb0cube2dfem

check_PROGRAMS += testopb2cube2dfem
testopb2cube2dfem_SOURCES = testfemwithgfs.cc
testopb2cube2dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=2 -DFEM_CUBE -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=2
testopb2cube2dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb2cube2dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb0cube2dfem


check_PROGRAMS += testopb3cube2dfem
testopb3cube2dfem_SOURCES = testfemwithgfs.cc
testopb3cube2dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=2 -DFEM_CUBE -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=2
testopb3cube2dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb3cube2dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb0cube2dfem

check_PROGRAMS += testopb0simplex2dfem
testopb0simplex2dfem_SOURCES = testfemwithgfs.cc
testopb0simplex2dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=2 -DFEM_SIMPLEX -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=0
testopb0simplex2dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb0simplex2dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb0simplex2dfem

check_PROGRAMS += testopb1simplex2dfem
testopb1simplex2dfem_SOURCES = testfemwithgfs.cc
testopb1simplex2dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=2 -DFEM_SIMPLEX -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=1
testopb1simplex2dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb1simplex2dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb1simplex2dfem

check_PROGRAMS += testopb2simplex2dfem
testopb2simplex2dfem_SOURCES = testfemwithgfs.cc
testopb2simplex2dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=2 -DFEM_SIMPLEX -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=2
testopb2simplex2dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb2simplex2dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb2simplex2dfem

check_PROGRAMS += testopb3simplex2dfem
testopb3simplex2dfem_SOURCES = testfemwithgfs.cc
testopb3simplex2dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=2 -DFEM_SIMPLEX -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=3
testopb3simplex2dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb3simplex2dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb3simplex2dfem



check_PROGRAMS += testopb0cube3dfem
testopb0cube3dfem_SOURCES = testfemwithgfs.cc
testopb0cube3dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=3 -DFEM_CUBE -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=0
testopb0cube3dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb0cube3dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb0cube3dfem

check_PROGRAMS += testopb1cube3dfem
testopb1cube3dfem_SOURCES = testfemwithgfs.cc
testopb1cube3dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=3 -DFEM_CUBE -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=1
testopb1cube3dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb1cube3dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb1cube3dfem

check_PROGRAMS += testopb2cube3dfem
testopb2cube3dfem_SOURCES = testfemwithgfs.cc
testopb2cube3dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=3 -DFEM_CUBE -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=2
testopb2cube3dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb2cube3dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb2cube3dfem

check_PROGRAMS += testopb3cube3dfem
testopb3cube3dfem_SOURCES = testfemwithgfs.cc
testopb3cube3dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=3 -DFEM_CUBE -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=2
testopb3cube3dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb3cube3dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb3cube3dfem


check_PROGRAMS += testopb0simplex3dfem
testopb0simplex3dfem_SOURCES = testfemwithgfs.cc
testopb0simplex3dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=3 -DFEM_SIMPLEX -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=0
testopb0simplex3dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb0simplex3dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb0simplex3dfem

check_PROGRAMS += testopb1simplex3dfem
testopb1simplex3dfem_SOURCES = testfemwithgfs.cc
testopb1simplex3dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=3 -DFEM_SIMPLEX -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=1
testopb1simplex3dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb1simplex3dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb1simplex3dfem

check_PROGRAMS += testopb2simplex3dfem
testopb2simplex3dfem_SOURCES = testfemwithgfs.cc
testopb2simplex3dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=3 -DFEM_SIMPLEX -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=2
testopb2simplex3dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb2simplex3dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb2simplex3dfem

check_PROGRAMS += testopb3simplex3dfem
testopb3simplex3dfem_SOURCES = testfemwithgfs.cc
testopb3simplex3dfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) -DFEM_DIM=3 -DFEM_SIMPLEX -DUSE_OPB_FEM_FACTORY -DFEM_FACTORY=OPBFEMFactory -DFEM_FACTORY_ORDER=2
testopb3simplex3dfem_LDFLAGS = $(AM_LDFLAGS) $(GMP_LDFLAGS)
testopb3simplex3dfem_LDADD = $(GMP_LIBS) $(LDADD)
TESTS += testopb3simplex3dfem


check_PROGRAMS += testp0simplex2dfem
testp0simplex2dfem_SOURCES = testfemwithgfs.cc
testp0simplex2dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=2 -DFEM_SIMPLEX -DUSE_PK_FEM_FACTORY -DFEM_FACTORY=PKFEMFactory -DFEM_FACTORY_ORDER=0
TESTS += testp0simplex2dfem

check_PROGRAMS += testp1simplex2dfem
testp1simplex2dfem_SOURCES = testfemwithgfs.cc
testp1simplex2dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=2 -DFEM_SIMPLEX -DUSE_PK_FEM_FACTORY -DFEM_FACTORY=PKFEMFactory -DFEM_FACTORY_ORDER=1
TESTS += testp1simplex2dfem

check_PROGRAMS += testp2simplex2dfem
testp2simplex2dfem_SOURCES = testfemwithgfs.cc
testp2simplex2dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=2 -DFEM_SIMPLEX -DUSE_PK_FEM_FACTORY -DFEM_FACTORY=PKFEMFactory -DFEM_FACTORY_ORDER=2
TESTS += testp2simplex2dfem

check_PROGRAMS += testp3simplex2dfem
testp3simplex2dfem_SOURCES = testfemwithgfs.cc
testp3simplex2dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=2 -DFEM_SIMPLEX -DUSE_PK_FEM_FACTORY -DFEM_FACTORY=PKFEMFactory -DFEM_FACTORY_ORDER=3
TESTS += testp3simplex2dfem


check_PROGRAMS += testp0simplex3dfem
testp0simplex3dfem_SOURCES = testfemwithgfs.cc
testp0simplex3dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=3 -DFEM_SIMPLEX -DUSE_PK_FEM_FACTORY -DFEM_FACTORY=PKFEMFactory -DFEM_FACTORY_ORDER=0
TESTS += testp0simplex3dfem

check_PROGRAMS += testp1simplex3dfem
testp1simplex3dfem_SOURCES = testfemwithgfs.cc
testp1simplex3dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=3 -DFEM_SIMPLEX -DUSE_PK_FEM_FACTORY -DFEM_FACTORY=PKFEMFactory -DFEM_FACTORY_ORDER=1
TESTS += testp1simplex3dfem

check_PROGRAMS += testp2simplex3dfem
testp2simplex3dfem_SOURCES = testfemwithgfs.cc
testp2simplex3dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=3 -DFEM_SIMPLEX -DUSE_PK_FEM_FACTORY -DFEM_FACTORY=PKFEMFactory -DFEM_FACTORY_ORDER=2
TESTS += testp2simplex3dfem

check_PROGRAMS += testp3simplex3dfem
testp3simplex3dfem_SOURCES = testfemwithgfs.cc
testp3simplex3dfem_CPPFLAGS = $(AM_CPPFLAGS) -DFEM_DIM=3 -DFEM_SIMPLEX -DUSE_PK_FEM_FACTORY -DFEM_FACTORY=PKFEMFactory -DFEM_FACTORY_ORDER=3
TESTS += testp3simplex3dfem


include $(top_srcdir)/am/global-rules

DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"
