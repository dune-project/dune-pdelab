# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([dune/pdelab/test/testopbfem.cc])
AC_CONFIG_HEADERS([config.h])

# add configure flags needed to create log files for dune-autobuild
DUNE_AUTOBUILD_FLAGS
#
# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid dune-localfunctions dune-istl
# this implies checking for [dune-common], [dune-grid], [dune-localfunctions], [dune-istl]
DUNE_CHECK_ALL

# preset variable to path such that #include <dune/...> works
#AC_SUBST([DUNE_PDELAB_ROOT], '$(top_builddir)')
#AC_SUBST([AM_CPPFLAGS], '-I$(top_srcdir)')
#AC_SUBST([LOCAL_LIBS], '$(top_builddir)/lib/libpdelab.la')

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
AC_SUBST([LDADD], '$(top_builddir)/lib/libdunepdelab.la $(DUNE_LIBS)')

AC_CONFIG_FILES([
  cmake/Makefile
  cmake/modules/Makefile
  doc/doxygen/Doxyfile
  doc/doxygen/Makefile
  doc/Makefile
  dune-pdelab.pc
  dune/Makefile
  dune/pdelab/adaptivity/Makefile
  dune/pdelab/backend/Makefile
  dune/pdelab/backend/common/Makefile
  dune/pdelab/backend/eigen/Makefile
  dune/pdelab/backend/istl/Makefile
  dune/pdelab/backend/simple/Makefile
  dune/pdelab/boilerplate/Makefile
  dune/pdelab/common/Makefile
  dune/pdelab/constraints/Makefile
  dune/pdelab/constraints/common/Makefile
  dune/pdelab/finiteelement/Makefile
  dune/pdelab/finiteelementmap/Makefile
  dune/pdelab/function/Makefile
  dune/pdelab/gridfunctionspace/Makefile
  dune/pdelab/gridoperator/Makefile
  dune/pdelab/gridoperator/common/Makefile
  dune/pdelab/gridoperator/default/Makefile
  dune/pdelab/gridoperator/onestep/Makefile
  dune/pdelab/gridoperator/tbb/Makefile
  dune/pdelab/instationary/Makefile
  dune/pdelab/linearsolver/Makefile
  dune/pdelab/localoperator/Makefile
  dune/pdelab/Makefile
  dune/pdelab/multistep/Makefile
  dune/pdelab/newton/Makefile
  dune/pdelab/ordering/Makefile
  dune/pdelab/stationary/Makefile
  dune/pdelab/test/Makefile
  dune/pdelab/test/fem/Makefile
  lib/Makefile
  m4/Makefile
  Makefile
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
