# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([dune/pdelab/test/testcountingptr.cc])
AM_CONFIG_HEADER([config.h])

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
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  m4/Makefile
  dune/Makefile
  dune/pdelab/Makefile
  dune/pdelab/test/Makefile
  dune/pdelab/common/Makefile
  dune/pdelab/finiteelementmap/Makefile
  dune/pdelab/gridfunctionspace/Makefile
  dune/pdelab/gridoperatorspace/Makefile
  dune/pdelab/backend/Makefile
  dune/pdelab/localoperator/Makefile
  dune-pdelab.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
