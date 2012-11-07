dnl DUNE_EIGEN
dnl ------------------------------------------------------
dnl use pkgconfig to check if EIGEN is usable and working
AC_DEFUN([DUNE_EIGEN], [
  AC_REQUIRE([PKG_PROG_PKG_CONFIG])

  # use pkg-config to search for eigen
  PKG_CHECK_MODULES([EIGEN], [eigen3], [
    HAVE_EIGEN=yes
    # we use C++ ... rename CFLAGS to CXXFLAGS
    EIGEN_CXXFLAGS="$EIGEN_CFLAGS -DENABLE_EIGEN=1"
    EIGEN_CFLAGS=""
    AC_SUBST(EIGEN_CXXFLAGS)
    AC_SUBST(EIGEN_CFLAGS)
    ], [
    AC_MSG_WARN($EIGEN_PKG_ERRORS)
    HAVE_EIGEN=no
    ])

  # use the ENABLE_* trick
  AC_DEFINE([HAVE_EIGEN], ENABLE_EIGEN,
    [This is only true if EIGEN was found by configure 
     _and_ if the application uses the EIGEN_CXXFLAGS])

  # add automake conditionsl
  AM_CONDITIONAL(EIGEN, test x$HAVE_EIGEN = xyes)

  # print summary
  DUNE_ADD_SUMMARY_ENTRY([Eigen],["$HAVE_EIGEN"])

])
