dnl DUNE_EIGEN
dnl ------------------------------------------------------
dnl use pkgconfig to check if EIGEN is usable and working
AC_DEFUN([DUNE_EIGEN], [
  AC_REQUIRE([PKG_PROG_PKG_CONFIG])

  PKG_CHECK_MODULES([EIGEN], [eigen3], [
    with_eigen=yes
    AC_DEFINE([HAVE_EIGEN],
      [1],
      [Define wether the eigen includes were found.])
    # we use C++ ... rename CFLAGS to CXXFLAGS
    EIGEN_CXXFLAGS="$EIGEN_CFLAGS"
    EIGEN_CFLAGS=""
    AC_SUBST(EIGEN_CPPFLAGS)
    AC_SUBST(EIGEN_CFLAGS)
    ], [
    AC_MSG_WARN($EIGEN_PKG_ERRORS)
    with_eigen=no
    ])

  DUNE_ADD_SUMMARY_ENTRY([Eigen],["$with_eigen"])

])
