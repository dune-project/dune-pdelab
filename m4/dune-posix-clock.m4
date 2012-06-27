dnl DUNE_FUNC_POSIX_CLOCK
dnl ------------------------------------------------------
dnl Check whether the clock_*()-functions are available.  The result is
dnl recorded as follows:
dnl
dnl shell variables:
dnl   dune_cv_defs_posix_clock
dnl     "yes" or "no"
dnl   dune_cv_lib_posix_clock
dnl     linker options to get the required library
dnl     "no" if no suitable library was found
dnl   dune_cv_func_posix_clock
dnl     "yes" or "no"
dnl
dnl defines:
dnl   HAVE_POSIX_CLOCK
dnl     undef or ENABLE_POSIX_CLOCK, equivalent to _POSIX_TIMERS >= 0
dnl
dnl automake conditionals:
dnl   POSIX_CLOCK
dnl
dnl Makefile variables:
dnl   POSIX_CLOCK_CPPFLAGS
dnl     -DENABLE_POSIX_CLOCK
dnl   POSIX_CLOCK_LDFLAGS
dnl   POSIX_CLOCK_LIBS
dnl
dnl NOTE: This only tests for the presence of the clock_*() functions and for
dnl       the flags necessary to compile and link programs that reference them.
dnl       It does not test whether these functions are actually functional.
dnl       You have to check _POSIX_TIMERS and maybe sysconf(_SC_TIMERS) to
dnl       figure that out.
dnl
dnl NOTE: This will not make these functions available automatically, you still
dnl       have to make sure _POSIX_C_SOURCE is defined >= 199309L at the top of
dnl       your compilation unit before including "config.h".
AC_DEFUN([DUNE_FUNC_POSIX_CLOCK], [
  AC_LANG_PUSH([C])

  AC_CACHE_CHECK(
    [whether clock_gettime()'s definition is available],
    [dune_cv_defs_posix_clock],
    [AC_COMPILE_IFELSE([_DUNE_POSIX_CLOCK_TESTPROG],
      [dune_cv_defs_posix_clock=yes],
      [dune_cv_defs_posix_clock=no])])
  AS_CASE(["$dune_cv_defs_posix_clock"],
    [yes], [
      AC_CACHE_CHECK(
        [for library required for clock_gettime()],
        [dune_cv_lib_posix_clock],
        [
          dune_cv_lib_posix_clock=no
          _DUNE_POSIX_CLOCK_CHECK_LIB([""])
          _DUNE_POSIX_CLOCK_CHECK_LIB([-lrt])
        ])
     ])
  AC_MSG_CHECKING([for clock_gettime()])
  AS_CASE(["$dune_cv_defs_posix_clock"],
    [yes], [
      AS_CASE(["$dune_cv_lib_posix_clock"],
        [no], [
          dune_cv_func_posix_clock=no
          AC_MSG_RESULT([no])],
        [""], [
          dune_cv_func_posix_clock=yes
          AC_MSG_RESULT([yes (no lib required)])],
        [
          dune_cv_func_posix_clock=yes
          AC_MSG_RESULT([yes (in $dune_cv_lib_posix_clock)])])
     ],
     [
       dune_cv_func_posix_clock=no
       AC_MSG_RESULT([no])
     ])
  AS_CASE(["$dune_cv_func_posix_clock"],
    [yes], [
      AC_DEFINE([HAVE_POSIX_CLOCK], [ENABLE_POSIX_CLOCK],
        [Define if clock_gettime() is available, given _POSIX_C_SOURCE >= 199309L])

      AC_SUBST([POSIX_CLOCK_CPPFLAGS], [-DENABLE_POSIX_CLOCK])
      AC_SUBST([POSIX_CLOCK_LDFLAGS],  [])
      AC_SUBST([POSIX_CLOCK_LIBS],     ["$dune_cv_lib_posix_clock"])
    ],
    [
      AC_SUBST([POSIX_CLOCK_CPPFLAGS], [])
      AC_SUBST([POSIX_CLOCK_LDFLAGS],  [])
      AC_SUBST([POSIX_CLOCK_LIBS], [])
    ])
  AM_CONDITIONAL([POSIX_CLOCK], [test x"$dune_cv_func_posix_clock" = xyes])

  AC_LANG_POP([C])
])

AC_DEFUN([_DUNE_POSIX_CLOCK_CHECK_LIB], [
  AS_CASE(["$dune_cv_lib_posix_clock"],
    [no], [
      dune_save_LIBS="$LIBS"
      LIBS=$1" $LIBS"
      AC_LINK_IFELSE([_DUNE_POSIX_CLOCK_TESTPROG],
        [dune_cv_lib_posix_clock=$1],
        [dune_cv_lib_posix_clock=no])
      LIBS="$dune_save_LIBS"
    ])
])

dnl Generate test
AC_DEFUN([_DUNE_POSIX_CLOCK_TESTPROG], [dnl
    AC_LANG_PROGRAM(
[[#if defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE < 199309L
#undef _POSIX_C_SOURCE
#endif
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif

#include <time.h>
#include <unistd.h>

#if _POSIX_TIMERS < 0
#error _POSIX_TIMERS < 0
#endif
]],
[[struct timespec result;
clock_gettime(CLOCK_REALTIME, &result);
]])])
