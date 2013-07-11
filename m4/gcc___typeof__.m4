# tests for the __typeof__ keyword that provides decltype-like functionality
# for older versions of GCC
# the associated macro is called HAVE_GCC___TYPEOF__

AC_DEFUN([GCC___TYPEOF___CHECK],[
  AC_CACHE_CHECK([whether the GCC extension __typeof__ is supported as a decltype fallback], dune_cv_gcc___typeof___support, [
    AC_REQUIRE([AC_PROG_CXX])
    AC_REQUIRE([GXX0X])
    AC_LANG_PUSH([C++])
    AC_RUN_IFELSE([
      AC_LANG_PROGRAM([

        template<typename A, typename B>
        struct check_equal;

        template<typename A>
        struct check_equal<A,A>
        {
          static const int result = 0;
        };

        struct A {};

        A foo();

        ],
        [
          return check_equal<__typeof__(foo()),A>::result;
        ])],
      dune_cv_gcc___typeof___support=yes,
      dune_cv_gcc___typeof___support=no)
    AC_LANG_POP
  ])
  if test "x$dune_cv_gcc___typeof___support" = xyes; then
    AC_DEFINE(HAVE_GCC___TYPEOF__, 1, [Define to 1 if the GCC extension __typeof__ (decltype fallback) is supported])
  fi
])
