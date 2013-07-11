# tests for C++11 decltype support
# the associated macro is called HAVE_STD_DECLTYPE

AC_DEFUN([CXX11_DECLTYPE_CHECK],[
  AC_CACHE_CHECK([whether the keyword decltype is supported], dune_cv_cxx11_decltype_support, [
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

        struct A;

        const A& foo();

        ],
        [
          return check_equal<decltype(foo()),const A&>::result;
        ])],
      dune_cv_cxx11_decltype_support=yes,
      dune_cv_cxx11_decltype_support=no)
    AC_LANG_POP
  ])
  if test "x$dune_cv_cxx11_decltype_support" = xyes; then
    AC_DEFINE(HAVE_STD_DECLTYPE, 1, [Define to 1 if the C++11 keyword decltype is supported])
  fi
])
