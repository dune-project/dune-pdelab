# tests for C++0x template alias support
# the associated macro is called HAVE_TEMPLATE_ALIASES

AC_DEFUN([TEMPLATE_ALIASES_CHECK],[
  AC_CACHE_CHECK([whether template aliases are supported], dune_cv_template_aliases_support, [
    AC_REQUIRE([AC_PROG_CXX])
    AC_REQUIRE([GXX0X])
    AC_LANG_PUSH([C++])
    AC_RUN_IFELSE([
      AC_LANG_PROGRAM([

        template<typename T, typename U>
        struct A
        {};

        template<typename T>
        using A1 = A<T,int>;

        template<typename U>
        using A2 = A<int,U>;

        template<typename T, typename U>
        struct assert_equal;

        template<typename T>
        struct assert_equal<T,T>
        {};
        ],
        [
          assert_equal<A1<int>,A2<int> >();
          assert_equal<A<bool,int>,A1<bool> >();
          return 0;
        ])],
      dune_cv_template_aliases_support=yes,
      dune_cv_template_aliases_support=no)
    AC_LANG_POP
  ])
  if test "x$dune_cv_template_aliases_support" = xyes; then
    AC_DEFINE(HAVE_TEMPLATE_ALIASES, 1, [Define to 1 if template aliases are supported])
  fi
])
