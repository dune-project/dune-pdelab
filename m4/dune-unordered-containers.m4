## -*- autoconf -*-
AC_DEFUN([DUNE_UNORDERED_HEADERS], [
  # Only run the check if we also run the basic tr1 test from dune-common.
  AC_LANG_PUSH([C++])
  AS_IF([test "x$enable_tr1_headers" != "xno"],
    [AC_CHECK_HEADERS([unordered_map tr1/unordered_map unordered_set tr1/unordered_set])
  ])
  AC_LANG_POP([C++])
])
