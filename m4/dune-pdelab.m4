# Additional checks needed to build the module
AC_DEFUN([DUNE_PDELAB_CHECKS],
[
  AC_REQUIRE([RVALUE_REFERENCES_CHECK])
  AC_REQUIRE([VARIADIC_TEMPLATES_CHECK])
  AC_REQUIRE([VARIADIC_CONSTRUCTOR_SFINAE_CHECK])
  AC_REQUIRE([DUNE_PATH_PETSC])
])

# Additional checks needed to find the module
AC_DEFUN([DUNE_PDELAB_CHECK_MODULE],[
  AC_MSG_NOTICE([Searching for dune-pdelab...])
  DUNE_CHECK_MODULES([dune-pdelab], [pdelab/common/function.hh])
])
