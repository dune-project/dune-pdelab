# Additional checks needed to build the module
AC_DEFUN([DUNE_PDELAB_CHECKS])

# Additional checks needed to find the module
AC_DEFUN([DUNE_PDELAB_CHECK_MODULE],[
  AC_MSG_NOTICE([Searching for dune-pdelab...])
  DUNE_CHECK_MODULES([dune-pdelab], [pdelab/common/function.hh])
])
