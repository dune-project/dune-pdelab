#ifndef DUNE_PDELAB_COMMON_CHECKS_HH
#define DUNE_PDELAB_COMMON_CHECKS_HH

#include <dune/pdelab/common/exceptions.hh>

/**
 * \addtogroup checks Optional code correctness checks
 * \brief Flags for enabling or disabling certain potentially expensive assertions in PDELab.
 * \{
 */

/** \file
 * Handling of PDELab-specific correctness verification flags.
 *
 * This file contains processing for a number of PDELab-specific verification flags that check the
 * correctness of your code with respect to different problems.
 *
 * All of these flags are enabled by default if your code is compiled in debug mode, i.e. when the
 * macro NDEBUG is not defined. If you want to override this default for a check, just set its value
 * before including this file.
 */

#ifndef DUNE_PDELAB_ENABLE_CHECK_ASSEMBLY
#ifndef NDEBUG
#define DUNE_PDELAB_ENABLE_CHECK_ASSEMBLY 1
#else
#define DUNE_PDELAB_ENABLE_CHECK_ASSEMBLY 0
#endif
#endif

/**
 * \} checks
 */

#endif // DUNE_PDELAB_COMMON_CHECKS_HH
