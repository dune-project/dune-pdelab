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

#ifndef DUNE_PDELAB_CHECK_FORMAT_STRINGS
#ifdef NDEBUG
#define DUNE_PDELAB_CHECK_FORMAT_STRINGS 0
#else
/**
 * \brief Flag for enabling compile-time verification of `{fmt}` format strings. This flag controls
 * whether `{fmt}` will parse your format strings at compile time and make sure that all of the
 * referenced formatting arguments are passed to the formatting call.
 *
 * \sa fmt
 */
#define DUNE_PDELAB_CHECK_FORMAT_STRINGS 1
#endif
#endif


#ifdef DOXYGEN
/**
 * \brief Flag for enabling internal consistency checks of the Logger class, in particular whether a
 * Logger is attached to a LoggerBackend. If this flag is enabled, logging a message to a Logger
 * without attached LoggerBackend will throw an exception of type LoggingError; if it is disabled,
 * it results in undefined behavior.
 *
 * \note After including checks.hh, you cannot check whether this flag has been enabled or not.
 *
 * \sa Logger
 */
#define DUNE_PDELAB_CHECK_LOGGER 1
#else
#ifndef DUNE_PDELAB_CHECK_LOGGER
#ifdef NDEBUG
#define DUNE_PDELAB_CHECK_LOGGER(check)
#else
#define DUNE_PDELAB_CHECK_LOGGER(check)                                        \
  if (not (check)) {                                                           \
    DUNE_THROW(LoggingError,                                                   \
               "Error while using Logger: " << #check );                       \
  }
#endif
#endif
#endif

/**
 * \} checks
 */

#endif // DUNE_PDELAB_COMMON_CHECKS_HH
