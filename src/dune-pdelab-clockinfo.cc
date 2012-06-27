// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#if defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE < 199309L
#undef _POSIX_C_SOURCE
#endif
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <ostream>

#include <unistd.h>

#define stringify(arg) #arg
#define estringify(arg) stringify(arg)

int main() {
#if _POSIX_TIMERS == -1
  std::cout << "_POSIX_TIMERS == -1:" << std::endl
            << "    clock_*() and CLOCK_REALTIME are not supported"
            << std::endl;
#elif _POSIX_TIMERS >= 0
# if _POSIX_TIMERS > 0
  std::cout << "_POSIX_TIMERS > 0:" << std::endl
            << "    clock_*() and CLOCK_REALTIME are supported and functional"
            << std::endl;
# elif defined(_POSIX_TIMERS)
  std::cout << "_POSIX_TIMERS == 0:" << std::endl
            << "    clock_*() and CLOCK_REALTIME are present, but "
            << "functionality is determined via sysconf()" << std::endl;
# else // !defined(_POSIX_TIMERS)
  std::cout << "_POSIX_TIMERS undefined:" << std::endl
            << "    clock_*() and CLOCK_REALTIME are present, but "
            << "functionality is determined via sysconf()" << std::endl;
# endif
  if(sysconf(_SC_TIMERS) > 0)
    std::cout << "sysconf(_SC_TIMERS) > 0:" << std::endl
              << "    clock_*() and CLOCK_REALTIME are functional"
              << std::endl;
  else
    std::cout << "sysconf(_SC_TIMERS) <= 0:" << std::endl
              << "    clock_*() and CLOCK_REALTIME are not functional"
              << std::endl;
#else // _POSIX_TIMERS < -1
  std::cout << "Error: strange value for _POSIX_TIMERS: "
            << estringify(_POSIX_TIMERS) << std::endl;
#endif

  std::cout << std::endl;

#if _POSIX_CPUTIME == -1
  std::cout << "_POSIX_CPUTIME == -1:" << std::endl
            << "    CLOCK_PROCESS_CPUTIME_ID is not supported" << std::endl;
#elif _POSIX_CPUTIME >= 0
# if _POSIX_CPUTIME > 0
  std::cout << "_POSIX_CPUTIME > 0:" << std::endl
            << "    CLOCK_PROCESS_CPUTIME_ID is supported and functional"
            << std::endl;
# elif defined(_POSIX_CPUTIME)
  std::cout << "_POSIX_CPUTIME == 0:" << std::endl
            << "    CLOCK_PROCESS_CPUTIME_ID is present, but functionality is "
            << "determined via sysconf()" << std::endl;
# else // !defined(_POSIX_CPUTIME)
  std::cout << "_POSIX_CPUTIME undefined:" << std::endl
            << "    CLOCK_PROCESS_CPUTIME_ID is present, but functionality is "
            << "determined via sysconf()" << std::endl;
# endif
  if(sysconf(_SC_CPUTIME) > 0)
    std::cout << "sysconf(_SC_CPUTIME) > 0:" << std::endl
              << "    CLOCK_PROCESS_CPUTIME_ID is functional" << std::endl;
  else
    std::cout << "sysconf(_SC_CPUTIME) <= 0:" << std::endl
              << "    CLOCK_PROCESS_CPUTIME_ID is not functional" << std::endl;
#else // _POSIX_CPUTIME < -1
  std::cout << "Error: strange value for _POSIX_CPUTIME: "
            << estringify(_POSIX_CPUTIME) << std::endl;
#endif
}
