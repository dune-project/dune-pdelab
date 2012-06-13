// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

//make sure clock_gettime is available
#if defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE < 199309L
#undef _POSIX_C_SOURCE
#endif
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <cerrno>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>

#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>

#include <dune/common/exceptions.hh>

#include "clock.hh"

namespace Dune {
  namespace PDELab {

    std::ostream &operator<<(std::ostream &s, const TimeSpec &t) {
      std::ostringstream tmp;
      tmp << t.tv_sec << '.' << std::setfill('0') << std::setw(9) << t.tv_nsec;
      std::string tmpstr = tmp.str();
      if(s.precision() < 9)
        tmpstr.resize(tmpstr.size() - ( 9 - s.precision() ));
      if(s.precision() == 0)
        tmpstr.resize(tmpstr.size() - 1);
      s << tmpstr;
      return s;
    }

    TimeSpec getWallTime() {
      timespec result;
      if(clock_gettime(CLOCK_REALTIME, &result) < 0)
        DUNE_THROW(ClockError, "clock_gettime(CLOCK_REALTIME, ...) failed: "
                   "errno = " << errno);
      TimeSpec tmp = { result.tv_sec, result.tv_nsec };
      return tmp;
    }

    TimeSpec getProcessTime() {
      // Use clock_gettime(CLOCK_PROCESS_CPUTIME_ID, ...) even though that may
      // be problematic in the context of process migration between cores.  In
      // practice, it appears to still be far better then the next best
      // alternative, getrusage(), which will only update the clock every
      // jiffy.
      timespec result;
      if(clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &result) < 0)
        DUNE_THROW(ClockError, "clock_gettime(CLOCK_PROCESS_CPUTIME_ID, ...) "
                   "failed: errno = " << errno);
      TimeSpec tmp = { result.tv_sec, result.tv_nsec };
      return tmp;
    }

  } // namespace PDELab
} // namespace Dune
