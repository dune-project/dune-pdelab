// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

//make sure clock_gettime is available
#if defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE < 199309L
#undef _POSIX_C_SOURCE
#endif
#if !defined(_POSIX_C_SOURCE) && !defined(__APPLE__)
#define _POSIX_C_SOURCE 199309L
#endif

//make sure gettimeofday is available
#ifndef _BSD_SOURCE
#define _BSD_SOURCE
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
#include <unistd.h>

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

    //////////////////////////////////////////////////////////////////////
    //
    //  Wall time
    //

#if HAVE_POSIX_CLOCK
    TimeSpec posixGetWallTime() {
      timespec result;
      if(clock_gettime(CLOCK_REALTIME, &result) < 0)
        DUNE_THROW(ClockError, "clock_gettime(CLOCK_REALTIME, ...) failed: "
                   "errno = " << errno);
      TimeSpec tmp = { result.tv_sec, result.tv_nsec };
      return tmp;
    }

    TimeSpec posixGetWallTimeResolution() {
      timespec result;
      if(clock_getres(CLOCK_REALTIME, &result) < 0)
        DUNE_THROW(ClockError, "clock_getres(CLOCK_REALTIME, ...) failed: "
                   "errno = " << errno);
      TimeSpec tmp = { result.tv_sec, result.tv_nsec };
      return tmp;
    }

    bool checkPOSIXGetWallTime() {
# if _POSIX_TIMERS == 0
      return sysconf(_SC_TIMERS) > 0;
# else // _POSIX_TIMERS > 0
      return true;
# endif // _POSIX_TIMERS > 0
    }
#endif // HAVE_POSIX_CLOCK

    TimeSpec gettimeofdayWallTime() {
      timeval result;
      if(gettimeofday(&result, NULL) < 0)
        DUNE_THROW(ClockError, "gettimeofday() failed: errno = " << errno);
      TimeSpec tmp = { result.tv_sec, 1000*result.tv_usec };
      return tmp;
    }

    const TimeSpec &gettimeofdayWallTimeResolution() {
      static const TimeSpec res = { 0, 1000 };
      return res;
    }

    struct WallTimeClock {
      TimeSpec (*clock)();
      TimeSpec resolution;
      std::string clockName;

      static const WallTimeClock &instance() {
        static const WallTimeClock clock;
        return clock;
      }

    private:
      WallTimeClock() {
#if HAVE_POSIX_CLOCK
        if(checkPOSIXGetWallTime()) {
          clock = posixGetWallTime;
          resolution = posixGetWallTimeResolution();
          clockName = "clock_gettime(CLOCK_REALTIME, ...)";
          return;
        }
#endif // HAVE_POSIX_CLOCK
        {
          clock = gettimeofdayWallTime;
          resolution = gettimeofdayWallTimeResolution();
          clockName = "gettimeofday(...)";
          return;
        }
      }
    };
    TimeSpec getWallTime() { return WallTimeClock::instance().clock(); }
    TimeSpec getWallTimeResolution()
    { return WallTimeClock::instance().resolution; }
    const std::string &getWallTimeImp()
    { return WallTimeClock::instance().clockName; }

    //////////////////////////////////////////////////////////////////////
    //
    //  Process Time
    //

#if HAVE_POSIX_CLOCK && _POSIX_CPUTIME >= 0
    TimeSpec posixGetProcessTime() {
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

    TimeSpec posixGetProcessTimeResolution() {
      timespec result;
      if(clock_getres(CLOCK_PROCESS_CPUTIME_ID, &result) < 0)
        DUNE_THROW(ClockError, "clock_getres(CLOCK_PROCESS_CPUTIME_ID, ...) "
                   "failed: errno = " << errno);
      TimeSpec tmp = { result.tv_sec, result.tv_nsec };
      return tmp;
    }

    bool checkPOSIXGetProcessTime() {
# if _POSIX_CPUTIME == 0
      return sysconf(_SC_CPUTIME) > 0;
# else // _POSIX_CPUTIME > 0
      return true;
# endif // _POSIX_CPUTIME > 0
    }
#endif // HAVE_POSIX_CLOCK && _POSIX_CPUTIME >= 0

    TimeSpec getrusageProcessTime() {
      rusage ru;
      if(getrusage(RUSAGE_SELF, &ru) < 0)
        DUNE_THROW(ClockError, "getrusage(RUSAGE_SELF, ...) failed: errno = "
                   << errno);
      TimeSpec time = { ru.ru_utime.tv_sec, 1000*ru.ru_utime.tv_usec };
      TimeSpec tmp = { ru.ru_stime.tv_sec, 1000*ru.ru_stime.tv_usec };
      time += tmp;
      return time;
    }

    const TimeSpec &getrusageProcessTimeResolution() {
      static const TimeSpec res = { 0, 1000 };
      return res;
    }

    struct ProcessTimeClock {
      TimeSpec (*clock)();
      TimeSpec resolution;
      std::string clockName;

      static const ProcessTimeClock &instance() {
        static const ProcessTimeClock clock;
        return clock;
      }

    private:
      ProcessTimeClock() {
#if HAVE_POSIX_CLOCK && _POSIX_CPUTIME >= 0
        if(checkPOSIXGetProcessTime())
        {
          clock = posixGetProcessTime;
          resolution = posixGetProcessTimeResolution();
          clockName = "clock_gettime(CLOCK_PROCESS_CPUTIME_ID, ...)";
          return;
        }
#endif // HAVE_POSIX_CLOCK && _POSIX_CPUTIME
        {
          clock = getrusageProcessTime;
          resolution = getrusageProcessTimeResolution();
          clockName = "getrusage(RUSAGE_SELF, ...)";
        }
      }
    };
    TimeSpec getProcessTime() { return ProcessTimeClock::instance().clock(); }
    TimeSpec getProcessTimeResolution()
    { return ProcessTimeClock::instance().resolution; }
    const std::string &getProcessTimeImp()
    { return ProcessTimeClock::instance().clockName; }

  } // namespace PDELab
} // namespace Dune
