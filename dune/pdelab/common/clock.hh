// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_CLOCK_HH
#define DUNE_PDELAB_COMMON_CLOCK_HH

#include <ostream>

#include <sys/types.h>

#include <dune/common/exceptions.hh>

namespace Dune {
  namespace PDELab {

    //! struct to store temporal values
    struct TimeSpec {
      //! seconds part
      time_t tv_sec;
      //! nanoseconds part
      /**
       * \note 0 <= tv_nsec < 1e9 should always hold.
       */
      long tv_nsec;

      inline TimeSpec &operator+=(const TimeSpec &o) {
        tv_sec += o.tv_sec;
        tv_nsec += o.tv_nsec;
        if(tv_nsec >= 1000000000L) {
          ++tv_sec;
          tv_nsec -= 1000000000L;
        }
        return *this;
      }
      inline TimeSpec operator+(const TimeSpec &o) const {
        TimeSpec tmp(*this);
        tmp += o;
        return tmp;
      }
      inline TimeSpec &operator-=(const TimeSpec &o) {
        tv_sec -= o.tv_sec;
        tv_nsec -= o.tv_nsec;
        if(tv_nsec < 0L) {
          --tv_sec;
          tv_nsec += 1000000000L;
        }
        return *this;
      }
      inline TimeSpec operator-(const TimeSpec &o) const {
        TimeSpec tmp(*this);
        tmp -= o;
        return tmp;
      }

    };

    //! insert a timespec into an output stream
    std::ostream &operator<<(std::ostream &s, const TimeSpec &t);

    //! exception thrown by clock functions
    struct ClockError : Exception {};

    //! get the wall time in seconds since the epoch
    TimeSpec getWallTime();
    //! get resolution of the wall time in seconds
    TimeSpec getWallTimeResolution();
    //! \brief return a string describing which implementation is used to get
    //!        the wall time
    const std::string &getWallTimeImp();
    //! get the process time in seconds used by the current process
    TimeSpec getProcessTime();
    //! get resolution of the process time in seconds
    TimeSpec getProcessTimeResolution();
    //! \brief return a string describing which implementation is used to get
    //!        the process time
    const std::string &getProcessTimeImp();

  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_CLOCK_HH

