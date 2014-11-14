// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_NULLMUTEX_HH
#define DUNE_PDELAB_COMMON_NULLMUTEX_HH

#include <dune/common/deprecated.hh>

namespace Dune {
  namespace PDELab {

    //! Dummy mutex
    /**
     * This is a type that fulfills the BasicLockable, Lockable, and
     * TimedLockable concepts but will always grant the lock immediately.
     *
     * \note The purpose is mainly for benchmarking to be able to see the the
     *       effect of locking on timing.  Since there isn't actually any
     *       locking happening, the result will likely be bogus.
     */
    class NullMutex
    {
    public:
#if ! SILENCE_NULLMUTEX_WARNING
      DUNE_DEPRECATED_MSG("Warning: Dune::PDELab::NullMutex should only be "
                          "used for debugging/benchmarking purposes, any "
                          "actual results will probably be bogus.")
#endif
        NullMutex();
      inline void lock() {}
      inline void unlock() {}
      inline bool try_lock() { return true; }
      template<class Duration>
      bool try_lock_for(const Duration &) { return true; }
      template<class Time>
      bool try_lock_until(const Time &) { return true; }
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_NULLMUTEX_HH
