// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_UTILITY_HH
#define DUNE_PDELAB_COMMON_UTILITY_HH

#include <dune/common/shared_ptr.hh>

namespace Dune {
  namespace PDELab {


    /** \addtogroup common Common Utilities
     *  \ingroup PDELab
     *  \{
     */

    //! Ensures that t is wrapped in a shared_ptr<T>
    /**
     * You have to consider three situations:
     *
     * a) t is of type T&
     *  t is a stack object and must not be deleted.
     *  You create a shared_ptr<T>(&t) with a null_deleter.
     *  b) t is of type T*
     *  t is a raw pointer and the user os assumed to own this pointer.
     *  You create a shared_ptr<T>(t) with a null_deleter.
     *  c) t is of type shared_ptr<T>
     *  t is already a shared_ptr<T>.
     *  You don't have to do anything.
     */
    template<typename T>
    shared_ptr<T> ensure_shared_ptr(T & t)
    {
      return shared_ptr<T>(&t, null_deleter<T>());
    }

#ifndef DOXYGEN

    template<typename T>
    shared_ptr<T> ensure_shared_ptr(T * t)
    {
      return shared_ptr<T>(t, null_deleter<T>());
    }

    template<typename T>
    shared_ptr<T> & ensure_shared_ptr(shared_ptr<T> & t)
    {
      return t;
    }

#endif // DOXYGEN

    //! \}

  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_UTILITY_HH
