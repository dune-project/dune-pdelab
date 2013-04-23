// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=8 sw=2 et sts=2:

#ifndef DUNE_PDELAB_COMMON_UNORDERED_SET_HH
#define DUNE_PDELAB_COMMON_UNORDERED_SET_HH

/** \file
    \brief Provide common name for std::unordered_set and std::unordered_multiset classes in Dune::PDELab namespace.
*/

#include <dune/common/static_assert.hh>

// Try to find an unordered_set implementation
#ifdef HAVE_UNORDERED_SET

#include <unordered_set>

#elif HAVE_TR1_UNORDERED_SET

#include <tr1/unordered_set>

#endif

namespace Dune {
  namespace PDELab {

    // import implementation into Dune::PDELab namespace if there is one.
#ifdef HAVE_UNORDERED_SET

    using std::unordered_set;

#elif HAVE_TR1_UNORDERED_SET

    using std::tr1::unordered_set;

#else

    // Dummy implementations to feed the user an explanation of what went wrong here.

    template<typename Key,
             typename Hash = int,
             typename Pred = int,
             typename Allocator = int
             >
    class unordered_set
    {
      dune_static_assert(Dune::AlwaysFalse<Key>::value,"Unable to find implementation for unordered_set.");
    };

    template<typename Key,
             typename Hash = int,
             typename Pred = int,
             typename Allocator = int
             >
    class unordered_multiset
    {
      dune_static_assert(Dune::AlwaysFalse<Key>::value,"Unable to find implementation for unordered_multiset.");
    };

#endif

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_COMMON_UNORDERED_SET_HH
