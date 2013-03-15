// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=8 sw=2 et sts=2:

#ifndef DUNE_PDELAB_COMMON_UNORDERED_MAP_HH
#define DUNE_PDELAB_COMMON_UNORDERED_MAP_HH

/** \file
    \brief Provide common name for std::unordered_map and std::unordered_multimap classes in Dune::PDELab namespace.
*/

#include <dune/common/static_assert.hh>

// Try to find an unordered_map implementation
#ifdef HAVE_UNORDERED_MAP

#include <unordered_map>

#elif HAVE_TR1_UNORDERED_MAP

#include <tr1/unordered_map>

#endif

namespace Dune {
  namespace PDELab {

    // import implementation into Dune::PDELab namespace if there is one.
#ifdef HAVE_UNORDERED_MAP

    using std::unordered_map;

#elif HAVE_TR1_UNORDERED_MAP

    using std::tr1::unordered_map;

#else

    // Dummy implementations to feed the user an explanation of what went wrong here.

    template<typename Key,
             typename T,
             typename Hash = int,
             typename Pred = int,
             typename Allocator = int
             >
    class unordered_map
    {
      dune_static_assert(Dune::AlwaysFalse<Key>::value,"Unable to find implementation for unordered_map.");
    };

    template<typename Key,
             typename T,
             typename Hash = int,
             typename Pred = int,
             typename Allocator = int
             >
    class unordered_multimap
    {
      dune_static_assert(Dune::AlwaysFalse<Key>::value,"Unable to find implementation for unordered_multimap.");
    };

#endif

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_COMMON_UNORDERED_MAP_HH
