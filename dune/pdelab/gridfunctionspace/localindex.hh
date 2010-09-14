// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALINDEX_HH
#define DUNE_PDELAB_LOCALINDEX_HH

#include "localfunctionspacetags.hh"

/** \file
    \author Christian Engwer
    A local Index class, which can be tagged by a tag from localfunctionspacetags.hh
 */

namespace Dune {
  namespace PDELab {
    
    /**
       \addtogroup PDELAB_StrictTrialAndTest Strict Trial and Test space handling
       \ingroup GridFunctionSpace
       \{
    */

    /**
       \brief hide an integer in a wrapper class. 

       This makes it possible to use strict type checking, even for integers.
     */
    template<typename I, typename TAG>
    struct LocalIndex
    {
      I i;
      LocalIndex(const I & _i) : i(_i) {}
    };

    /**
       \brief map a given tag from localfunctionspacetags.hh to the matching LocalIndex class
     */
    template<typename I, TAGENAME TAG>
    struct LocalIndexTraits
    {
      typedef Dune::PDELab::LocalIndex<I, TAG> LocalIndex;
    };

    /**
       For AnySpaceTag the LocalIndex is just the raw integer
     */
    template<typename I>
    struct LocalIndexTraits<AnySpaceTag>
    {
      typedef I LocalIndex;
    };

    /**
       \}
     */

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_LOCALINDEX_HH
