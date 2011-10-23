// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALOPERATOR_ZERO_HH
#define DUNE_PDELAB_LOCALOPERATOR_ZERO_HH

#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    //! Do-nothing local operator
    /**
     * This local operator does nothing.  It can be used as a
     * placeholder where an (instationary) local operator is semantically
     * required but not actually used.
     *
     * \tparam Time Type used for temporal values.
     */
    template<class Time>
    class ZeroLocalOperator :
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<Time>
    { };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_ZERO_HH
