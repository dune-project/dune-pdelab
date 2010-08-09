// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_MULTISTEP_PARAMETER_HH
#define DUNE_PDELAB_MULTISTEP_PARAMETER_HH

#include <string>

namespace Dune {
  namespace PDELab {

    //! \addtogroup MultiStepMethods Multi-Step Methods
    //! \ingroup PDELab
    //! \{

    //////////////////////////////////////////////////////////////////////
    //
    //  Base class
    //

    //! Base parameter class for multi step time schemes
    /**
     * \tparam value_type_ C++ type of the floating point parameters
     * \tparam order_      Order of the ODE's this scheme is apropriate for.
     */
    template<typename value_type_, unsigned order_>
    class MultiStepParameterInterface
    {
    public:
      //! export type of the parameters
      typedef value_type_ value_type;

      //! Order of the problems this method is apropriate for
      static const unsigned order = order_;

      //! Return number of steps of the method
      virtual unsigned steps () const = 0;

      //! Return alpha coefficients
      /**
       * Return \f$\alpha_{\text{\tt step}, \text{\tt deriv}}\f$.
       *
       * \note step ∈ [0,...,steps()] and deriv ∈ [0,...,order]
       *
       * \note If a coefficient is numerically zero (\f$\alpha_{\text{\tt
       *       step}, \text{\tt deriv}}=0\f$), the MultiStepGridOperatorSpace
       *       may skip certain loops.  To take advantage of this, a
       *       particular Parameters implementation should take care to force
       *       a parameter to exactly zero before returning it, if it is
       *       practically zero and if it is calculated in a way that may
       *       result in slightly off-zero values.  This way, the meaning of
       *       "practically zero" is up to the Parameters implementation.
       */
      virtual value_type alpha(int step, int deriv) const = 0;

      //! Return name of the scheme
      virtual std::string name () const = 0;

      //! every abstract base class has a virtual destructor
      virtual ~MultiStepParameterInterface () {}
    };

    //! \} group MultiStepMethods
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTISTEP_PARAMETER_HH
