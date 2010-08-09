// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_MULTISTEP_PARAMETER_HH
#define DUNE_PDELAB_MULTISTEP_PARAMETER_HH

#include <sstream>
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

    //////////////////////////////////////////////////////////////////////
    //
    //  Central Differences
    //

    //! Parameter class for the central differences scheme
    /**
     * \tparam value_type C++ type of the floating point parameters
     */
    template<typename value_type>
    class CentralDifferencesParameters
      : public MultiStepParameterInterface<value_type, 2>
    {
      static const value_type a[3][3];

    public:
      //! Return number of steps of the method
      /**
       * \returns 2 for central differences
       */
      virtual unsigned steps () const { return 2; }

      //! Return alpha coefficients
      /**
       * Return \f$\alpha_{\text{\tt step}, \text{\tt deriv}}\f$:
       * \f{align*}{
       *   \alpha_{00}&=0 & \alpha_{01}&=\frac12 & \alpha_{02}&=1 \\
       *   \alpha_{10}&=1 & \alpha_{11}&=0       & \alpha_{12}&=-2 \\
       *   \alpha_{20}&=0 & \alpha_{21}&=\frac12 & \alpha_{22}&=1
       * \f}
       *
       * \note step ∈ [0,...,steps()] and deriv ∈ [0,...,order]
       */
      virtual value_type alpha(int step, int deriv) const {
        return a[step][deriv];
      }

      //! Return name of the scheme
      virtual std::string name () const {
        return "Central Differences";
      }
    };

    template<typename value_type>
    const value_type CentralDifferencesParameters<value_type>::a[3][3] = {
      {0, 0.5,  1},
      {1, 0,   -2},
      {0, 0.5,  1}
    };

    //////////////////////////////////////////////////////////////////////
    //
    //  Newmark-β scheme
    //

    //! Parameter class for the Newmark-β scheme
    /**
     * \tparam value_type C++ type of the floating point parameters
     */
    template<typename value_type>
    class NewmarkBetaParameters
      : public MultiStepParameterInterface<value_type, 2>
    {
      value_type a[3][3];
      std::string name_;

    public:
      NewmarkBetaParameters(value_type beta) {
        a[0][0]=beta;     a[0][1]=0.5; a[0][2]=1;
        a[1][0]=1-2*beta; a[1][1]=0;   a[1][2]=-2;
        a[2][0]=beta;     a[2][1]=0.5; a[2][2]=1;

        std::ostringstream s;
        s << "Newmark-β (β=" << beta << ")";
        name_ = s.str();
      }

      //! Return number of steps of the method
      /**
       * \returns 2 for Newmark-β
       */
      virtual unsigned steps () const { return 2; }

      //! Return alpha coefficients
      /**
       * Return \f$\alpha_{\text{\tt step}, \text{\tt deriv}}\f$:
       * \f{align*}{
       *   \alpha_{00}&=\beta    & \alpha_{01}&=\frac12 & \alpha_{02}&=1 \\
       *   \alpha_{10}&=1-2\beta & \alpha_{11}&=0       & \alpha_{12}&=-2 \\
       *   \alpha_{20}&=\beta    & \alpha_{21}&=\frac12 & \alpha_{22}&=1
       * \f}
       *
       * \note step ∈ [0,...,steps()] and deriv ∈ [0,...,order]
       */
      virtual value_type alpha(int step, int deriv) const {
        return a[step][deriv];
      }

      //! Return name of the scheme
      virtual std::string name () const {
        return name_;
      }
    };

    //! \} group MultiStepMethods
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTISTEP_PARAMETER_HH
