// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_L2_HH
#define DUNE_PDELAB_LOCALOPERATOR_L2_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/common/quadraturerules.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/blockdiagonal.hh>

namespace Dune {
  namespace PDELab {

    namespace impl {

      // Scalar L2 operator. Only for internal use! Use the L2 class instead,
      // as that will also work for non-scalar spaces.
      class ScalarL2 :
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<double>
      {
      public:
        // Pattern assembly flags
        enum { doPatternVolume = true };

        // Residual assembly flags
        enum { doAlphaVolume = true };

        ScalarL2 (int intorderadd, double scaling)
          : _intorderadd(intorderadd)
          , _scaling(scaling)
        {}

        // Volume integral depending on test and ansatz functions
        template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
        void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
        {
          // Switches between local and global interface
          using FESwitch = FiniteElementInterfaceSwitch<
            typename LFSU::Traits::FiniteElementType>;
          using BasisSwitch = BasisInterfaceSwitch<
            typename FESwitch::Basis>;

          // Define types
          using RF = typename BasisSwitch::RangeField;
          using RangeType = typename BasisSwitch::Range;
          using size_type = typename LFSU::Traits::SizeType;

          // Get geometry
          auto geo = eg.geometry();

          // Initialize vectors outside for loop
          std::vector<RangeType> phi(lfsu.size());

          // determine integration order
          auto intorder = 2*FESwitch::basis(lfsu.finiteElement()).order() + _intorderadd;

          // Loop over quadrature points
          for (const auto& qp : quadratureRule(geo,intorder))
            {
              // Evaluate basis functions
              FESwitch::basis(lfsu.finiteElement()).evaluateFunction(qp.position(),phi);

              // Evaluate u
              RF u=0.0;
              for (size_type i=0; i<lfsu.size(); i++)
                u += RF(x(lfsu,i)*phi[i]);

              // u*phi_i
              auto factor = _scaling * qp.weight() * geo.integrationElement(qp.position());
              for (size_type i=0; i<lfsu.size(); i++)
                r.accumulate(lfsv,i, u*phi[i]*factor);
            }
        }

        // apply jacobian of volume term
        template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
        void jacobian_apply_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
        {
          alpha_volume(eg,lfsu,x,lfsv,y);
        }

        // Jacobian of volume term
        template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
        void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat) const
        {
          // Switches between local and global interface
          using FESwitch = FiniteElementInterfaceSwitch<
            typename LFSU::Traits::FiniteElementType>;
          using BasisSwitch = BasisInterfaceSwitch<
            typename FESwitch::Basis>;

          // Define types
          using RangeType = typename BasisSwitch::Range;
          using size_type = typename LFSU::Traits::SizeType;

          // Get geometry
          auto geo = eg.geometry();

          // Inititialize vectors outside for loop
          std::vector<RangeType> phi(lfsu.size());

          // determine integration order
          auto intorder = 2*FESwitch::basis(lfsu.finiteElement()).order() + _intorderadd;

          // Loop over quadrature points
          for (const auto& qp : quadratureRule(geo,intorder))
            {
              // Evaluate basis functions
              FESwitch::basis(lfsu.finiteElement()).evaluateFunction(qp.position(),phi);

              // Integrate phi_j*phi_i
              auto factor = _scaling * qp.weight() * geo.integrationElement(qp.position());
              for (size_type j=0; j<lfsu.size(); j++)
                for (size_type i=0; i<lfsu.size(); i++)
                  mat.accumulate(lfsv,i,lfsu,j, phi[j]*phi[i]*factor);
            }
        }

      private:
        int _intorderadd;
        double _scaling;
      };

    } // namespace impl

    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** A local operator for the mass operator (L_2 integral)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     *
     * This operator also works for trees of function spaces by applying
     * the L2 operator on the block diagonal.
     */
    class L2 :
      public BlockDiagonalLocalOperatorFullCoupling<impl::ScalarL2>
    {

    public:

      //! Constructs a new L2 operator.
      /**
       * This constructor creates a new L2 operator.
       *
       * \param intorderadd By default, the operator will use the sum of the degrees of the ansatz
       *                    and test space as its integration order. This parameter gets added to
       *                    that value and lets you modify the default.
       * \param scaling     The output of the operator will be scaled by this value.
       */
      L2 (int intorderadd = 0, double scaling = 1.0)
        : BlockDiagonalLocalOperatorFullCoupling<impl::ScalarL2>(intorderadd,scaling)
      {}

    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_L2_HH
