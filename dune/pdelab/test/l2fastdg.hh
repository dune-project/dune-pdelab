// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_L2FASTDG_HH
#define DUNE_PDELAB_LOCALOPERATOR_L2FASTDG_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

#include<dune/geometry/type.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/localfunctions/common/interfaceswitch.hh>

#include<dune/pdelab/common/quadraturerules.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** A local operator for the mass operator (L_2 integral)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     */
    class L2FastDG : public NumericalJacobianApplyVolume<L2FastDG>,
                     public FullVolumePattern,
                     public LocalOperatorDefaultFlags,
                     public InstationaryLocalOperatorDefaultMethods<double>
    {
    public:
      // Pattern assembly flags
      enum { doPatternVolume = true };

      // Residual assembly flags
      enum { doAlphaVolume = true };

      L2FastDG (int intorder_=2,double scaling=1.0)
        : intorder(intorder_)
        , _scaling(scaling)
      {}

      // Volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
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

        // Weight of time stepping
        auto t_weight = r.weight();

        // Loop over quadrature points
        for (const auto& qp : quadratureRule(geo,intorder))
          {
            // Evaluate basis functions
            FESwitch::basis(lfsu.finiteElement()).evaluateFunction(qp.position(),phi);

            // Evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += RF(x[i]*phi[i]);

            // u*phi_i
            auto factor = _scaling * qp.weight() * geo.integrationElement(qp.position()) * t_weight;
            for (size_type i=0; i<lfsu.size(); i++)
              r.data()[i] += (u*phi[i]*factor);
          }
      }

      // Jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M & mat) const
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

        // Weight of time stepping
        auto t_weight = mat.weight();

        // Loop over quadrature points
        for (const auto& qp : quadratureRule(geo,intorder))
          {
            // Evaluate basis functions
            FESwitch::basis(lfsu.finiteElement()).evaluateFunction(qp.position(),phi);

            // Integrate phi_j*phi_i
            auto factor = _scaling * qp.weight() * geo.integrationElement(qp.position()) * t_weight;
            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat.data()[lfsu.size()*i+j] += (phi[j]*phi[i]*factor);
          }
      }

    private:
      int intorder;
      const double _scaling;
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_L2FASTDG_HH
