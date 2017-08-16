// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_L2_HH
#define DUNE_PDELAB_LOCALOPERATOR_L2_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

#include<dune/geometry/type.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/localfunctions/common/interfaceswitch.hh>

#include<dune/pdelab/common/quadraturerules.hh>

#include"defaultimp.hh"
#include"pattern.hh"
#include"flags.hh"
#include"idefault.hh"

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
    class L2 :
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<double>
    {
    public:
      // Pattern assembly flags
      enum { doPatternVolume = true };

      // Residual assembly flags
      enum { doAlphaVolume = true };

      L2 (int intorder_=2,double scaling=1.0)
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
      void jacobian_apply_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
      {
        alpha_volume(eg,lfsu,x,lfsv,y);
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
      int intorder;
      const double _scaling;
    };

    /** A local operator for the mass operator of a vector valued lfs (L_2 integral)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     */
    class PowerL2 :
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<double>
    {
    public:
      // Pattern assembly flags
      enum { doPatternVolume = true };

      // Residual assembly flags
      enum { doAlphaVolume = true };

      PowerL2 (int intorder_=2)
        : scalar_operator(intorder_)
      {}

      // Volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        for(unsigned int i=0; i<TypeTree::degree(lfsu); ++i)
          scalar_operator.alpha_volume(eg,lfsu.child(i),x,lfsv.child(i),r);
      }

      // apply jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
      {
        alpha_volume(eg,lfsu,x,lfsv,y);
      }

      // Jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        for(unsigned int i=0; i<TypeTree::degree(lfsu); ++i)
          scalar_operator.jacobian_volume(eg,lfsu.child(i),x,lfsv.child(i),mat);
      }

    private:
      L2 scalar_operator;
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_L2_HH
