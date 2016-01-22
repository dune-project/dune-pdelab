// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_L2_HH
#define DUNE_PDELAB_L2_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

#include<dune/geometry/type.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>

#include<dune/localfunctions/common/interfaceswitch.hh>

#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/common/referenceelements.hh>

#include"defaultimp.hh"
#include"pattern.hh"
#include"flags.hh"
#include"idefault.hh"

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** a local operator for the mass operator (L_2 integral)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     */
    class L2 : public NumericalJacobianApplyVolume<L2>,
               public FullVolumePattern,
               public LocalOperatorDefaultFlags,
               public InstationaryLocalOperatorDefaultMethods<double>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      L2 (int intorder_=2,double scaling=1.0)
        : intorder(intorder_)
        , _scaling(scaling)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // Switches between local and global interface
        typedef FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType
          > FESwitch;
        typedef BasisInterfaceSwitch<
          typename FESwitch::Basis
          > BasisSwitch;

        // domain and range field type
        typedef typename BasisSwitch::RangeField RF;
        typedef typename BasisSwitch::Range RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // select quadrature rule
        auto geo = eg.geometry();

        // Initialize outside for loop
        std::vector<RangeType> phi(lfsu.size());

        // loop over quadrature points
        for (const auto& qp : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            FESwitch::basis(lfsu.finiteElement()).evaluateFunction(qp.position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += RF(x(lfsu,i)*phi[i]);

            // u*phi_i
            auto factor = _scaling * qp.weight() * geo.integrationElement(qp.position());
            for (size_type i=0; i<lfsu.size(); i++)
              r.accumulate(lfsv,i, u*phi[i]*factor);
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M & mat) const
      {
        // Switches between local and global interface
        typedef FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType
          > FESwitch;
        typedef BasisInterfaceSwitch<
          typename FESwitch::Basis
          > BasisSwitch;

        // domain and range field type
        typedef typename BasisSwitch::Range RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // select quadrature rule
        auto geo = eg.geometry();

        // Inititialize outside for loop
        std::vector<RangeType> phi(lfsu.size());

        // loop over quadrature points
        for (const auto& qp : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            FESwitch::basis(lfsu.finiteElement()).evaluateFunction(qp.position(),phi);

            // integrate phi_j*phi_i
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

    /** a local operator for the mass operator of a vector valued lfs (L_2 integral)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     */
    class PowerL2 : public NumericalJacobianApplyVolume<PowerL2>,
                    public FullVolumePattern,
                    public LocalOperatorDefaultFlags,
                    public InstationaryLocalOperatorDefaultMethods<double>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      PowerL2 (int intorder_=2)
        : scalar_operator(intorder_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        for(unsigned int i=0; i<LFSU::CHILDREN; ++i)
          scalar_operator.alpha_volume(eg,lfsu.child(i),x,lfsv.child(i),r);
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        for(unsigned int i=0; i<LFSU::CHILDREN; ++i)
          scalar_operator.jacobian_volume(eg,lfsu.child(i),x,lfsv.child(i),mat);
      }

    private:
      L2 scalar_operator;
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
