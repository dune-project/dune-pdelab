// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALOPERATOR_NAVIERSTOKESMASS_HH
#define DUNE_PDELAB_LOCALOPERATOR_NAVIERSTOKESMASS_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>

namespace Dune {
  namespace PDELab {

    /** \brief A local operator for the mass term corresponding to the
        instationary part in the Navier-Stokes equations.

        \f{align*}{
        \int_\Omega \rho u\cdot v dx
        \f}
    */
    template<typename PRM>
    class NavierStokesMass :
      public FullVolumePattern ,
      public LocalOperatorDefaultFlags ,
      public InstationaryLocalOperatorDefaultMethods<double>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      NavierStokesMass (const PRM & p_, int intorder_=4)
        : p(p_), intorder(intorder_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        typedef typename LFSV::template Child<0>::Type LFSV_PFS_V;
        const LFSV_PFS_V& lfsv_pfs_v = lfsv.template child<0>();

        for(unsigned int i=0; i<LFSV_PFS_V::CHILDREN; ++i)
          {
            scalar_alpha_volume(eg,lfsv_pfs_v.child(i),x,lfsv_pfs_v.child(i),r);
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        typedef typename LFSV::template Child<0>::Type LFSV_PFS_V;
        const LFSV_PFS_V& lfsv_pfs_v = lfsv.template child<0>();

        for(unsigned int i=0; i<LFSV_PFS_V::CHILDREN; ++i)
          {
            scalar_jacobian_volume(eg,lfsv_pfs_v.child(i),x,lfsv_pfs_v.child(i),mat);
          }
      }

    private:
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void scalar_alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                                R& r) const
      {

        // Switches between local and global interface
        typedef FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType
          > FESwitch;
        typedef BasisInterfaceSwitch<
          typename FESwitch::Basis
          > BasisSwitch;

        // domain and range field type
        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::RangeField RF;
        typedef typename BasisSwitch::Range RangeType;

        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            FESwitch::basis(lfsu.finiteElement()).evaluateFunction(ip.position(),phi);

            RF rho = p.rho(eg,ip.position());
            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // u*phi_i
            RF factor = ip.weight() * rho * eg.geometry().integrationElement(ip.position());

            for (size_type i=0; i<lfsu.size(); i++)
              r.accumulate(lfsv,i, u*phi[i]*factor);
          }
      }

      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void scalar_jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                                   M& mat) const
      {

        // Switches between local and global interface
        typedef FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType
          > FESwitch;
        typedef BasisInterfaceSwitch<
          typename FESwitch::Basis
          > BasisSwitch;

        // domain and range field type
        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::RangeField RF;
        typedef typename BasisSwitch::Range RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            FESwitch::basis(lfsu.finiteElement()).evaluateFunction(ip.position(),phi);

            // integrate phi_j*phi_i
            RF rho = p.rho(eg,ip.position());
            RF factor = ip.weight() * rho * eg.geometry().integrationElement(ip.position());
            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat.accumulate(lfsv,i,lfsu,j, phi[j]*phi[i]*factor);
          }
      }

      const PRM& p;
      int intorder;
    }; // end class NavierStokesMass

  } // end namespace PDELab
} // end namespace Dune
#endif
