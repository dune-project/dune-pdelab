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

      NavierStokesMass (const PRM & p_, int superintegration_order_ = 0)
        : p(p_), superintegration_order(superintegration_order_)
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
        const int v_order = FESwitch::basis(lfsu.finiteElement()).order();
        const int det_jac_order = gt.isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

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
        const int v_order = FESwitch::basis(lfsu.finiteElement()).order();
        const int det_jac_order = gt.isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

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
      const int superintegration_order;
    }; // end class NavierStokesMass

    /** \brief A local operator for the mass term corresponding to the
        instationary part in the Navier-Stokes equations
        using a vector-valued Finite Element map for the velocity grid function space.

        \f{align*}{
        \int_\Omega \rho u\cdot v dx
        \f}
    */
    template<typename PRM>
    class NavierStokesVelVecMass :
      public FullVolumePattern ,
      public LocalOperatorDefaultFlags ,
      public InstationaryLocalOperatorDefaultMethods<double>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      NavierStokesVelVecMass (const PRM & p_, int superintegration_order_ = 0)
        : p(p_), superintegration_order(superintegration_order_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // dimensions
        const int dim = EG::Geometry::dimension;

        // subspaces
        typedef typename LFSV::template Child<0>::Type LFSV_V;
        const LFSV_V& lfsv_v = lfsv.template child<0>();
        const LFSV_V& lfsu_v = lfsu.template child<0>();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef typename BasisSwitch_V::Range Range_V;
        typedef typename LFSV::Traits::SizeType size_type;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = gt.isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            const Dune::FieldVector<DF,dim> local = ip.position();
            const RF rho = p.rho(eg,local);

            // compute basis functions
            std::vector<Range_V> phi_v(lfsv_v.size());
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

            // compute u
            Range_V u(0.0);
            for(size_type i=0; i<lfsu_v.size(); i++)
              u.axpy(x(lfsu_v,i),phi_v[i]);

            const RF factor = ip.weight() * rho * eg.geometry().integrationElement(ip.position());

            for(size_type i=0; i<lfsv_v.size(); i++)
              r.accumulate(lfsv_v,i, (u*phi_v[i]) * factor);

          } // end loop quadrature points
      } // end alpha_volume

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // dimensions
        const int dim = EG::Geometry::dimension;

        // subspaces
        typedef typename LFSV::template Child<0>::Type LFSV_V;
        const LFSV_V& lfsv_v = lfsv.template child<0>();
        const LFSV_V& lfsu_v = lfsu.template child<0>();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef typename BasisSwitch_V::Range Range_V;
        typedef typename LFSV::Traits::SizeType size_type;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = gt.isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            const Dune::FieldVector<DF,dim> local = ip.position();
            const RF rho = p.rho(eg,local);

            // compute basis functions
            std::vector<Range_V> phi_v(lfsv_v.size());
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

            const RF factor = ip.weight() * rho * eg.geometry().integrationElement(ip.position());

            for(size_type i=0; i<lfsv_v.size(); i++)
              for(size_type j=0; j<lfsu_v.size(); j++)
                mat.accumulate(lfsv_v,i,lfsu_v,j, (phi_v[j]*phi_v[i]) * factor);
          } // end loop quadrature points
      } // end jacobian_volume

    private :
      const PRM& p;
      const int superintegration_order;
    }; // end class NavierStokesVelVecMass

  } // end namespace PDELab
} // end namespace Dune
#endif
