// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALOPERATOR_NAVIERSTOKESMASS_HH
#define DUNE_PDELAB_LOCALOPERATOR_NAVIERSTOKESMASS_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/utility.hh>

#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>

#include <dune/pdelab/common/quadraturerules.hh>

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
      public InstationaryLocalOperatorDefaultMethods<typename PRM::Traits::RangeField>
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
        using namespace Indices;
        const auto& lfsv_pfs_v = child(lfsv,_0);
        for(unsigned int i=0; i<TypeTree::degree(lfsv_pfs_v); ++i)
          {
            scalar_alpha_volume(eg,lfsv_pfs_v.child(i),x,lfsv_pfs_v.child(i),r);
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        using namespace Indices;
        const auto& lfsv_pfs_v = child(lfsv,_0);
        for(unsigned int i=0; i<TypeTree::degree(lfsv_pfs_v); ++i)
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
        using FESwitch = FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType>;
        using BasisSwitch = BasisInterfaceSwitch<
          typename FESwitch::Basis>;

        // Define types
        using RF = typename BasisSwitch::RangeField;
        using RangeType = typename BasisSwitch::Range;
        using size_type = typename LFSU::Traits::SizeType;

        // Dimensions
        const int dim = EG::Geometry::mydimension;

        // Get geometry
        auto geo = eg.geometry();

        // Determine quadrature order
        const int v_order = FESwitch::basis(lfsu.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        // Initialize vectors outside for loop
        std::vector<RangeType> phi(lfsu.size());

        // Loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            // evaluate basis functions
            FESwitch::basis(lfsu.finiteElement()).evaluateFunction(ip.position(),phi);

            auto rho = p.rho(eg,ip.position());
            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // u*phi_i
            auto factor = ip.weight() * rho * geo.integrationElement(ip.position());

            for (size_type i=0; i<lfsu.size(); i++)
              r.accumulate(lfsv,i, u*phi[i]*factor);
          }
      }

      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void scalar_jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                                   M& mat) const
      {

        // Switches between local and global interface
        using FESwitch = FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType>;
        using BasisSwitch = BasisInterfaceSwitch<
          typename FESwitch::Basis>;

        // Define types
        using RangeType = typename BasisSwitch::Range;
        using size_type = typename LFSU::Traits::SizeType;

        // Dimensions
        const int dim = EG::Geometry::mydimension;

        // Get geometry
        auto geo = eg.geometry();

        // Determine quadrature order
        const int v_order = FESwitch::basis(lfsu.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        // Initialize vectors outside for loop
        std::vector<RangeType> phi(lfsu.size());

        // Loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            // evaluate basis functions
            FESwitch::basis(lfsu.finiteElement()).evaluateFunction(ip.position(),phi);

            // integrate phi_j*phi_i
            auto rho = p.rho(eg,ip.position());
            auto factor = ip.weight() * rho * geo.integrationElement(ip.position());
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
      public InstationaryLocalOperatorDefaultMethods<typename PRM::Traits::RangeField>
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
        // subspaces
        using namespace Indices;
        using LFSV_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_v = child(lfsv,_0);
        const auto& lfsu_v = child(lfsu,_0);

        // Define types
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using Range_V = typename BasisSwitch_V::Range;
        using size_type = typename LFSV::Traits::SizeType;

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // Get geometry
        auto geo = eg.geometry();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        // Initialize vectors outside for loop
        std::vector<Range_V> phi_v(lfsv_v.size());
        Range_V u(0.0);

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            auto local = ip.position();
            auto rho = p.rho(eg,local);

            // compute basis functions
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

            // compute u
            u = 0.0;
            for(size_type i=0; i<lfsu_v.size(); i++)
              u.axpy(x(lfsu_v,i),phi_v[i]);

            auto factor = ip.weight() * rho * geo.integrationElement(ip.position());

            for(size_type i=0; i<lfsv_v.size(); i++)
              r.accumulate(lfsv_v,i, (u*phi_v[i]) * factor);

          } // end loop quadrature points
      } // end alpha_volume

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // subspaces
        using namespace Indices;
        using LFSV_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_v = child(lfsv,_0);
        const auto& lfsu_v = child(lfsu,_0);

        // Define types
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using Range_V = typename BasisSwitch_V::Range;
        using size_type = typename LFSV::Traits::SizeType;

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // Get geometry
        auto geo = eg.geometry();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        // Initialize vectors outside for loop
        std::vector<Range_V> phi_v(lfsv_v.size());

        // Loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            auto local = ip.position();
            auto rho = p.rho(eg,local);

            // compute basis functions
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

            auto factor = ip.weight() * rho * geo.integrationElement(ip.position());

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
#endif // DUNE_PDELAB_LOCALOPERATOR_NAVIERSTOKESMASS_HH
