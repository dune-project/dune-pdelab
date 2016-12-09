// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALOPERATOR_DGNAVIERSTOKESVELVECFEM_HH
#define DUNE_PDELAB_LOCALOPERATOR_DGNAVIERSTOKESVELVECFEM_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/pdelab/localoperator/idefault.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/dgnavierstokesparameter.hh>
#include <dune/pdelab/localoperator/navierstokesmass.hh>

#ifndef VBLOCK
#define VBLOCK 0
#endif
#define PBLOCK (- VBLOCK + 1)

namespace Dune {
  namespace PDELab {

  template<class Basis, class Dummy = void>
  struct VectorBasisInterfaceSwitch {
    //! export vector type of the local coordinates
    using DomainLocal = typename Basis::Traits::DomainLocal;
    //! export field type of the values
    using RangeField = typename Basis::Traits::RangeField;
    //! export dimension of the values
    static const std::size_t dimRange = Basis::Traits::dimRange;

    //! Compute global jacobian matrix for vector valued bases
    template<typename Geometry>
    static void jacobian(const Basis& basis, const Geometry& geometry,
                         const DomainLocal& xl,
                         std::vector<FieldMatrix<RangeField, dimRange,
                                          Geometry::coorddimension> >& jac)
    {
      jac.resize(basis.size());
      basis.evaluateJacobian(xl, jac);
    }
  };

  //! Switch for uniform treatment of local and global basis classes
  template<class Basis>
  struct VectorBasisInterfaceSwitch<
    Basis, typename std::enable_if<
             Std::to_true_type<
               std::integral_constant<
                 std::size_t,
                 Basis::Traits::dimDomain
                 >
               >::value
             >::type
    >
  {
    //! export vector type of the local coordinates
    using DomainLocal = typename Basis::Traits::DomainType;
    //! export field type of the values
    using RangeField = typename Basis::Traits::RangeFieldType;
    //! export dimension of the values
    static const std::size_t dimRange = Basis::Traits::dimRange;

    //! Compute global jacobian matrix for vector valued bases
    template<typename Geometry>
    static void jacobian(const Basis& basis, const Geometry& geometry,
                         const DomainLocal& xl,
                         std::vector<FieldMatrix<RangeField, dimRange,
                                          Geometry::coorddimension> >& jac)
    {
      jac.resize(basis.size());

      std::vector<FieldMatrix<
      RangeField, dimRange, Geometry::coorddimension> > ljac(basis.size());
      basis.evaluateJacobian(xl, ljac);

      const typename Geometry::JacobianInverseTransposed& geojac =
        geometry.jacobianInverseTransposed(xl);

      for(std::size_t i = 0; i < basis.size(); ++i)
        for(std::size_t row=0; row < dimRange; ++row)
          geojac.mv(ljac[i][row], jac[i][row]);
    }
  };

    /** \brief A local operator for solving the Navier-Stokes equations
        using a DG discretization and a vector-valued finite element map
        for the velocity grid function space.

        \tparam PRM Parameter class for this local operator.

    */
    template<typename PRM>
    class DGNavierStokesVelVecFEM :
      public LocalOperatorDefaultFlags,
      public FullSkeletonPattern, public FullVolumePattern,
      public InstationaryLocalOperatorDefaultMethods<typename PRM::Traits::RangeField>
    {
      using BC = StokesBoundaryCondition;
      using RF = typename PRM::Traits::RangeField;

      using InstatBase = InstationaryLocalOperatorDefaultMethods<typename PRM::Traits::RangeField>;
      using Real = typename InstatBase::RealType;

      static const bool navier = PRM::assemble_navier;
      static const bool full_tensor = PRM::assemble_full_tensor;

    public :

      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // call the assembler for each face only once
      enum { doSkeletonTwoSided = false };

      // residual assembly flags
      enum { doAlphaVolume    = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume   = true };

      /** \brief Constructor

          \param [in] _prm                        Parameter class for this local operator
          \param [in] _superintegration_order     This number will be added to the order of
                                                  quadrature in every integration. It is
                                                  only needed, when one of the parameters (e.g
                                                  rho, mu) is not constant or the mappings from
                                                  the reference elements to the cells are
                                                  nonlinear. Boundary conditions are assumed to
                                                  have the same order as the corresponding
                                                  finite element.
      */
      DGNavierStokesVelVecFEM (PRM& _prm, int _superintegration_order=0) :
        prm(_prm), superintegration_order(_superintegration_order),
        current_dt(1.0)
      {}

      // Store current dt
      void preStep (Real , Real dt, int )
      {
        current_dt = dt;
      }

      // set time in parameter class
      void setTime(Real t)
      {
        InstatBase::setTime(t);
        prm.setTime(t);
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        const unsigned int dim = EG::Geometry::mydimension;

        // subspaces
        using namespace Indices;
        using LFSV_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_v = child(lfsv,_0);
        const auto& lfsu_v = child(lfsu,_0);

        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_p = child(lfsv,_1);
        const auto& lfsu_p = child(lfsu,_1);

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using VectorBasisSwitch_V = VectorBasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;
        using BasisSwitch_P = BasisInterfaceSwitch<typename FESwitch_P::Basis >;
        using RF = typename BasisSwitch_V::RangeField;
        using Range_V = typename BasisSwitch_V::Range;
        using Range_P = typename BasisSwitch_P::Range;
        using size_type = typename LFSV::Traits::SizeType;

        // Get geometry
        auto geo = eg.geometry();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int jac_order = geo.type().isSimplex() ? 0 : 1;
        const int qorder = 3*v_order - 1 + jac_order + det_jac_order + superintegration_order;

        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            auto local = ip.position();
            auto mu = prm.mu(eg,local);
            auto rho = prm.rho(eg,local);

            // compute u (if Navier term enabled)
            std::vector<Range_V> phi_v(lfsv_v.size());
            Range_V val_u(0.0);
            if(navier) {
              FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);
              for (size_type i=0; i<lfsu_v.size(); i++)
                val_u.axpy(x(lfsu_v,i),phi_v[i]);
            }

            // values of pressure shape functions
            std::vector<Range_P> phi_p(lfsv_p.size());
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            // compute pressure value
            Range_P val_p(0.0);
            for (size_type i=0; i<lfsu_p.size(); i++)
              val_p.axpy(x(lfsu_p,i),phi_p[i]);

            // evaluate jacobian of velocity shape functions on reference element
            std::vector<Dune::FieldMatrix<RF,dim,dim> > jac_phi_v(lfsu_v.size());
            VectorBasisSwitch_V::jacobian
              (FESwitch_V::basis(lfsv_v.finiteElement()), geo, local, jac_phi_v);

            // compute divergence of test functions
            std::vector<RF> div_phi_v(lfsv_v.size(),0.0);
            for (size_type i=0; i<lfsv_v.size(); i++)
              for (size_type d=0; d<dim; d++)
                div_phi_v[i] += jac_phi_v[i][d][d];

            // compute velocity jacobian and divergence
            Dune::FieldMatrix<RF,dim,dim> jac_u(0.0);
            RF div_u(0.0);
            for (size_type i=0; i<lfsu_v.size(); i++){
              jac_u.axpy(x(lfsu_v,i),jac_phi_v[i]);
              div_u += x(lfsu_v,i) * div_phi_v[i];
            }

            auto detj = geo.integrationElement(ip.position());
            auto weight = ip.weight() * detj;

            for (size_type i=0; i<lfsv_v.size(); i++) {
              //================================================//
              // \int (mu*grad_u*grad_v)
              //================================================//
              RF dvdu(0);
              contraction(jac_u,jac_phi_v[i],dvdu);
              r.accumulate(lfsv_v, i, dvdu * mu * weight);

              //================================================//
              // \int -p \nabla\cdot v
              //================================================//
              r.accumulate(lfsv_v, i, - div_phi_v[i] * val_p * weight);

              //================================================//
              // \int \rho ((u\cdot\nabla ) u )\cdot v
              //================================================//
              if(navier) {
                // compute (grad u) u (matrix-vector product)
                Range_V nabla_u_u(0.0);
                jac_u.mv(val_u,nabla_u_u);
                r.accumulate(lfsv_v, i, rho * (nabla_u_u*phi_v[i]) * weight);
              } // end navier

            } // end i

            for (size_type i=0; i<lfsv_p.size(); i++) {
              //================================================//
              // \int -q \nabla\cdot u
              //================================================//
              r.accumulate(lfsv_p, i, - div_u * phi_p[i] * incomp_scaling * weight);
            }

          } // end loop quadrature points
      } // end alpha_volume

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            LocalMatrix& mat) const
      {
        const unsigned int dim = EG::Geometry::mydimension;

        // subspaces
        using namespace Indices;
        using LFSV_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_v = child(lfsv,_0);
        const auto& lfsu_v = child(lfsu,_0);

        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_p = child(lfsv,_1);
        const auto& lfsu_p = child(lfsu,_1);

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using VectorBasisSwitch_V = VectorBasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;
        using BasisSwitch_P = BasisInterfaceSwitch<typename FESwitch_P::Basis >;
        using RF = typename BasisSwitch_V::RangeField;
        using Range_V = typename BasisSwitch_V::Range;
        using Range_P = typename BasisSwitch_P::Range;
        using size_type = typename LFSV::Traits::SizeType;

        // Get geometry
        auto geo = eg.geometry();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int jac_order = geo.type().isSimplex() ? 0 : 1;
        const int qorder = 3*v_order - 1 + jac_order + det_jac_order + superintegration_order;

        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            auto local = ip.position();
            auto mu = prm.mu(eg,local);
            auto rho = prm.rho(eg,local);

            // compute u (if Navier term enabled)
            std::vector<Range_V> phi_v(lfsv_v.size());
            Range_V val_u(0.0);
            if(navier) {
              FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);
              for (size_type i=0; i<lfsu_v.size(); i++)
                val_u.axpy(x(lfsu_v,i),phi_v[i]);
            }

            // values of pressure shape functions
            std::vector<Range_P> phi_p(lfsv_p.size());
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            // evaluate jacobian of velocity shape functions on reference element
            std::vector<Dune::FieldMatrix<RF,dim,dim> > jac_phi_v(lfsu_v.size());
            VectorBasisSwitch_V::jacobian
              (FESwitch_V::basis(lfsv_v.finiteElement()), geo, local, jac_phi_v);

            assert(lfsu_v.size() == lfsv_v.size());
            // compute divergence of velocity shape functions
            std::vector<RF> div_phi_v(lfsv_v.size(),0.0);
            for (size_type i=0; i<lfsv_v.size(); i++)
              for (size_type d=0; d<dim; d++)
                div_phi_v[i] += jac_phi_v[i][d][d];

            // compute velocity jacobian (if Navier term enabled)
            Dune::FieldMatrix<RF,dim,dim> jac_u(0.0);
            if(navier) {
              for (size_type i=0; i<lfsu_v.size(); i++){
                jac_u.axpy(x(lfsu_v,i),jac_phi_v[i]);
              }
            }

            auto detj = geo.integrationElement(ip.position());
            auto weight = ip.weight() * detj;

            for(size_type i=0; i<lfsv_v.size(); i++) {

              for(size_type j=0; j<lfsu_v.size(); j++) {
                //================================================//
                // \int (mu*grad_u*grad_v)
                //================================================//
                RF dvdu(0.0);
                contraction(jac_phi_v[j],jac_phi_v[i],dvdu);
                mat.accumulate(lfsv_v,i,lfsu_v,j, mu * dvdu * weight);

                //================================================//
                // \int \rho ((u\cdot\nabla ) u )\cdot v
                //================================================//
                if(navier) {
                  // compute (grad u) phi_v_j (matrix-vector product)
                  Range_V nabla_u_phi_v_j(0.0);
                  jac_u.mv(phi_v[j],nabla_u_phi_v_j);
                  // compute (grad phi_v_j) u (matrix-vector product)
                  Range_V nabla_phi_v_j_u(0.0);
                  jac_phi_v[j].mv(val_u,nabla_phi_v_j_u);
                  mat.accumulate(lfsv_v,i,lfsu_v,j, rho * ((nabla_u_phi_v_j*phi_v[i]) + (nabla_phi_v_j_u*phi_v[i])) * weight);
                } // end navier

              } // end j

              for(size_type j=0; j<lfsv_p.size(); j++) {
                //================================================//
                // - p * div v
                // - q * div u
                //================================================//
                mat.accumulate(lfsv_v,i,lfsu_p,j, -phi_p[j] * div_phi_v[i] * weight);
                mat.accumulate(lfsv_p,j,lfsu_v,i, -phi_p[j] * div_phi_v[i] * incomp_scaling * weight);
              }
            } // end i

          } // end loop quadrature points
      } // end jacobian_volume

      // skeleton term, each face is only visited ONCE
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // dimensions
        const unsigned int dim = IG::dimension;
        const unsigned int dimw = IG::coorddimension;

        // subspaces
        using namespace Indices;
        using LFSV_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_v_s = child(lfsv_s,_0);
        const auto& lfsu_v_s = child(lfsu_s,_0);
        const auto& lfsv_v_n = child(lfsv_n,_0);
        const auto& lfsu_v_n = child(lfsu_n,_0);

        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_p_s = child(lfsv_s,_1);
        const auto& lfsu_p_s = child(lfsu_s,_1);
        const auto& lfsv_p_n = child(lfsv_n,_1);
        const auto& lfsu_p_n = child(lfsu_n,_1);

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using VectorBasisSwitch_V = VectorBasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;
        using BasisSwitch_P = BasisInterfaceSwitch<typename FESwitch_P::Basis >;
        using DF = typename BasisSwitch_V::DomainField;
        using RF = typename BasisSwitch_V::RangeField;
        using Range_V = typename BasisSwitch_V::Range;
        using Range_P = typename BasisSwitch_P::Range;
        using size_type = typename LFSV::Traits::SizeType;

        // References to inside and outside cells
        const auto& cell_inside = ig.inside();
        const auto& cell_outside = ig.outside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();
        auto geo_outside = cell_outside.geometry();

        // Get geometry of intersection in local coordinates of cell_inside and cell_outside
        auto geo_in_inside = ig.geometryInInside();
        auto geo_in_outside = ig.geometryInOutside();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v_s.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-2);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        const int epsilon = prm.epsilonIPSymmetryFactor();
        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        auto penalty_factor = prm.getFaceIP(geo,geo_inside,geo_outside);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : quadratureRule(geo,qorder))
          {

            // position of quadrature point in local coordinates of element
            auto local_s = geo_in_inside.global(ip.position());
            auto local_n = geo_in_outside.global(ip.position());

            // values of velocity shape functions
            std::vector<Range_V> phi_v_s(lfsv_v_s.size());
            std::vector<Range_V> phi_v_n(lfsv_v_n.size());
            FESwitch_V::basis(lfsv_v_s.finiteElement()).evaluateFunction(local_s,phi_v_s);
            FESwitch_V::basis(lfsv_v_n.finiteElement()).evaluateFunction(local_n,phi_v_n);

            // values of pressure shape functions
            std::vector<Range_P> phi_p_s(lfsv_p_s.size());
            std::vector<Range_P> phi_p_n(lfsv_p_n.size());
            FESwitch_P::basis(lfsv_p_s.finiteElement()).evaluateFunction(local_s,phi_p_s);
            FESwitch_P::basis(lfsv_p_n.finiteElement()).evaluateFunction(local_n,phi_p_n);

            // compute pressure value
            Range_P val_p_s(0.0);
            Range_P val_p_n(0.0);
            for (size_type i=0; i<lfsu_p_s.size(); i++)
              val_p_s.axpy(x_s(lfsu_p_s,i),phi_p_s[i]);
            for (size_type i=0; i<lfsu_p_n.size(); i++)
              val_p_n.axpy(x_n(lfsu_p_n,i),phi_p_n[i]);

            // evaluate jacobian of velocity shape functions on reference element
            std::vector<Dune::FieldMatrix<RF,dim,dim> > jac_phi_v_s(lfsu_v_s.size());
            std::vector<Dune::FieldMatrix<RF,dim,dim> > jac_phi_v_n(lfsu_v_n.size());
            VectorBasisSwitch_V::jacobian
              (FESwitch_V::basis(lfsv_v_s.finiteElement()), geo_inside, local_s, jac_phi_v_s);
            VectorBasisSwitch_V::jacobian
              (FESwitch_V::basis(lfsv_v_n.finiteElement()), geo_outside, local_n, jac_phi_v_n);

            // compute velocity value, jacobian, and divergence
            Range_V val_u_s(0.0);
            Range_V val_u_n(0.0);
            Dune::FieldMatrix<RF,dim,dim> jac_u_s(0.0);
            Dune::FieldMatrix<RF,dim,dim> jac_u_n(0.0);
            for (size_type i=0; i<lfsu_v_s.size(); i++){
              val_u_s.axpy(x_s(lfsu_v_s,i),phi_v_s[i]);
              jac_u_s.axpy(x_s(lfsu_v_s,i),jac_phi_v_s[i]);
            }
            for (size_type i=0; i<lfsu_v_n.size(); i++){
              val_u_n.axpy(x_n(lfsu_v_n,i),phi_v_n[i]);
              jac_u_n.axpy(x_n(lfsu_v_n,i),jac_phi_v_n[i]);
            }

            auto normal = ig.unitOuterNormal(ip.position());
            auto weight = ip.weight()*geo.integrationElement(ip.position());
            auto mu = prm.mu(ig,ip.position());

            auto factor = mu * weight;

            // compute jump in velocity
            auto jump = val_u_s - val_u_n;

            // compute mean in pressure
            auto mean_p = 0.5*(val_p_s + val_p_n);

            // compute flux of velocity jacobian
            Dune::FieldVector<DF,dimw> flux_jac_u(0.0);
            add_compute_flux(jac_u_s,normal,flux_jac_u);
            add_compute_flux(jac_u_n,normal,flux_jac_u);
            flux_jac_u *= 0.5;

            // loop over test functions, same element
            for (size_t i=0; i<lfsv_v_s.size(); i++) {
              //================================================//
              // diffusion term
              //================================================//
              r_s.accumulate(lfsv_v_s, i, -(flux_jac_u * phi_v_s[i]) * factor);

              //================================================//
              // (non-)symmetric IP term
              //================================================//
              Dune::FieldVector<DF,dimw> flux_jac_phi(0.0);
              add_compute_flux(jac_phi_v_s[i],normal,flux_jac_phi);
              r_s.accumulate(lfsv_v_s, i, epsilon * 0.5 * (flux_jac_phi * jump) * factor);

              //================================================//
              // standard IP term integral
              //================================================//
              r_s.accumulate(lfsv_v_s,i, penalty_factor * (jump*phi_v_s[i]) * factor);

              //================================================//
              // pressure-velocity-coupling in momentum equation
              //================================================//
              r_s.accumulate(lfsv_v_s,i, mean_p * (phi_v_s[i]*normal) * weight);
            }

            // loop over test functions, neighbour element
            for (size_t i=0; i<lfsv_v_n.size(); i++) {
              //================================================//
              // diffusion term
              //================================================//
              r_n.accumulate(lfsv_v_n, i,  (flux_jac_u * phi_v_n[i]) * factor);

              //================================================//
              // (non-)symmetric IP term
              //================================================//
              Dune::FieldVector<DF,dimw> flux_jac_phi(0.0);
              add_compute_flux(jac_phi_v_n[i],normal,flux_jac_phi);
              r_n.accumulate(lfsv_v_n, i, epsilon * 0.5 * (flux_jac_phi * jump) * factor);

              //================================================//
              // standard IP term integral
              //================================================//
              r_n.accumulate(lfsv_v_n,i, -penalty_factor * (jump*phi_v_n[i]) * factor);

              //================================================//
              // pressure-velocity-coupling in momentum equation
              //================================================//
              r_n.accumulate(lfsv_v_n,i, -mean_p * (phi_v_n[i]*normal) * weight);
            }

            //================================================//
            // incompressibility constraint
            //================================================//
            for (size_t i=0; i<lfsv_p_s.size(); i++)
              r_s.accumulate(lfsv_p_s,i, 0.5*phi_p_s[i] * (jump*normal) * incomp_scaling * weight);
            for (size_t i=0; i<lfsv_p_n.size(); i++)
              r_n.accumulate(lfsv_p_n,i, 0.5*phi_p_n[i] * (jump*normal) * incomp_scaling * weight);

          } // end loop quadrature points
      } // end alpha_skeleton

      // jacobian of skeleton term, each face is only visited ONCE
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X&, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X&, const LFSV& lfsv_n,
                              LocalMatrix& mat_ss, LocalMatrix& mat_sn,
                              LocalMatrix& mat_ns, LocalMatrix& mat_nn) const
      {
        // dimensions
        const unsigned int dim = IG::dimension;
        const unsigned int dimw = IG::coorddimension;

        // subspaces
        using namespace Indices;
        using LFSV_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_v_s = child(lfsv_s,_0);
        const auto& lfsu_v_s = child(lfsu_s,_0);
        const auto& lfsv_v_n = child(lfsv_n,_0);
        const auto& lfsu_v_n = child(lfsu_n,_0);

        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_p_s = child(lfsv_s,_1);
        const auto& lfsu_p_s = child(lfsu_s,_1);
        const auto& lfsv_p_n = child(lfsv_n,_1);
        const auto& lfsu_p_n = child(lfsu_n,_1);

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using VectorBasisSwitch_V = VectorBasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;
        using BasisSwitch_P = BasisInterfaceSwitch<typename FESwitch_P::Basis >;
        using DF = typename BasisSwitch_V::DomainField;
        using RF = typename BasisSwitch_V::RangeField;
        using Range_V = typename BasisSwitch_V::Range;
        using Range_P = typename BasisSwitch_P::Range;
        using size_type = typename LFSV::Traits::SizeType;

        // References to inside and outside cells
        auto const& cell_inside = ig.inside();
        auto const& cell_outside = ig.outside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();
        auto geo_outside = cell_outside.geometry();

        // Get geometry of intersection in local coordinates of cell_inside and cell_outside
        auto geo_in_inside = ig.geometryInInside();
        auto geo_in_outside = ig.geometryInOutside();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v_s.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-2);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        auto epsilon = prm.epsilonIPSymmetryFactor();
        auto incomp_scaling = prm.incompressibilityScaling(current_dt);

        auto penalty_factor = prm.getFaceIP(geo,geo_inside,geo_outside);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : quadratureRule(geo,qorder))
          {

            // position of quadrature point in local coordinates of element
            auto local_s = geo_in_inside.global(ip.position());
            auto local_n = geo_in_outside.global(ip.position());

            // values of velocity shape functions
            std::vector<Range_V> phi_v_s(lfsv_v_s.size());
            std::vector<Range_V> phi_v_n(lfsv_v_n.size());
            FESwitch_V::basis(lfsv_v_s.finiteElement()).evaluateFunction(local_s,phi_v_s);
            FESwitch_V::basis(lfsv_v_n.finiteElement()).evaluateFunction(local_n,phi_v_n);

            // values of pressure shape functions
            std::vector<Range_P> phi_p_s(lfsv_p_s.size());
            std::vector<Range_P> phi_p_n(lfsv_p_n.size());
            FESwitch_P::basis(lfsv_p_s.finiteElement()).evaluateFunction(local_s,phi_p_s);
            FESwitch_P::basis(lfsv_p_n.finiteElement()).evaluateFunction(local_n,phi_p_n);

            // evaluate jacobian of velocity shape functions on reference element
            std::vector<Dune::FieldMatrix<RF,dim,dim> > jac_phi_v_s(lfsu_v_s.size());
            std::vector<Dune::FieldMatrix<RF,dim,dim> > jac_phi_v_n(lfsu_v_n.size());
            VectorBasisSwitch_V::jacobian
              (FESwitch_V::basis(lfsv_v_s.finiteElement()), geo_inside, local_s, jac_phi_v_s);
            VectorBasisSwitch_V::jacobian
              (FESwitch_V::basis(lfsv_v_n.finiteElement()), geo_outside, local_n, jac_phi_v_n);

            auto normal = ig.unitOuterNormal(ip.position());
            auto weight = ip.weight()*geo.integrationElement(ip.position());
            auto mu = prm.mu(ig,ip.position());

            auto factor = mu * weight;

            //============================================
            // loop over test functions, same element
            //============================================
            for(size_type i=0; i<lfsv_v_s.size(); i++) {

              // compute flux
              Dune::FieldVector<DF,dimw> flux_jac_phi_i(0.0);
              add_compute_flux(jac_phi_v_s[i],normal,flux_jac_phi_i);

              //============================================
              // diffusion
              // (non-)symmetric IP-Term
              // standard IP integral
              //============================================
              for(size_type j=0; j<lfsu_v_s.size(); j++) {
                Dune::FieldVector<DF,dimw> flux_jac_phi_j(0.0);
                add_compute_flux(jac_phi_v_s[j],normal,flux_jac_phi_j);

                mat_ss.accumulate(lfsv_v_s,i,lfsu_v_s,j, -0.5 * (flux_jac_phi_j*phi_v_s[i]) * factor);
                mat_ss.accumulate(lfsv_v_s,i,lfsu_v_s,j, epsilon * 0.5 * (flux_jac_phi_i*phi_v_s[j]) * factor);
                mat_ss.accumulate(lfsv_v_s,i,lfsu_v_s,j, penalty_factor * (phi_v_s[j]*phi_v_s[i]) * factor);
              }

              for(size_type j=0; j<lfsu_v_n.size(); j++) {
                Dune::FieldVector<DF,dimw> flux_jac_phi_j(0.0);
                add_compute_flux(jac_phi_v_n[j],normal,flux_jac_phi_j);

                mat_sn.accumulate(lfsv_v_s,i,lfsu_v_n,j, -0.5 * (flux_jac_phi_j*phi_v_s[i]) * factor);
                mat_sn.accumulate(lfsv_v_s,i,lfsu_v_n,j, -epsilon * 0.5 * (flux_jac_phi_i*phi_v_n[j]) * factor);
                mat_sn.accumulate(lfsv_v_s,i,lfsu_v_n,j, -penalty_factor * (phi_v_n[j]*phi_v_s[i]) * factor);
              }

              //============================================
              // pressure-velocity coupling in momentum equation
              //============================================
              for(size_type j=0; j<lfsu_p_s.size(); j++) {
                mat_ss.accumulate(lfsv_v_s,i,lfsu_p_s,j, 0.5*phi_p_s[j] * (phi_v_s[i]*normal) * weight);
              }

              for(size_type j=0; j<lfsu_p_n.size(); j++) {
                mat_sn.accumulate(lfsv_v_s,i,lfsu_p_n,j, 0.5*phi_p_n[j] * (phi_v_s[i]*normal) * weight);
              }
            } // end i (same)

            //============================================
            // loop over test functions, neighbour element
            //============================================
            for(size_type i=0; i<lfsv_v_n.size(); i++) {

              // compute flux
              Dune::FieldVector<DF,dimw> flux_jac_phi_i(0.0);
              add_compute_flux(jac_phi_v_n[i],normal,flux_jac_phi_i);

              //============================================
              // diffusion
              // (non-)symmetric IP-Term
              // standard IP integral
              //============================================
              for(size_type j=0; j<lfsu_v_s.size(); j++) {
                Dune::FieldVector<DF,dimw> flux_jac_phi_j(0.0);
                add_compute_flux(jac_phi_v_s[j],normal,flux_jac_phi_j);

                mat_ns.accumulate(lfsv_v_n,i,lfsu_v_s,j, 0.5 * (flux_jac_phi_j*phi_v_n[i]) * factor);
                mat_ns.accumulate(lfsv_v_n,i,lfsu_v_s,j, epsilon * 0.5 * (flux_jac_phi_i*phi_v_s[j]) * factor);
                mat_ns.accumulate(lfsv_v_n,i,lfsu_v_s,j, -penalty_factor * (phi_v_s[j]*phi_v_n[i]) * factor);
              }

              for(size_type j=0; j<lfsu_v_n.size(); j++) {
                Dune::FieldVector<DF,dimw> flux_jac_phi_j(0.0);
                add_compute_flux(jac_phi_v_n[j],normal,flux_jac_phi_j);

                mat_nn.accumulate(lfsv_v_n,i,lfsu_v_n,j, 0.5 * (flux_jac_phi_j*phi_v_n[i]) * factor);
                mat_nn.accumulate(lfsv_v_n,i,lfsu_v_n,j, -epsilon * 0.5 * (flux_jac_phi_i*phi_v_n[j]) * factor);
                mat_nn.accumulate(lfsv_v_n,i,lfsu_v_n,j, penalty_factor * (phi_v_n[j]*phi_v_n[i]) * factor);
              }

              //============================================
              // pressure-velocity coupling in momentum equation
              //============================================
              for(size_type j=0; j<lfsu_p_s.size(); j++) {
                mat_ns.accumulate(lfsv_v_n,i,lfsu_p_s,j, -0.5*phi_p_s[j] * (phi_v_n[i]*normal) * weight);
              }

              for(size_type j=0; j<lfsu_p_n.size(); j++) {
                mat_nn.accumulate(lfsv_v_n,i,lfsu_p_n,j, -0.5*phi_p_n[j] * (phi_v_n[i]*normal) * weight);
              }
            } // end i (neighbour)

            //================================================//
            // \int <q> [u] n
            //================================================//
            for(size_type i=0; i<lfsv_p_s.size(); i++) {
              for(size_type j=0; j<lfsu_v_s.size(); j++)
                mat_ss.accumulate(lfsv_p_s,i,lfsu_v_s,j, 0.5*phi_p_s[i] * (phi_v_s[j]*normal) * incomp_scaling * weight);

              for(size_type j=0; j<lfsu_v_n.size(); j++)
                mat_sn.accumulate(lfsv_p_s,i,lfsu_v_n,j, -0.5*phi_p_s[i] * (phi_v_n[j]*normal) * incomp_scaling * weight);
            }

            for(size_type i=0; i<lfsv_p_n.size(); i++) {
              for(size_type j=0; j<lfsu_v_s.size(); j++)
                mat_ns.accumulate(lfsv_p_n,i,lfsu_v_s,j, 0.5*phi_p_n[i] * (phi_v_s[j]*normal) * incomp_scaling * weight);

              for(size_type j=0; j<lfsu_v_n.size(); j++)
                mat_nn.accumulate(lfsv_p_n,i,lfsu_v_n,j, -0.5*phi_p_n[i] * (phi_v_n[j]*normal) * incomp_scaling * weight);
            }

          } // end loop quadrature points
      } // end jacobian_skeleton

      // boundary term
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu, const X& x, const LFSV& lfsv,
                           R& r) const
      {
        // dimensions
        const unsigned int dim = IG::dimension;
        const unsigned int dimw = IG::coorddimension;

        // subspaces
        using namespace Indices;
        using LFSV_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_v = child(lfsv,_0);
        const auto& lfsu_v = child(lfsu,_0);

        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_p = child(lfsv,_1);
        const auto& lfsu_p = child(lfsu,_1);

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using VectorBasisSwitch_V = VectorBasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;
        using BasisSwitch_P = BasisInterfaceSwitch<typename FESwitch_P::Basis >;
        using DF = typename BasisSwitch_V::DomainField;
        using RF = typename BasisSwitch_V::RangeField;
        using Range_V = typename BasisSwitch_V::Range;
        using Range_P = typename BasisSwitch_P::Range;
        using size_type = typename LFSV::Traits::SizeType;

        // References to inside cell
        const auto& cell_inside = ig.inside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // Get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        auto epsilon = prm.epsilonIPSymmetryFactor();
        auto incomp_scaling = prm.incompressibilityScaling(current_dt);

        auto penalty_factor = prm.getFaceIP(geo,geo_inside);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            // position of quadrature point in local coordinates of element
            auto local = geo_in_inside.global(ip.position());

            // values of velocity shape functions
            std::vector<Range_V> phi_v(lfsv_v.size());
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

            // values of pressure shape functions
            std::vector<Range_P> phi_p(lfsv_p.size());
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            // evaluate jacobian of basis functions on reference element
            std::vector<Dune::FieldMatrix<RF,dim,dim> > jac_phi_v(lfsu_v.size());
            VectorBasisSwitch_V::jacobian
              (FESwitch_V::basis(lfsv_v.finiteElement()), geo_inside, local, jac_phi_v);

            // compute pressure value
            Range_P val_p(0.0);
            for (size_type i=0; i<lfsu_p.size(); i++)
              val_p.axpy(x(lfsu_p,i),phi_p[i]);

            // compute u and velocity jacobian
            Range_V val_u(0.0);
            Dune::FieldMatrix<RF,dim,dim> jac_u(0.0);
            for (size_type i=0; i<lfsu_v.size(); i++){
              val_u.axpy(x(lfsu_v,i),phi_v[i]);
              jac_u.axpy(x(lfsu_v,i),jac_phi_v[i]);
            }

            auto normal = ig.unitOuterNormal(ip.position());
            auto weight = ip.weight()*geo.integrationElement(ip.position());
            auto mu = prm.mu(ig,ip.position());

            // evaluate boundary condition type
            auto bctype(prm.bctype(ig,ip.position()));

            if (bctype == BC::VelocityDirichlet) {
              // compute jump relative to Dirichlet value
              auto u0(prm.g(cell_inside,local));
              auto jump = val_u - u0;

              // compute flux of velocity jacobian
              Dune::FieldVector<DF,dimw> flux_jac_u(0.0);
              add_compute_flux(jac_u,normal,flux_jac_u);

              for (size_t i=0; i<lfsv_v.size(); i++) {
                //================================================//
                // diffusion term
                //================================================//
                r.accumulate(lfsv_v,i, -mu * (flux_jac_u * phi_v[i]) * weight);

                //================================================//
                // (non-)symmetric IP term
                //================================================//
                Dune::FieldVector<DF,dimw> flux_jac_phi(0.0);
                add_compute_flux(jac_phi_v[i],normal,flux_jac_phi);
                r.accumulate(lfsv_v,i, mu * epsilon * (flux_jac_phi*jump) * weight);

                //================================================//
                // standard IP term integral
                //================================================//
                r.accumulate(lfsv_v,i, mu * (jump*phi_v[i]) * penalty_factor * weight);

                //================================================//
                // pressure-velocity-coupling in momentum equation
                //================================================//
                r.accumulate(lfsv_v,i, val_p * (phi_v[i]*normal) * weight);
              } // end i

              //================================================//
              // incompressibility constraint
              //================================================//
              for(size_type i=0; i<lfsv_p.size(); i++) {
                r.accumulate(lfsv_p,i, phi_p[i] * (jump*normal) * incomp_scaling * weight);
              }
            } // Velocity Dirichlet

            if (bctype == BC::StressNeumann) {
              auto stress(prm.j(ig,ip.position(),normal));

              for(size_type i=0; i<lfsv_v.size(); i++) {
                r.accumulate(lfsv_v,i, (stress*phi_v[i]) * weight);
              }
            } // Pressure Dirichlet

          } // end loop quadrature points
      } // end alpha_boundary

      // jacobian of boundary term
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu, const X& x, const LFSV& lfsv,
                              LocalMatrix& mat) const
      {
        // dimensions
        const unsigned int dim = IG::dimension;
        const unsigned int dimw = IG::coorddimension;

        // subspaces
        using namespace Indices;
        using LFSV_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_v = child(lfsv,_0);
        const auto& lfsu_v = child(lfsu,_0);

        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_p = child(lfsv,_1);
        const auto& lfsu_p = child(lfsu,_1);

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using VectorBasisSwitch_V = VectorBasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;
        using BasisSwitch_P = BasisInterfaceSwitch<typename FESwitch_P::Basis >;
        using DF = typename BasisSwitch_V::DomainField;
        using RF = typename BasisSwitch_V::RangeField;
        using Range_V = typename BasisSwitch_V::Range;
        using Range_P = typename BasisSwitch_P::Range;
        using size_type = typename LFSV::Traits::SizeType;

        // References to inside cell
        const auto& cell_inside = ig.inside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // Get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        auto epsilon = prm.epsilonIPSymmetryFactor();
        auto incomp_scaling = prm.incompressibilityScaling(current_dt);

        auto penalty_factor = prm.getFaceIP(geo,geo_inside);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            // position of quadrature point in local coordinates of element
            auto local = geo_in_inside.global(ip.position());

            // values of velocity shape functions
            std::vector<Range_V> phi_v(lfsv_v.size());
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

            // values of pressure shape functions
            std::vector<Range_P> phi_p(lfsv_p.size());
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            // evaluate jacobian of basis functions on reference element
            std::vector<Dune::FieldMatrix<RF,dim,dim> > jac_phi_v(lfsu_v.size());
            VectorBasisSwitch_V::jacobian
              (FESwitch_V::basis(lfsv_v.finiteElement()), geo_inside, local, jac_phi_v);

            auto normal = ig.unitOuterNormal(ip.position());
            auto weight = ip.weight()*geo.integrationElement(ip.position());
            auto mu = prm.mu(ig,ip.position());

            // evaluate boundary condition type
            auto bctype(prm.bctype(ig,ip.position()));

            if (bctype == BC::VelocityDirichlet) {

              for(size_type i=0; i<lfsv_v.size(); i++) {
                // compute flux
                Dune::FieldVector<DF,dimw> flux_jac_phi_i(0.0);
                add_compute_flux(jac_phi_v[i],normal,flux_jac_phi_i);

                for(size_type j=0; j<lfsu_v.size(); j++) {
                  //================================================//
                  // diffusion term
                  // (non-)symmetric IP term
                  //================================================//
                  Dune::FieldVector<DF,dimw> flux_jac_phi_j(0.0);
                  add_compute_flux(jac_phi_v[j],normal,flux_jac_phi_j);

                  mat.accumulate(lfsv_v,i,lfsu_v,j, -mu * (flux_jac_phi_j*phi_v[i]) * weight);
                  mat.accumulate(lfsv_v,i,lfsu_v,j, mu * epsilon * (flux_jac_phi_i*phi_v[j])  *weight);

                  //================================================//
                  // standard IP term integral
                  //================================================//
                  mat.accumulate(lfsv_v,i,lfsu_v,j, mu * (phi_v[j]*phi_v[i]) * penalty_factor * weight);
                }

                //================================================//
                // pressure-velocity-coupling in momentum equation
                //================================================//
                for(size_type j=0; j<lfsu_p.size(); j++) {
                  mat.accumulate(lfsv_v,i,lfsu_p,j, phi_p[j] * (phi_v[i]*normal) * weight);
                }
              } // end i

              //================================================//
              // incompressibility constraint
              //================================================//
              for(size_type i=0; i<lfsv_p.size(); i++) {
                for(size_type j=0; j<lfsu_v.size(); j++) {
                  mat.accumulate(lfsv_p,i,lfsu_v,j, phi_p[i] * (phi_v[j]*normal) * incomp_scaling * weight);
                }
              }

            } // Velocity Dirichlet

          } // end loop quadrature points
      } // end jacobian_boundary

      // volume integral depending only on test functions,
      // contains f on the right hand side
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        const unsigned int dim = EG::Geometry::mydimension;

        // subspaces
        using namespace Indices;
        using LFSV_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_v = child(lfsv,_0);

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using Range_V = typename BasisSwitch_V::Range;
        using size_type = typename LFSV::Traits::SizeType;

        // Get cell
        const auto& cell = eg.entity();

        // Get geometries
        auto geo = eg.geometry();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            auto local = ip.position();
            //const Dune::FieldVector<DF,dimw> global = eg.geometry().global(local);

            // values of velocity shape functions
            std::vector<Range_V> phi_v(lfsv_v.size());
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

            auto weight = ip.weight() * geo.integrationElement(ip.position());

            // evaluate source term
            auto fval(prm.f(cell,local));

            //================================================//
            // \int (f*v)
            //================================================//
            for(size_type i=0; i<lfsv_v.size(); i++)
              r.accumulate(lfsv_v,i, -(fval*phi_v[i]) * weight);

          } // end loop quadrature points
      } // end lambda_volume

    private :

      template<class M, class RF>
      static void contraction(const M & a, const M & b, RF & v)
      {
        v = 0;
        const int an = a.N();
        const int am = a.M();
        for(int r=0; r<an; ++r)
          for(int c=0; c<am; ++c)
            v += a[r][c] * b[r][c];
      }

      template<class DU, class R>
      static void add_compute_flux(const DU & du, const R & n, R & result)
      {
        const int N = du.N();
        const int M = du.M();
        for(int r=0; r<N; ++r)
          for(int c=0; c<M; ++c)
            result[r] += du[r][c] * n[c];
      }

      PRM& prm;
      const int superintegration_order;
      Real current_dt;
    }; // end class DGNavierStokesVelVecFEM

  } // end namespace PDELab
} // end namespace Dune
#endif // DUNE_PDELAB_LOCALOPERATOR_DGNAVIERSTOKESVELVECFEM_HH
