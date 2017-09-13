// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALOPERATOR_DGNAVIERSTOKES_HH
#define DUNE_PDELAB_LOCALOPERATOR_DGNAVIERSTOKES_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/pdelab/localoperator/idefault.hh>

#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/dgnavierstokesparameter.hh>
#include <dune/pdelab/localoperator/navierstokesmass.hh>

namespace Dune {
  namespace PDELab {

    /** \brief A local operator for solving the Navier-Stokes equations using a DG discretization

        \tparam PRM Parameter class for this local operator.

    */
    template<typename PRM>
    class DGNavierStokes :
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

    public:

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
      enum { doLambdaBoundary = true };

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
      DGNavierStokes (PRM& _prm, int _superintegration_order=0) :
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
        // dimensions
        const unsigned int dim = EG::Geometry::mydimension;

        // subspaces
        using namespace Indices;
        using LFSV_PFS_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_pfs_v = child(lfsv,_0);

        // ... we assume all velocity components are the same
        using LFSV_V = TypeTree::Child<LFSV_PFS_V,_0>;
        const auto& lfsv_v = child(lfsv_pfs_v,_0);
        const unsigned int vsize = lfsv_v.size();
        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_p = child(lfsv,_1);
        const unsigned int psize = lfsv_p.size();

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using RT = typename BasisSwitch_V::Range;
        using RF = typename BasisSwitch_V::RangeField;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;
        using size_type = typename LFSV::Traits::SizeType;

        // Get geometry
        auto geo = eg.geometry();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int jac_order = geo.type().isSimplex() ? 0 : 1;
        const int qorder = 3*v_order - 1 + jac_order + det_jac_order + superintegration_order;

        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // Initialize vectors outside for loop
        std::vector<RT> phi_v(vsize);
        Dune::FieldVector<RF,dim> vu(0.0);
        std::vector<RT> phi_p(psize);
        std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);
        Dune::FieldMatrix<RF,dim,dim> jacu(0.0);

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            auto local = ip.position();
            auto mu = prm.mu(eg,local);
            auto rho = prm.rho(eg,local);

            // compute u (if Navier term enabled)
            if(navier) {
              FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

              for(unsigned int d=0; d<dim; ++d) {
                vu[d] = 0.0;
                const auto& lfsu_v = lfsv_pfs_v.child(d);
                for(size_type i=0; i<vsize; i++)
                  vu[d] += x(lfsu_v,i) * phi_v[i];
              }
            } // end navier

            // and value of pressure shape functions
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            // compute pressure
            RF p(0.0);
            for(size_type i=0; i<psize; i++)
              p += x(lfsv_p,i) * phi_p[i];

            // compute gradients
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                    geo, local, grad_phi_v);

            // compute velocity jacobian
            for(unsigned int d = 0; d<dim; ++d) {
              jacu[d] = 0.0;
              const auto& lfsv_v = lfsv_pfs_v.child(d);
              for(size_type i=0; i<vsize; i++)
                jacu[d].axpy(x(lfsv_v,i), grad_phi_v[i][0]);
            }

            auto detj = geo.integrationElement(ip.position());
            auto weight = ip.weight() * detj;

            for(unsigned int d = 0; d<dim; ++d) {
              const auto& lfsv_v = lfsv_pfs_v.child(d);

              for(size_type i=0; i<vsize; i++) {
                //================================================//
                // \int (mu*grad_u*grad_v)
                //================================================//
                r.accumulate(lfsv_v,i, mu * (jacu[d]*grad_phi_v[i][0]) * weight);

                // Assemble symmetric part for (grad u)^T
                if(full_tensor)
                  for(unsigned int dd = 0; dd<dim; ++dd)
                    r.accumulate(lfsv_v,i, mu * jacu[dd][d] * grad_phi_v[i][0][dd] * weight);

                //================================================//
                // \int -p \nabla\cdot v
                //================================================//
                r.accumulate(lfsv_v,i, -p * grad_phi_v[i][0][d] * weight);

                //================================================//
                // \int \rho ((u\cdot\nabla ) u )\cdot v
                //================================================//
                if(navier) {
                  // compute u * grad u_d
                  auto u_nabla_u = vu * jacu[d];

                  r.accumulate(lfsv_v,i, rho * u_nabla_u * phi_v[i] * weight);
                } // end navier

              } // end i
            } // end d

            //================================================//
            // \int -q \nabla\cdot u
            //================================================//
            for(size_type i=0; i<psize; i++)
              for(unsigned int d = 0; d < dim; ++d)
                // divergence of u is the trace of the velocity jacobian
                r.accumulate(lfsv_p,i, -jacu[d][d] * phi_p[i] * incomp_scaling * weight);

          } // end loop quadrature points
      } // end alpha_volume

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            LocalMatrix& mat) const
      {
        // dimensions
        const unsigned int dim = EG::Geometry::mydimension;

        // subspaces
        using namespace Indices;
        using LFSV_PFS_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_pfs_v = child(lfsv,_0);

        // ... we assume all velocity components are the same
        using LFSV_V = TypeTree::Child<LFSV_PFS_V,_0>;
        const auto& lfsv_v = child(lfsv_pfs_v,_0);
        const unsigned int vsize = lfsv_v.size();
        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_p = child(lfsv,_1);
        const unsigned int psize = lfsv_p.size();

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using RT = typename BasisSwitch_V::Range;
        using RF = typename BasisSwitch_V::RangeField;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;
        using size_type = typename LFSV::Traits::SizeType;

         // Get geometry
        auto geo = eg.geometry();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int jac_order = geo.type().isSimplex() ? 0 : 1;
        const int qorder = 3*v_order - 1 + jac_order + det_jac_order + superintegration_order;

        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // Initialize vectors outside for loop
        std::vector<RT> phi_v(vsize);
        Dune::FieldVector<RF,dim> vu(0.0);
        std::vector<RT> phi_p(psize);
        std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);
        Dune::FieldVector<RF,dim> gradu_dv(0.0);

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            auto local = ip.position();
            auto mu = prm.mu(eg,local);
            auto rho = prm.rho(eg,local);

            // compute u (if Navier term enabled)
            if(navier) {
              FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

              for(unsigned int d=0; d<dim; ++d) {
                vu[d] = 0.0;
                const auto& lfsu_v = lfsv_pfs_v.child(d);
                for(size_type i=0; i<vsize; i++)
                  vu[d] += x(lfsu_v,i) * phi_v[i];
              }
            } // end navier

            // and value of pressure shape functions
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            // compute gradients
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                    geo, local, grad_phi_v);

            auto detj = geo.integrationElement(ip.position());
            auto weight = ip.weight() * detj;

            for(unsigned int dv = 0; dv<dim; ++dv) {
              const LFSV_V& lfsv_v = lfsv_pfs_v.child(dv);

              // gradient of dv-th velocity component
              gradu_dv = 0.0;
              if(navier)
                for(size_type l=0; l<vsize; ++l)
                  gradu_dv.axpy(x(lfsv_v,l), grad_phi_v[l][0]);

              for(size_type i=0; i<vsize; i++) {

                for(size_type j=0; j<vsize; j++) {
                  //================================================//
                  // \int (mu*grad_u*grad_v)
                  //================================================//
                  mat.accumulate(lfsv_v,i,lfsv_v,j, mu * (grad_phi_v[j][0]*grad_phi_v[i][0]) * weight);

                  // Assemble symmetric part for (grad u)^T
                  if(full_tensor)
                    for(unsigned int du = 0; du<dim; ++du) {
                      const auto& lfsu_v = lfsv_pfs_v.child(du);
                      mat.accumulate(lfsv_v,i,lfsu_v,j, mu * (grad_phi_v[j][0][dv]*grad_phi_v[i][0][du]) * weight);
                    }
                }

                //================================================//
                // - q * div u
                // - p * div v
                //================================================//
                for(size_type j=0; j<psize; j++) {
                  mat.accumulate(lfsv_p,j,lfsv_v,i, -phi_p[j] * grad_phi_v[i][0][dv] * incomp_scaling * weight);
                  mat.accumulate(lfsv_v,i,lfsv_p,j, -phi_p[j] * grad_phi_v[i][0][dv] * weight);
                }

                //================================================//
                // \int \rho ((u\cdot\nabla ) u )\cdot v
                //================================================//
                if(navier) {

                  // block diagonal contribution
                  for(size_type j=0; j<vsize; j++)
                    mat.accumulate(lfsv_v,i,lfsv_v,j, rho * (vu * grad_phi_v[j][0]) * phi_v[i] * weight);

                  // remaining contribution
                  for(unsigned int du = 0; du < dim; ++du) {
                    const auto& lfsu_v = lfsv_pfs_v.child(du);
                    for(size_type j=0; j<vsize; j++)
                      mat.accumulate(lfsv_v,i,lfsu_v,j, rho * phi_v[j] * gradu_dv[du] * phi_v[i] * weight);
                  }

                } // end navier

              } // end i
            } // end dv

          } // end loop quadrature points
      } // end jacobian_volume

      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // dimensions
        const unsigned int dim = IG::dimension;

        // subspaces
        using namespace Indices;
        using LFSV_PFS_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_s_pfs_v = child(lfsv_s,_0);
        const auto& lfsv_n_pfs_v = child(lfsv_n,_0);

        // ... we assume all velocity components are the same
        using LFSV_V = TypeTree::Child<LFSV_PFS_V,_0>;
        const auto& lfsv_s_v = child(lfsv_s_pfs_v,_0);
        const auto& lfsv_n_v = child(lfsv_n_pfs_v,_0);
        const unsigned int vsize_s = lfsv_s_v.size();
        const unsigned int vsize_n = lfsv_n_v.size();
        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_s_p = child(lfsv_s,_1);
        const auto& lfsv_n_p = child(lfsv_n,_1);
        const unsigned int psize_s = lfsv_s_p.size();
        const unsigned int psize_n = lfsv_n_p.size();

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using RT = typename BasisSwitch_V::Range;
        using RF = typename BasisSwitch_V::RangeField;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;

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
        const int v_order = FESwitch_V::basis(lfsv_s_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-2);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        const int epsilon = prm.epsilonIPSymmetryFactor();
        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // Initialize vectors outside for loop
        std::vector<RT> phi_v_s(vsize_s);
        std::vector<RT> phi_v_n(vsize_n);
        Dune::FieldVector<RF,dim> u_s(0.0), u_n(0.0);
        std::vector<RT> phi_p_s(psize_s);
        std::vector<RT> phi_p_n(psize_n);
        std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v_s(vsize_s);
        std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v_n(vsize_n);
        Dune::FieldMatrix<RF,dim,dim> jacu_s(0.0), jacu_n(0.0);

        auto penalty_factor = prm.getFaceIP(geo,geo_inside,geo_outside);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : quadratureRule(geo,qorder))
          {

            // position of quadrature point in local coordinates of element
            auto local_s = geo_in_inside.global(ip.position());
            auto local_n = geo_in_outside.global(ip.position());

            // value of velocity shape functions
            FESwitch_V::basis(lfsv_s_v.finiteElement()).evaluateFunction(local_s,phi_v_s);
            FESwitch_V::basis(lfsv_n_v.finiteElement()).evaluateFunction(local_n,phi_v_n);

            // evaluate u
            assert(vsize_s == vsize_n);
            for(unsigned int d=0; d<dim; ++d) {
              u_s[d] = 0.0;
              u_n[d] = 0.0;
              const auto& lfsv_s_v = lfsv_s_pfs_v.child(d);
              const auto& lfsv_n_v = lfsv_n_pfs_v.child(d);
              for(unsigned int i=0; i<vsize_s; i++) {
                u_s[d] += x_s(lfsv_s_v,i) * phi_v_s[i];
                u_n[d] += x_n(lfsv_n_v,i) * phi_v_n[i];
              }
            }

            // value of pressure shape functions
            FESwitch_P::basis(lfsv_s_p.finiteElement()).evaluateFunction(local_s,phi_p_s);
            FESwitch_P::basis(lfsv_n_p.finiteElement()).evaluateFunction(local_n,phi_p_n);

            // evaluate pressure
            assert(psize_s == psize_n);
            RF p_s(0.0), p_n(0.0);
            for(unsigned int i=0; i<psize_s; i++) {
              p_s += x_s(lfsv_s_p,i) * phi_p_s[i];
              p_n += x_n(lfsv_n_p,i) * phi_p_n[i];
            }

            // compute gradients
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_s_v.finiteElement()),
                                    geo_inside, local_s, grad_phi_v_s);

            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_n_v.finiteElement()),
                                    geo_outside, local_n, grad_phi_v_n);

            // evaluate velocity jacobian
            for(unsigned int d=0; d<dim; ++d) {
              jacu_s[d] = 0.0;
              jacu_n[d] = 0.0;
              const auto& lfsv_s_v = lfsv_s_pfs_v.child(d);
              const auto& lfsv_n_v = lfsv_n_pfs_v.child(d);
              for(unsigned int i=0; i<vsize_s; i++) {
                jacu_s[d].axpy(x_s(lfsv_s_v,i), grad_phi_v_s[i][0]);
                jacu_n[d].axpy(x_n(lfsv_n_v,i), grad_phi_v_n[i][0]);
              }
            }

            auto normal = ig.unitOuterNormal(ip.position());
            auto weight = ip.weight()*geo.integrationElement(ip.position());
            auto mu = prm.mu(ig,ip.position());

            auto factor = mu * weight;

            for(unsigned int d=0; d<dim; ++d) {
              const auto& lfsv_s_v = lfsv_s_pfs_v.child(d);
              const auto& lfsv_n_v = lfsv_n_pfs_v.child(d);

              //================================================//
              // diffusion term
              //================================================//
              auto val = 0.5 * ((jacu_s[d] * normal) + (jacu_n[d] * normal)) * factor;
              for(unsigned int i=0; i<vsize_s; i++) {
                r_s.accumulate(lfsv_s_v,i, -val * phi_v_s[i]);
                r_n.accumulate(lfsv_n_v,i, val * phi_v_n[i]);

                if(full_tensor) {
                  for(unsigned int dd=0; dd<dim; ++dd) {
                    auto Tval = 0.5 * (jacu_s[dd][d] + jacu_n[dd][d]) * normal[dd] * factor;
                    r_s.accumulate(lfsv_s_v,i, -Tval * phi_v_s[i]);
                    r_n.accumulate(lfsv_n_v,i, Tval * phi_v_n[i]);
                  }
                }
              } // end i

              //================================================//
              // (non-)symmetric IP term
              //================================================//
              auto jumpu_d  = u_s[d] - u_n[d];
              for(unsigned int i=0; i<vsize_s; i++) {
                r_s.accumulate(lfsv_s_v,i, epsilon * 0.5 * (grad_phi_v_s[i][0] * normal) * jumpu_d * factor);
                r_n.accumulate(lfsv_n_v,i, epsilon * 0.5 * (grad_phi_v_n[i][0] * normal) * jumpu_d * factor);

                if(full_tensor) {
                  for(unsigned int dd=0; dd<dim; ++dd) {
                    r_s.accumulate(lfsv_s_v,i, epsilon * 0.5 * grad_phi_v_s[i][0][dd] * (u_s[dd] - u_n[dd]) * normal[d] * factor);
                    r_n.accumulate(lfsv_n_v,i, epsilon * 0.5 * grad_phi_v_n[i][0][dd] * (u_s[dd] - u_n[dd]) * normal[d] * factor);
                  }
                }
              } // end i

              //================================================//
              // standard IP term integral
              //================================================//
              for(unsigned int i=0; i<vsize_s; i++) {
                r_s.accumulate(lfsv_s_v,i, penalty_factor * jumpu_d * phi_v_s[i] * factor);
                r_n.accumulate(lfsv_n_v,i, -penalty_factor * jumpu_d * phi_v_n[i] * factor);
              } // end i

              //================================================//
              // pressure-velocity-coupling in momentum equation
              //================================================//
              auto mean_p = 0.5*(p_s + p_n);
              for(unsigned int i=0; i<vsize_s; i++) {
                r_s.accumulate(lfsv_s_v,i, mean_p * phi_v_s[i] * normal[d] * weight);
                r_n.accumulate(lfsv_n_v,i, -mean_p * phi_v_n[i] * normal[d] * weight);
              } // end i
            } // end d

            //================================================//
            // incompressibility constraint
            //================================================//
            auto jumpu_n = (u_s*normal) - (u_n*normal);
            for(unsigned int i=0; i<psize_s; i++) {
              r_s.accumulate(lfsv_s_p,i, 0.5 * phi_p_s[i] * jumpu_n * incomp_scaling * weight);
              r_n.accumulate(lfsv_n_p,i, 0.5 * phi_p_n[i] * jumpu_n * incomp_scaling * weight);
            } // end i

          } // end loop quadrature points
      } // end alpha_skeleton

      // jacobian of skeleton term
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

        // subspaces
        using namespace Indices;
        using LFSV_PFS_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_s_pfs_v = child(lfsv_s,_0);
        const auto& lfsv_n_pfs_v = child(lfsv_n,_0);

        // ... we assume all velocity components are the same
        using LFSV_V = TypeTree::Child<LFSV_PFS_V,_0>;
        const auto& lfsv_s_v = child(lfsv_s_pfs_v,_0);
        const auto& lfsv_n_v = child(lfsv_n_pfs_v,_0);
        const unsigned int vsize_s = lfsv_s_v.size();
        const unsigned int vsize_n = lfsv_n_v.size();
        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_s_p = child(lfsv_s,_1);
        const auto& lfsv_n_p = child(lfsv_n,_1);
        const unsigned int psize_s = lfsv_s_p.size();
        const unsigned int psize_n = lfsv_n_p.size();

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using RT = typename BasisSwitch_V::Range;
        using RF = typename BasisSwitch_V::RangeField;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;

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
        const int v_order = FESwitch_V::basis(lfsv_s_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-2);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        const int epsilon = prm.epsilonIPSymmetryFactor();
        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // Initialize vectors outside for loop
        std::vector<RT> phi_v_s(vsize_s);
        std::vector<RT> phi_v_n(vsize_n);
        std::vector<RT> phi_p_s(psize_s);
        std::vector<RT> phi_p_n(psize_n);
        std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v_s(vsize_s);
        std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v_n(vsize_n);

        auto penalty_factor = prm.getFaceIP(geo,geo_inside,geo_outside);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : quadratureRule(geo,qorder))
          {

            // position of quadrature point in local coordinates of element
            auto local_s = geo_in_inside.global(ip.position());
            auto local_n = geo_in_outside.global(ip.position());

            // value of velocity shape functions
            FESwitch_V::basis(lfsv_s_v.finiteElement()).evaluateFunction(local_s,phi_v_s);
            FESwitch_V::basis(lfsv_n_v.finiteElement()).evaluateFunction(local_n,phi_v_n);
            // and value of pressure shape functions
            FESwitch_P::basis(lfsv_s_p.finiteElement()).evaluateFunction(local_s,phi_p_s);
            FESwitch_P::basis(lfsv_n_p.finiteElement()).evaluateFunction(local_n,phi_p_n);

            // compute gradients
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_s_v.finiteElement()),
                                    geo_inside, local_s, grad_phi_v_s);

            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_n_v.finiteElement()),
                                    geo_outside, local_n, grad_phi_v_n);

            auto normal = ig.unitOuterNormal(ip.position());
            auto weight = ip.weight()*geo.integrationElement(ip.position());
            auto mu = prm.mu(ig,ip.position());

            assert(vsize_s == vsize_n);
            auto factor = mu * weight;

            for(unsigned int d = 0; d < dim; ++d) {
              const auto& lfsv_s_v = lfsv_s_pfs_v.child(d);
              const auto& lfsv_n_v = lfsv_n_pfs_v.child(d);

              //================================================//
              // - (\mu \int < \nabla u > . normal . [v])
              // \mu \int \frac{\sigma}{|\gamma|^\beta} [u] \cdot [v]
              //================================================//
              for(unsigned int i=0; i<vsize_s; ++i) {

                for(unsigned int j=0; j<vsize_s; ++j) {
                  auto val = (0.5*(grad_phi_v_s[i][0]*normal)*phi_v_s[j]) * factor;
                  mat_ss.accumulate(lfsv_s_v,j,lfsv_s_v,i, -val);
                  mat_ss.accumulate(lfsv_s_v,i,lfsv_s_v,j, epsilon * val);
                  mat_ss.accumulate(lfsv_s_v,i,lfsv_s_v,j, phi_v_s[i] * phi_v_s[j] * penalty_factor * factor);

                  // Assemble symmetric part for (grad u)^T
                  if(full_tensor) {
                    for(unsigned int dd = 0; dd < dim; ++dd) {
                      auto Tval = (0.5*(grad_phi_v_s[i][0][d]*normal[dd])*phi_v_s[j]) * factor;
                      const auto& lfsv_s_v_dd = lfsv_s_pfs_v.child(dd);
                      mat_ss.accumulate(lfsv_s_v,j,lfsv_s_v_dd,i, - Tval);
                      mat_ss.accumulate(lfsv_s_v_dd,i,lfsv_s_v,j, epsilon*Tval );
                    }
                  }
                }

                for(unsigned int j=0; j<vsize_n; ++j) {
                  // the normal vector flipped, thus the sign flips
                  auto val = (-0.5*(grad_phi_v_s[i][0]*normal)*phi_v_n[j]) * factor;
                  mat_ns.accumulate(lfsv_n_v,j,lfsv_s_v,i,- val);
                  mat_sn.accumulate(lfsv_s_v,i,lfsv_n_v,j, epsilon*val);
                  mat_ns.accumulate(lfsv_n_v,j,lfsv_s_v,i, -phi_v_s[i] * phi_v_n[j] * penalty_factor * factor);

                  // Assemble symmetric part for (grad u)^T
                  if(full_tensor) {
                    for (unsigned int dd=0;dd<dim;++dd) {
                      auto Tval = (-0.5*(grad_phi_v_s[i][0][d]*normal[dd])*phi_v_n[j]) * factor;
                      const auto& lfsv_s_v_dd = lfsv_s_pfs_v.child(dd);
                      mat_ns.accumulate(lfsv_n_v,j,lfsv_s_v_dd,i,- Tval);
                      mat_sn.accumulate(lfsv_s_v_dd,i,lfsv_n_v,j, epsilon*Tval);
                    }
                  }
                }

                for(unsigned int j=0; j<vsize_s; ++j) {
                  auto val = (0.5*(grad_phi_v_n[i][0]*normal)*phi_v_s[j]) * factor;
                  mat_sn.accumulate(lfsv_s_v,j,lfsv_n_v,i, - val);
                  mat_ns.accumulate(lfsv_n_v,i,lfsv_s_v,j, epsilon*val );
                  mat_sn.accumulate(lfsv_s_v,j,lfsv_n_v,i, -phi_v_n[i] * phi_v_s[j] * penalty_factor * factor);

                  // Assemble symmetric part for (grad u)^T
                  if(full_tensor) {
                    for (unsigned int dd=0;dd<dim;++dd) {
                      auto Tval = (0.5*(grad_phi_v_n[i][0][d]*normal[dd])*phi_v_s[j]) * factor;
                      const auto& lfsv_n_v_dd = lfsv_n_pfs_v.child(dd);
                      mat_sn.accumulate(lfsv_s_v,j,lfsv_n_v_dd,i, - Tval);
                      mat_ns.accumulate(lfsv_n_v_dd,i,lfsv_s_v,j, epsilon*Tval );
                    }
                  }
                }

                for(unsigned int j=0; j<vsize_n; ++j) {
                  // the normal vector flipped, thus the sign flips
                  auto val = (-0.5*(grad_phi_v_n[i][0]*normal)*phi_v_n[j]) * factor;
                  mat_nn.accumulate(lfsv_n_v,j,lfsv_n_v,i, - val);
                  mat_nn.accumulate(lfsv_n_v,i,lfsv_n_v,j, epsilon*val);
                  mat_nn.accumulate(lfsv_n_v,j,lfsv_n_v,i, phi_v_n[i] * phi_v_n[j] * penalty_factor * factor);

                  // Assemble symmetric part for (grad u)^T
                  if(full_tensor) {
                    for (unsigned int dd=0;dd<dim;++dd) {
                      auto Tval = (-0.5*(grad_phi_v_n[i][0][d]*normal[dd])*phi_v_n[j]) * factor;
                      const auto& lfsv_n_v_dd = lfsv_n_pfs_v.child(dd);
                      mat_nn.accumulate(lfsv_n_v,j,lfsv_n_v_dd,i,- Tval);
                      mat_nn.accumulate(lfsv_n_v_dd,i,lfsv_n_v,j, epsilon*Tval);
                    }
                  }
                }

                //================================================//
                // \int <q> [u] n
                // \int <p> [v] n
                //================================================//
                for(unsigned int j=0; j<psize_s; ++j) {
                  auto val = 0.5*(phi_p_s[j]*normal[d]*phi_v_s[i]) * weight;
                  mat_ss.accumulate(lfsv_s_v,i,lfsv_s_p,j, val);
                  mat_ss.accumulate(lfsv_s_p,j,lfsv_s_v,i, val * incomp_scaling);
                }

                for(unsigned int j=0; j<psize_n; ++j) {
                  auto val = 0.5*(phi_p_n[j]*normal[d]*phi_v_s[i]) * weight;
                  mat_sn.accumulate(lfsv_s_v,i,lfsv_n_p,j, val);
                  mat_ns.accumulate(lfsv_n_p,j,lfsv_s_v,i, val * incomp_scaling);
                }

                for (unsigned int j=0; j<psize_s;++j) {
                  // the normal vector flipped, thus the sign flips
                  auto val = -0.5*(phi_p_s[j]*normal[d]*phi_v_n[i]) * weight;
                  mat_ns.accumulate(lfsv_n_v,i,lfsv_s_p,j, val);
                  mat_sn.accumulate(lfsv_s_p,j,lfsv_n_v,i, val * incomp_scaling);
                }

                for (unsigned int j=0; j<psize_n;++j) {
                  // the normal vector flipped, thus the sign flips
                  auto val = -0.5*(phi_p_n[j]*normal[d]*phi_v_n[i]) * weight;
                  mat_nn.accumulate(lfsv_n_v,i,lfsv_n_p,j, val);
                  mat_nn.accumulate(lfsv_n_p,j,lfsv_n_v,i, val * incomp_scaling);
                }
              } // end i
            } // end d

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

        // subspaces
        using namespace Indices;
        using LFSV_PFS_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_pfs_v = child(lfsv,_0);

        // ... we assume all velocity components are the same
        using LFSV_V = TypeTree::Child<LFSV_PFS_V,_0>;
        const auto& lfsv_v = child(lfsv_pfs_v,_0);
        const unsigned int vsize = lfsv_v.size();
        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_p = child(lfsv,_1);
        const unsigned int psize = lfsv_p.size();

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using RT = typename BasisSwitch_V::Range;
        using RF = typename BasisSwitch_V::RangeField;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;

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

        // Initialize vectors outside for loop
        std::vector<RT> phi_v(vsize);
        Dune::FieldVector<RF,dim> u(0.0);
        std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);
        Dune::FieldMatrix<RF,dim,dim> jacu(0.0);

        auto penalty_factor = prm.getFaceIP(geo,geo_inside);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            // position of quadrature point in local coordinates of element
            auto local = geo_in_inside.global(ip.position());

            // value of velocity shape functions
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

            // evaluate u
            for(unsigned int d=0; d<dim; ++d) {
              u[d] = 0.0;
              const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);
              for(unsigned int i=0; i<vsize; i++)
                u[d] += x(lfsv_v,i) * phi_v[i];
            }

            // value of pressure shape functions
            std::vector<RT> phi_p(psize);
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            // evaluate pressure
            RF p(0.0);
            for(unsigned int i=0; i<psize; i++)
              p += x(lfsv_p,i) * phi_p[i];

            // compute gradients
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                    geo_inside, local, grad_phi_v);

            // evaluate velocity jacobian
            for(unsigned int d=0; d<dim; ++d) {
              jacu[d] = 0.0;
              const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);
              for(unsigned int i=0; i<vsize; i++)
                jacu[d].axpy(x(lfsv_v,i), grad_phi_v[i][0]);
            }

            auto normal = ig.unitOuterNormal(ip.position());
            auto weight = ip.weight()*geo.integrationElement(ip.position());
            auto mu = prm.mu(ig,ip.position());

            // evaluate boundary condition type
            auto bctype(prm.bctype(ig,ip.position()));

            // Slip factor smoothly switching between slip and no slip conditions.
            RF slip_factor = 0.0;
            using BoundarySlipSwitch = NavierStokesDGImp::VariableBoundarySlipSwitch<PRM>;
            if (bctype == BC::SlipVelocity)
              // Calls boundarySlip(..) function of parameter
              // class if available, i.e. if
              // enable_variable_slip is defined. Otherwise
              // returns 1.0;
              slip_factor = BoundarySlipSwitch::boundarySlip(prm,ig,ip.position());

            // velocity boundary condition
            if (bctype == BC::VelocityDirichlet || bctype == BC::SlipVelocity)
              {
                // on BC::VelocityDirichlet: 1.0 - slip_factor = 1.0
                auto factor = weight * (1.0 - slip_factor);

                for(unsigned int d = 0; d < dim; ++d) {
                  const auto& lfsv_v = lfsv_pfs_v.child(d);

                  auto val = (jacu[d] * normal) * factor * mu;
                  for(unsigned int i=0; i<vsize; i++) {
                    //================================================//
                    // - (\mu \int \nabla u. normal . v)
                    //================================================//
                    r.accumulate(lfsv_v,i, -val * phi_v[i]);
                    r.accumulate(lfsv_v,i, epsilon * mu * (grad_phi_v[i][0] * normal) * u[d] * factor);

                    if(full_tensor) {
                      for(unsigned int dd=0; dd<dim; ++dd) {
                        r.accumulate(lfsv_v,i, -mu * jacu[dd][d] * normal[dd] * phi_v[i] * factor);
                        r.accumulate(lfsv_v,i, epsilon * mu * grad_phi_v[i][0][dd] * u[dd] * normal[d] * factor);
                      }
                    }

                    //================================================//
                    // \mu \int \sigma / |\gamma|^\beta v u
                    //================================================//
                    r.accumulate(lfsv_v,i, u[d] * phi_v[i] * mu * penalty_factor * factor);

                    //================================================//
                    // \int p v n
                    //================================================//
                    r.accumulate(lfsv_v,i, p * phi_v[i] * normal[d] * weight);
                  } // end i
                } // end d

                //================================================//
                // \int q u n
                //================================================//
                for(unsigned int i=0; i<psize; i++) {
                  r.accumulate(lfsv_p,i, phi_p[i] * (u * normal) * incomp_scaling * weight);
                }
              } // Velocity Dirichlet

            if(bctype == BC::SlipVelocity)
              {
                auto factor = weight * (slip_factor);

                RF ten_sum = 1.0;
                if(full_tensor)
                  ten_sum = 2.0;

                for(unsigned int d = 0; d < dim; ++d) {
                  const auto& lfsv_v = lfsv_pfs_v.child(d);

                  for(unsigned int i=0; i<vsize; i++) {
                    //================================================//
                    // - (\mu \int \nabla u. normal . v)
                    //================================================//
                    for(unsigned int dd = 0; dd < dim; ++dd)
                      r.accumulate(lfsv_v,i, -ten_sum * mu * (jacu[dd] * normal) * normal[dd] *phi_v[i] * normal[d] * factor);
                    r.accumulate(lfsv_v,i, epsilon * ten_sum * mu * (u * normal) * (grad_phi_v[i][0] * normal) * normal[d] * factor);

                    //================================================//
                    // \mu \int \sigma / |\gamma|^\beta v u
                    //================================================//
                    r.accumulate(lfsv_v,i, mu * penalty_factor * (u * normal) * phi_v[i] * normal[d] * factor);
                  } // end i
                } // end d

              } // Slip Velocity
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

        // subspaces
        using namespace Indices;
        using LFSV_PFS_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_pfs_v = child(lfsv,_0);

        // ... we assume all velocity components are the same
        using LFSV_V = TypeTree::Child<LFSV_PFS_V,_0>;
        const auto& lfsv_v = child(lfsv_pfs_v,_0);
        const unsigned int vsize = lfsv_v.size();
        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_p = child(lfsv,_1);
        const unsigned int psize = lfsv_p.size();

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using RT = typename BasisSwitch_V::Range;
        using RF = typename BasisSwitch_V::RangeField;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;

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

        // Initialize vectors outside for loop
        std::vector<RT> phi_v(vsize);
        std::vector<RT> phi_p(psize);
        std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);

        auto penalty_factor = prm.getFaceIP(geo,geo_inside);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            // position of quadrature point in local coordinates of element
            auto local = geo_in_inside.global(ip.position());

            // value of velocity shape functions
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);
            // and value of pressure shape functions
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                    geo_inside, local, grad_phi_v);

            auto normal = ig.unitOuterNormal(ip.position());
            auto weight = ip.weight()*geo.integrationElement(ip.position());
            auto mu = prm.mu(ig,ip.position());

            // evaluate boundary condition type
            auto bctype(prm.bctype(ig,ip.position()));

            // Slip factor smoothly switching between slip and no slip conditions.
            RF slip_factor = 0.0;
            using BoundarySlipSwitch = NavierStokesDGImp::VariableBoundarySlipSwitch<PRM>;
            if (bctype == BC::SlipVelocity)
              // Calls boundarySlip(..) function of parameter
              // class if available, i.e. if
              // enable_variable_slip is defined. Otherwise
              // returns 1.0;
              slip_factor = BoundarySlipSwitch::boundarySlip(prm,ig,ip.position());

            // velocity boundary condition
            if (bctype == BC::VelocityDirichlet || bctype == BC::SlipVelocity)
              {
                // on BC::VelocityDirichlet: 1.0 - slip_factor = 1.0
                auto factor = weight * (1.0 - slip_factor);

                for(unsigned int d = 0; d < dim; ++d) {
                  const auto& lfsv_v = lfsv_pfs_v.child(d);

                  for(unsigned int i=0; i<vsize; i++) {

                    for(unsigned int j=0; j<vsize; j++) {
                      //================================================//
                      // - (\mu \int \nabla u. normal . v)
                      //================================================//
                      auto val = ((grad_phi_v[j][0]*normal)*phi_v[i]) * factor * mu;
                      mat.accumulate(lfsv_v,i,lfsv_v,j, - val);
                      mat.accumulate(lfsv_v,j,lfsv_v,i, epsilon*val);

                      // Assemble symmetric part for (grad u)^T
                      if(full_tensor) {
                        for(unsigned int dd = 0; dd < dim; ++dd) {
                          const auto& lfsv_v_dd = lfsv_pfs_v.child(dd);
                          auto Tval = ((grad_phi_v[j][0][d]*normal[dd])*phi_v[i]) * factor * mu;
                          mat.accumulate(lfsv_v,i,lfsv_v_dd,j, - Tval);
                          mat.accumulate(lfsv_v_dd,j,lfsv_v,i, epsilon*Tval);
                        }
                      }
                      //================================================//
                      // \mu \int \sigma / |\gamma|^\beta v u
                      //================================================//
                      mat.accumulate(lfsv_v,j,lfsv_v,i, phi_v[i] * phi_v[j] * mu * penalty_factor * factor);
                    }

                    //================================================//
                    // \int q u n
                    // \int p v n
                    //================================================//
                    for(unsigned int j=0; j<psize; j++) {
                      mat.accumulate(lfsv_p,j,lfsv_v,i, (phi_p[j]*normal[d]*phi_v[i]) * weight * incomp_scaling);
                      mat.accumulate(lfsv_v,i,lfsv_p,j, (phi_p[j]*normal[d]*phi_v[i]) * weight);
                    }
                  } // end i
                } // end d

              } // Velocity Dirichlet

            if (bctype == BC::SlipVelocity)
              {
                auto factor = weight * (slip_factor);

                //================================================//
                // - (\mu \int \nabla u. normal . v)
                //================================================//

                for (unsigned int i=0;i<vsize;++i) // ansatz
                  {
                    for (unsigned int j=0;j<vsize;++j) // test
                      {
                        RF ten_sum = 1.0;

                        // Assemble symmetric part for (grad u)^T
                        if(full_tensor)
                          ten_sum = 2.0;

                        auto val = ten_sum * ((grad_phi_v[j][0]*normal)*phi_v[i]) * factor * mu;
                        for (unsigned int d=0;d<dim;++d)
                          {
                            const auto& lfsv_v_d = lfsv_pfs_v.child(d);

                            for (unsigned int dd=0;dd<dim;++dd)
                              {
                                const auto& lfsv_v_dd = lfsv_pfs_v.child(dd);

                                mat.accumulate(lfsv_v_dd,i,lfsv_v_d,j, -val*normal[d]*normal[dd]);
                                mat.accumulate(lfsv_v_d,j,lfsv_v_dd,i, epsilon*val*normal[d]*normal[dd]);
                              }
                          }
                      }
                  }

                //================================================//
                // \mu \int \sigma / |\gamma|^\beta v u
                //================================================//
                auto p_factor = mu * penalty_factor * factor;
                for (unsigned int i=0;i<vsize;++i)
                  {
                    for (unsigned int j=0;j<vsize;++j)
                      {
                        auto val = phi_v[i]*phi_v[j] * p_factor;
                        for (unsigned int d=0;d<dim;++d)
                          {
                            const auto& lfsv_v_d = lfsv_pfs_v.child(d);
                            for (unsigned int dd=0;dd<dim;++dd)
                              {
                                const auto& lfsv_v_dd = lfsv_pfs_v.child(dd);
                                mat.accumulate(lfsv_v_d,j,lfsv_v_dd,i, val*normal[d]*normal[dd]);
                              }
                          }
                      }
                  }

              } // Slip Velocity
          } // end loop quadrature points
      } // end jacobian_boundary

      // volume integral depending only on test functions,
      // contains f on the right hand side
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // dimensions
        static const unsigned int dim = EG::Geometry::mydimension;

        // subspaces
        using namespace Indices;
        using LFSV_PFS_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_pfs_v = child(lfsv,_0);

        // we assume all velocity components are the same type
        using LFSV_V = TypeTree::Child<LFSV_PFS_V,_0>;
        const auto& lfsv_v = child(lfsv_pfs_v,_0);
        const unsigned int vsize = lfsv_v.size();
        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_p = child(lfsv,_1);
        const unsigned int psize = lfsv_p.size();

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using RT = typename BasisSwitch_V::Range;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;
        using size_type = typename LFSV::Traits::SizeType;

        // Get cell
        const auto& cell = eg.entity();

        // Get geometries
        auto geo = eg.geometry();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ?  0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        // Initialize vectors outside for loop
        std::vector<RT> phi_v(vsize);
        std::vector<RT> phi_p(psize);

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            auto local = ip.position();
            //const Dune::FieldVector<DF,dimw> global = eg.geometry().global(local);

            // values of velocity shape functions
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

            // values of pressure shape functions
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            auto weight = ip.weight() * geo.integrationElement(ip.position());

            // evaluate source term
            auto fval(prm.f(cell,local));

            //================================================//
            // \int (f*v)
            //================================================//
            auto factor = weight;
            for (unsigned int d=0; d<dim; d++) {
              const auto& lfsv_v = lfsv_pfs_v.child(d);

              // and store for each velocity component
              for (size_type i=0; i<vsize; i++) {
                auto val = phi_v[i]*factor;
                r.accumulate(lfsv_v,i, -fval[d] * val);
              }
            }

            auto g2 = prm.g2(eg,ip.position());

            // integrate div u * psi_i
            for (size_t i=0; i<lfsv_p.size(); i++) {
              r.accumulate(lfsv_p,i, g2 * phi_p[i] * factor);
            }

          } // end loop quadrature points
      } // end lambda_volume

      // boundary integral independent of ansatz functions,
      // Neumann and Dirichlet boundary conditions, DG penalty term's right hand side
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
        // dimensions
        static const unsigned int dim = IG::dimension;

        // subspaces
        using namespace Indices;
        using LFSV_PFS_V = TypeTree::Child<LFSV,_0>;
        const auto& lfsv_pfs_v = child(lfsv,_0);

        // ... we assume all velocity components are the same
        using LFSV_V = TypeTree::Child<LFSV_PFS_V,_0>;
        const auto& lfsv_v = child(lfsv_pfs_v,_0);
        const unsigned int vsize = lfsv_v.size();
        using LFSV_P = TypeTree::Child<LFSV,_1>;
        const auto& lfsv_p = child(lfsv,_1);
        const unsigned int psize = lfsv_p.size();

        // domain and range field type
        using FESwitch_V = FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType >;
        using BasisSwitch_V = BasisInterfaceSwitch<typename FESwitch_V::Basis >;
        using DF = typename BasisSwitch_V::DomainField;
        using RT = typename BasisSwitch_V::Range;
        using RF = typename BasisSwitch_V::RangeField;
        using FESwitch_P = FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType >;

        // References to inside cell
        const auto& cell_inside = ig.inside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // Get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // Determine quadrature order
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-2);
        const int jac_order = geo.type().isSimplex() ? 0 : 1;
        const int qorder = 2*v_order + det_jac_order + jac_order + superintegration_order;

        auto epsilon = prm.epsilonIPSymmetryFactor();
        auto incomp_scaling = prm.incompressibilityScaling(current_dt);

        // Initialize vectors outside for loop
        std::vector<RT> phi_v(vsize);
        std::vector<RT> phi_p(psize);
        std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);

        auto penalty_factor = prm.getFaceIP(geo,geo_inside);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            // position of quadrature point in local coordinates of element
            Dune::FieldVector<DF,dim-1> flocal = ip.position();
            Dune::FieldVector<DF,dim> local = geo_in_inside.global(flocal);
            //Dune::FieldVector<DF,dimw> global = ig.geometry().global(flocal);

            // value of velocity shape functions
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);
            // and value of pressure shape functions
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            // compute gradients
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                    geo_inside, local, grad_phi_v);

            auto normal = ig.unitOuterNormal(ip.position());
            auto weight = ip.weight()*geo.integrationElement(ip.position());
            auto mu = prm.mu(ig,flocal);

            // evaluate boundary condition type
            auto bctype(prm.bctype(ig,flocal));

            if (bctype == BC::VelocityDirichlet)
              {
                auto u0(prm.g(cell_inside,local));

                auto factor = mu * weight;
                for(unsigned int d = 0; d < dim; ++d) {
                  const auto& lfsv_v = lfsv_pfs_v.child(d);

                  for(unsigned int i=0; i<vsize; i++) {
                    //================================================//
                    // \mu \int \nabla v \cdot u_0 \cdot n
                    //================================================//
                    r.accumulate(lfsv_v,i, -epsilon * (grad_phi_v[i][0] * normal) * factor * u0[d]);

                    // Assemble symmetric part for (grad u)^T
                    if(full_tensor) {
                      for(unsigned int dd = 0; dd < dim; ++dd) {
                        const auto& lfsv_v_dd = lfsv_pfs_v.child(dd);
                        auto Tval = (grad_phi_v[i][0][d]*normal[dd]) * factor;
                        r.accumulate(lfsv_v_dd,i, -epsilon * Tval * u0[d]);
                      }
                    }
                    //================================================//
                    // \int \sigma / |\gamma|^\beta v u_0
                    //================================================//
                    r.accumulate(lfsv_v,i, -phi_v[i] * penalty_factor * u0[d] * factor);

                  } // end i
                } // end d

                //================================================//
                // \int q u_0 n
                //================================================//
                for (unsigned int i=0;i<psize;++i) // test
                  {
                    auto val = phi_p[i]*(u0 * normal) * weight;
                    r.accumulate(lfsv_p,i, - val * incomp_scaling);
                  }
              } // end BC velocity
            if (bctype == BC::StressNeumann)
              {
                auto stress(prm.j(ig,flocal,normal));

                //std::cout << "Pdirichlet\n";
                //================================================//
                // \int p u n
                //================================================//
                for(unsigned int d = 0; d < dim; ++d) {
                  const auto& lfsv_v = lfsv_pfs_v.child(d);

                  for(unsigned int i=0; i<vsize; i++)
                    r.accumulate(lfsv_v,i, stress[d] * phi_v[i] * weight);
                }
              }
          } // end loop quadrature points
      } // end lambda_boundary

    private :
      PRM& prm;
      const int superintegration_order;
      Real current_dt;
    }; // end class DGNavierStokes

  } // end namespace PDELab
} // end namespace Dune
#endif // DUNE_PDELAB_LOCALOPERATOR_DGNAVIERSTOKES_HH
