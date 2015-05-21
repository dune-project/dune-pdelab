// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALOPERATOR_DGNAVIERSTOKES_HH
#define DUNE_PDELAB_LOCALOPERATOR_DGNAVIERSTOKES_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/pdelab/localoperator/idefault.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/dgnavierstokesparameter.hh>

#ifndef VBLOCK
#define VBLOCK 0
#endif
#define PBLOCK (- VBLOCK + 1)

namespace Dune {
  namespace PDELab {

    /** \brief A local operator for solving the Navier-Stokes equations using a DG discretization

        \tparam PRM Parameter class for this local operator.

    */
    template<typename PRM>
    class DGNavierStokes :
      public LocalOperatorDefaultFlags,
      public FullSkeletonPattern, public FullVolumePattern,
      public InstationaryLocalOperatorDefaultMethods<double>
    {
      typedef StokesBoundaryCondition BC;
      typedef typename PRM::Traits::RangeField RF;

      typedef InstationaryLocalOperatorDefaultMethods<double> InstatBase;
      typedef typename InstatBase::RealType Real;

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
      DGNavierStokes (const PRM& _prm, int _superintegration_order=0) :
        prm(_prm), superintegration_order(_superintegration_order),
        current_dt(1.0)
      {}

      // Store current dt
      void preStep (RealType , RealType dt, int )
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
        const unsigned int dim = EG::Geometry::dimension;

        // subspaces
        static_assert
          ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");
        typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
        const LFSV_PFS_V& lfsv_pfs_v = lfsv.template child<VBLOCK>();
        static_assert
          ((LFSV_PFS_V::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");

        // ... we assume all velocity components are the same
        typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
        const LFSV_V& lfsv_v = lfsv_pfs_v.template child<0>();
        const unsigned int vsize = lfsv_v.size();
        typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
        const LFSV_P& lfsv_p = lfsv.template child<PBLOCK>();
        const unsigned int psize = lfsv_p.size();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::Range RT;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;
        typedef typename LFSV::Traits::SizeType size_type;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = gt.isSimplex() ? 0 : (dim-1);
        const int jac_order = gt.isSimplex() ? 0 : 1;
        const int qorder = 3*v_order - 1 + jac_order + det_jac_order + superintegration_order;
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            const Dune::FieldVector<DF,dim> local = ip.position();
            const RF mu = prm.mu(eg,local);
            const RF rho = prm.rho(eg,local);

            // compute u (if Navier term enabled)
            std::vector<RT> phi_v(vsize);
            Dune::FieldVector<RF,dim> vu(0.0);
            if(navier) {
              FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

              for(unsigned int d=0; d<dim; ++d) {
                const LFSV_V& lfsu_v = lfsv_pfs_v.child(d);
                for(size_type i=0; i<vsize; i++)
                  vu[d] += x(lfsu_v,i) * phi_v[i];
              }
            } // end navier

            // and value of pressure shape functions
            std::vector<RT> phi_p(psize);
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            // compute pressure
            RF p(0.0);
            for(size_type i=0; i<psize; i++)
              p += x(lfsv_p,i) * phi_p[i];

            // compute gradients
            std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                    eg.geometry(), local, grad_phi_v);

            // compute velocity jacobian
            Dune::FieldMatrix<RF,dim,dim> jacu(0.0);
            for(unsigned int d = 0; d<dim; ++d) {
              const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);
              for(size_type i=0; i<vsize; i++)
                jacu[d].axpy(x(lfsv_v,i), grad_phi_v[i][0]);
            }

            const RF detj = eg.geometry().integrationElement(ip.position());
            const RF weight = ip.weight() * detj;

            for(unsigned int d = 0; d<dim; ++d) {
              const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);

              for(size_type i=0; i<vsize; i++) {
                //================================================//
                // \int (mu*grad_u*grad_v)
                //================================================//
                r.accumulate(lfsv_v,i, mu * (jacu[d]*grad_phi_v[i][0]) * weight);

                // Assemble symmetric part for (grad u)^T
                if(full_tensor)
                  for(unsigned int dd = 0; dd<dim; ++dd)
                    r.accumulate(lfsv_v,i, mu * jacu[dd][d] * grad_phi_v[i][dd] * weight);

                //================================================//
                // \int -p \nabla\cdot v
                //================================================//
                r.accumulate(lfsv_v,i, -p * grad_phi_v[i][0][d] * weight);

                //================================================//
                // \int \rho ((u\cdot\nabla ) u )\cdot v
                //================================================//
                if(navier) {
                  // compute u * grad u_d
                  const RF u_nabla_u = vu * jacu[d];

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
        const unsigned int dim = EG::Geometry::dimension;

        // subspaces
        static_assert
          ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");
        typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
        const LFSV_PFS_V& lfsv_pfs_v = lfsv.template child<VBLOCK>();
        static_assert
          ((LFSV_PFS_V::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");

        // ... we assume all velocity components are the same
        typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
        const LFSV_V& lfsv_v = lfsv_pfs_v.template child<0>();
        const unsigned int vsize = lfsv_v.size();
        typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
        const LFSV_P& lfsv_p = lfsv.template child<PBLOCK>();
        const unsigned int psize = lfsv_p.size();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::Range RT;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;
        typedef typename LFSV::Traits::SizeType size_type;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = gt.isSimplex() ? 0 : (dim-1);
        const int jac_order = gt.isSimplex() ? 0 : 1;
        const int qorder = 3*v_order - 1 + jac_order + det_jac_order + superintegration_order;
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            const Dune::FieldVector<DF,dim> local = ip.position();
            const RF mu = prm.mu(eg,local);
            const RF rho = prm.rho(eg,local);

            // compute u (if Navier term enabled)
            std::vector<RT> phi_v(vsize);
            Dune::FieldVector<RF,dim> vu(0.0);
            if(navier) {
              FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

              for(unsigned int d=0; d<dim; ++d) {
                const LFSV_V& lfsu_v = lfsv_pfs_v.child(d);
                for(size_type i=0; i<vsize; i++)
                  vu[d] += x(lfsu_v,i) * phi_v[i];
              }
            } // end navier

            // and value of pressure shape functions
            std::vector<RT> phi_p(psize);
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            // compute gradients
            std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                    eg.geometry(), local, grad_phi_v);

            const RF detj = eg.geometry().integrationElement(ip.position());
            const RF weight = ip.weight() * detj;

            for(unsigned int dv = 0; dv<dim; ++dv) {
              const LFSV_V& lfsv_v = lfsv_pfs_v.child(dv);

              // gradient of dv-th velocity component
              Dune::FieldVector<RF,dim> gradu_dv(0.0);
              if(navier)
                for(size_type l=0; l<vsize; ++l)
                  gradu_dv.axpy(x(lfsv_v,l), grad_phi_v[l]);

              for(size_type i=0; i<vsize; i++) {

                for(size_type j=0; j<vsize; j++) {
                  //================================================//
                  // \int (mu*grad_u*grad_v)
                  //================================================//
                  mat.accumulate(lfsv_v,i,lfsv_v,j, mu * (grad_phi_v[j][0]*grad_phi_v[i][0]) * weight);

                  // Assemble symmetric part for (grad u)^T
                  if(full_tensor)
                    for(unsigned int du = 0; du<dim; ++du) {
                      const LFSV_V& lfsu_v = lfsv_pfs_v.child(du);
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
                    const LFSV_V& lfsu_v = lfsv_pfs_v.child(du);
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
        const unsigned int dim = IG::Geometry::dimension;
        const unsigned int dimw = IG::Geometry::dimensionworld;

        // subspaces
        static_assert
          ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");

        typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
        const LFSV_PFS_V& lfsv_s_pfs_v = lfsv_s.template child<VBLOCK>();
        const LFSV_PFS_V& lfsv_n_pfs_v = lfsv_n.template child<VBLOCK>();
        static_assert
          ((LFSV_PFS_V::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");

        // ... we assume all velocity components are the same
        typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
        const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.template child<0>();
        const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.template child<0>();
        const unsigned int vsize_s = lfsv_s_v.size();
        const unsigned int vsize_n = lfsv_n_v.size();
        typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
        const LFSV_P& lfsv_s_p = lfsv_s.template child<PBLOCK>();
        const LFSV_P& lfsv_n_p = lfsv_n.template child<PBLOCK>();
        const unsigned int psize_s = lfsv_s_p.size();
        const unsigned int psize_n = lfsv_n_p.size();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::Range RT;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;

        // make copy of inside and outside cell w.r.t. the intersection
        auto inside_cell = ig.inside();
        auto outside_cell = ig.outside();

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometry().type();
        const int v_order = FESwitch_V::basis(lfsv_s_v.finiteElement()).order();
        const int det_jac_order = gtface.isSimplex() ? 0 : (dim-2);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        const int epsilon = prm.epsilonIPSymmetryFactor();
        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : rule)
          {

            // position of quadrature point in local coordinates of element
            Dune::FieldVector<DF,dim> local_s = ig.geometryInInside().global(ip.position());
            Dune::FieldVector<DF,dim> local_n = ig.geometryInOutside().global(ip.position());

            const RF penalty_factor = prm.getFaceIP(ig,ip.position());

            // value of velocity shape functions
            std::vector<RT> phi_v_s(vsize_s);
            std::vector<RT> phi_v_n(vsize_n);
            FESwitch_V::basis(lfsv_s_v.finiteElement()).evaluateFunction(local_s,phi_v_s);
            FESwitch_V::basis(lfsv_n_v.finiteElement()).evaluateFunction(local_n,phi_v_n);

            // evaluate u
            assert(vsize_s == vsize_n);
            Dune::FieldVector<RF,dim> u_s(0.0), u_n(0.0);
            for(unsigned int d=0; d<dim; ++d) {
              const LFSV_V & lfsv_s_v = lfsv_s_pfs_v.child(d);
              const LFSV_V & lfsv_n_v = lfsv_n_pfs_v.child(d);
              for(unsigned int i=0; i<vsize_s; i++) {
                u_s[d] += x_s(lfsv_s_v,i) * phi_v_s[i];
                u_n[d] += x_n(lfsv_n_v,i) * phi_v_n[i];
              }
            }

            // value of pressure shape functions
            std::vector<RT> phi_p_s(psize_s);
            std::vector<RT> phi_p_n(psize_n);
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
            std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v_s(vsize_s);
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_s_v.finiteElement()),
                                    inside_cell.geometry(), local_s, grad_phi_v_s);

            std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v_n(vsize_n);
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_n_v.finiteElement()),
                                    outside_cell.geometry(), local_n, grad_phi_v_n);

            // evaluate velocity jacobian
            Dune::FieldMatrix<RF,dim,dim> jacu_s(0.0), jacu_n(0.0);
            for(unsigned int d=0; d<dim; ++d) {
              const LFSV_V & lfsv_s_v = lfsv_s_pfs_v.child(d);
              const LFSV_V & lfsv_n_v = lfsv_n_pfs_v.child(d);
              for(unsigned int i=0; i<vsize_s; i++) {
                jacu_s[d].axpy(x_s(lfsv_s_v,i), grad_phi_v_s[i][0]);
                jacu_n[d].axpy(x_n(lfsv_n_v,i), grad_phi_v_n[i][0]);
              }
            }

            const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(ip.position());
            const RF weight = ip.weight()*ig.geometry().integrationElement(ip.position());
            const RF mu = prm.mu(ig,ip.position());

            const RF factor = mu * weight;

            for(unsigned int d=0; d<dim; ++d) {
              const LFSV_V & lfsv_s_v = lfsv_s_pfs_v.child(d);
              const LFSV_V & lfsv_n_v = lfsv_n_pfs_v.child(d);

              //================================================//
              // diffusion term
              //================================================//
              RF val = 0.5 * ((jacu_s[d] * normal) + (jacu_n[d] * normal)) * factor;
              for(unsigned int i=0; i<vsize_s; i++) {
                r_s.accumulate(lfsv_s_v,i, -val * phi_v_s[i]);
                r_n.accumulate(lfsv_n_v,i, val * phi_v_n[i]);

                if(full_tensor) {
                  for(unsigned int dd=0; dd<dim; ++dd) {
                    RF Tval = 0.5 * (jacu_s[dd][d] + jacu_n[dd][d]) * normal[dd] * factor;
                    r_s.accumulate(lfsv_s_v,i, -Tval * phi_v_s[i]);
                    r_n.accumulate(lfsv_n_v,i, Tval * phi_v_n[i]);
                  }
                }
              } // end i

              //================================================//
              // (non-)symmetric IP term
              //================================================//
              RF jumpu_d  = u_s[d] - u_n[d];
              for(unsigned int i=0; i<vsize_s; i++) {
                r_s.accumulate(lfsv_s_v,i, epsilon * 0.5 * (grad_phi_v_s[i][0] * normal) * jumpu_d * factor);
                r_n.accumulate(lfsv_n_v,i, epsilon * 0.5 * (grad_phi_v_n[i][0] * normal) * jumpu_d * factor);

                if(full_tensor) {
                  for(unsigned int dd=0; dd<dim; ++dd) {
                    const LFSV_V& lfsv_s_v_dd = lfsv_s_pfs_v.child(dd);
                    const LFSV_V& lfsv_n_v_dd = lfsv_n_pfs_v.child(dd);

                    r_s.accumulate(lfsv_s_v_dd,i, epsilon * 0.5 * grad_phi_v_s[i][0][d] * normal[dd] * jumpu_d * factor);
                    r_n.accumulate(lfsv_n_v_dd,i, epsilon * 0.5 * grad_phi_v_n[i][0][d] * normal[dd] * jumpu_d * factor);
                  }
                }
              } // end i

              //================================================//
              // standard IP term integral
              //================================================//
              for(unsigned int i=0; i<vsize_s; i++) {
                r_s.accumulate(lfsv_s_v,i, penalty_factor * jumpu_d * phi_v_s[i] * weight);
                r_n.accumulate(lfsv_n_v,i, -penalty_factor * jumpu_d * phi_v_n[i] * weight);
              } // end i

              //================================================//
              // pressure-velocity-coupling in momentum equation
              //================================================//
              RF mean_p = 0.5*(p_s + p_n);
              for(unsigned int i=0; i<vsize_s; i++) {
                r_s.accumulate(lfsv_s_v,i, mean_p * phi_v_s[i] * normal[d] * weight);
                r_n.accumulate(lfsv_n_v,i, -mean_p * phi_v_n[i] * normal[d] * weight);
              } // end i
            } // end d

            //================================================//
            // incompressibility constraint
            //================================================//
            RF jumpu_n = (u_s*normal) - (u_n*normal);
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
        const unsigned int dim = IG::Geometry::dimension;
        const unsigned int dimw = IG::Geometry::dimensionworld;

        // subspaces
        static_assert
          ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");

        typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
        const LFSV_PFS_V& lfsv_s_pfs_v = lfsv_s.template child<VBLOCK>();
        const LFSV_PFS_V& lfsv_n_pfs_v = lfsv_n.template child<VBLOCK>();
        static_assert
          ((LFSV_PFS_V::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");

        // ... we assume all velocity components are the same
        typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
        const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.template child<0>();
        const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.template child<0>();
        const unsigned int vsize_s = lfsv_s_v.size();
        const unsigned int vsize_n = lfsv_n_v.size();
        typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
        const LFSV_P& lfsv_s_p = lfsv_s.template child<PBLOCK>();
        const LFSV_P& lfsv_n_p = lfsv_n.template child<PBLOCK>();
        const unsigned int psize_s = lfsv_s_p.size();
        const unsigned int psize_n = lfsv_n_p.size();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::Range RT;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;

        // make copy of inside and outside cell w.r.t. the intersection
        auto inside_cell = ig.inside();
        auto outside_cell = ig.outside();

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometry().type();
        const int v_order = FESwitch_V::basis(lfsv_s_v.finiteElement()).order();
        const int det_jac_order = gtface.isSimplex() ? 0 : (dim-2);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        const int epsilon = prm.epsilonIPSymmetryFactor();
        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : rule)
          {

            // position of quadrature point in local coordinates of element
            Dune::FieldVector<DF,dim> local_s = ig.geometryInInside().global(ip.position());
            Dune::FieldVector<DF,dim> local_n = ig.geometryInOutside().global(ip.position());

            const RF penalty_factor = prm.getFaceIP(ig,ip.position());

            // value of velocity shape functions
            std::vector<RT> phi_v_s(vsize_s);
            std::vector<RT> phi_v_n(vsize_n);
            FESwitch_V::basis(lfsv_s_v.finiteElement()).evaluateFunction(local_s,phi_v_s);
            FESwitch_V::basis(lfsv_n_v.finiteElement()).evaluateFunction(local_n,phi_v_n);
            // and value of pressure shape functions
            std::vector<RT> phi_p_s(psize_s);
            std::vector<RT> phi_p_n(psize_n);
            FESwitch_P::basis(lfsv_s_p.finiteElement()).evaluateFunction(local_s,phi_p_s);
            FESwitch_P::basis(lfsv_n_p.finiteElement()).evaluateFunction(local_n,phi_p_n);

            // compute gradients
            std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v_s(vsize_s);
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_s_v.finiteElement()),
                                    inside_cell.geometry(), local_s, grad_phi_v_s);

            std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v_n(vsize_n);
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_n_v.finiteElement()),
                                    outside_cell.geometry(), local_n, grad_phi_v_n);

            const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(ip.position());
            const RF weight = ip.weight()*ig.geometry().integrationElement(ip.position());
            const RF mu = prm.mu(ig,ip.position());

            assert(vsize_s == vsize_n);
            const RF factor = mu * weight;

            for(unsigned int d = 0; d < dim; ++d) {
              const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.child(d);
              const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.child(d);

              //================================================//
              // - (\mu \int < \nabla u > . normal . [v])
              // \mu \int \frac{\sigma}{|\gamma|^\beta} [u] \cdot [v]
              //================================================//
              for(unsigned int i=0; i<vsize_s; ++i) {

                for(unsigned int j=0; j<vsize_s; ++j) {
                  RF val = (0.5*(grad_phi_v_s[i][0]*normal)*phi_v_s[j]) * factor;
                  mat_ss.accumulate(lfsv_s_v,j,lfsv_s_v,i, -val);
                  mat_ss.accumulate(lfsv_s_v,i,lfsv_s_v,j, epsilon * val);
                  mat_ss.accumulate(lfsv_s_v,i,lfsv_s_v,j, phi_v_s[i] * phi_v_s[j] * penalty_factor * weight);

                  // Assemble symmetric part for (grad u)^T
                  if(full_tensor) {
                    for(unsigned int dd = 0; dd < dim; ++dd) {
                      RF Tval = (0.5*(grad_phi_v_s[i][0][d]*normal[dd])*phi_v_s[j]) * factor;
                      const LFSV_V& lfsv_s_v_dd = lfsv_s_pfs_v.child(dd);
                      mat_ss.accumulate(lfsv_s_v,j,lfsv_s_v_dd,i, - Tval);
                      mat_ss.accumulate(lfsv_s_v_dd,i,lfsv_s_v,j, epsilon*Tval );
                    }
                  }
                }

                for(unsigned int j=0; j<vsize_n; ++j) {
                  // the normal vector flipped, thus the sign flips
                  RF val = (-0.5*(grad_phi_v_s[i][0]*normal)*phi_v_n[j]) * factor;
                  mat_ns.accumulate(lfsv_n_v,j,lfsv_s_v,i,- val);
                  mat_sn.accumulate(lfsv_s_v,i,lfsv_n_v,j, epsilon*val);
                  mat_ns.accumulate(lfsv_n_v,j,lfsv_s_v,i, -phi_v_s[i] * phi_v_n[j] * penalty_factor * weight);

                  // Assemble symmetric part for (grad u)^T
                  if(full_tensor) {
                    for (unsigned int dd=0;dd<dim;++dd) {
                      RF Tval = (-0.5*(grad_phi_v_s[i][0][d]*normal[dd])*phi_v_n[j]) * factor;
                      const LFSV_V& lfsv_s_v_dd = lfsv_s_pfs_v.child(dd);
                      mat_ns.accumulate(lfsv_n_v,j,lfsv_s_v_dd,i,- Tval);
                      mat_sn.accumulate(lfsv_s_v_dd,i,lfsv_n_v,j, epsilon*Tval);
                    }
                  }
                }
              } // end i

              for(unsigned int i=0; i<vsize_n; ++i) {

                for(unsigned int j=0; j<vsize_s; ++j) {
                  RF val = (0.5*(grad_phi_v_n[i][0]*normal)*phi_v_s[j]) * factor;
                  mat_sn.accumulate(lfsv_s_v,j,lfsv_n_v,i, - val);
                  mat_ns.accumulate(lfsv_n_v,i,lfsv_s_v,j, epsilon*val );
                  mat_sn.accumulate(lfsv_s_v,j,lfsv_n_v,i, -phi_v_n[i] * phi_v_s[j] * penalty_factor * weight);

                  // Assemble symmetric part for (grad u)^T
                  if(full_tensor) {
                    for (unsigned int dd=0;dd<dim;++dd) {
                      RF Tval = (0.5*(grad_phi_v_n[i][0][d]*normal[dd])*phi_v_s[j]) * factor;
                      const LFSV_V& lfsv_n_v_dd = lfsv_n_pfs_v.child(dd);
                      mat_sn.accumulate(lfsv_s_v,j,lfsv_n_v_dd,i, - Tval);
                      mat_ns.accumulate(lfsv_n_v_dd,i,lfsv_s_v,j, epsilon*Tval );
                    }
                  }
                }

                for(unsigned int j=0; j<vsize_n; ++j) {
                  // the normal vector flipped, thus the sign flips
                  RF val = (-0.5*(grad_phi_v_n[i][0]*normal)*phi_v_n[j]) * factor;
                  mat_nn.accumulate(lfsv_n_v,j,lfsv_n_v,i, - val);
                  mat_nn.accumulate(lfsv_n_v,i,lfsv_n_v,j, epsilon*val);
                  mat_nn.accumulate(lfsv_n_v,j,lfsv_n_v,i, phi_v_n[i] * phi_v_n[j] * penalty_factor * weight);

                  // Assemble symmetric part for (grad u)^T
                  if(full_tensor) {
                    for (unsigned int dd=0;dd<dim;++dd) {
                      RF Tval = (-0.5*(grad_phi_v_n[i][0][d]*normal[dd])*phi_v_n[j]) * factor;
                      const LFSV_V& lfsv_n_v_dd = lfsv_n_pfs_v.child(dd);
                      mat_nn.accumulate(lfsv_n_v,j,lfsv_n_v_dd,i,- Tval);
                      mat_nn.accumulate(lfsv_n_v_dd,i,lfsv_n_v,j, epsilon*Tval);
                    }
                  }
                }
              } // end i

              //================================================//
              // \int <q> [u] n
              // \int <p> [v] n
              //================================================//
              for(unsigned int i=0; i<vsize_s; ++i) {

                for(unsigned int j=0; j<psize_s; ++j) {
                  RF val = 0.5*(phi_p_s[j]*normal[d]*phi_v_s[i]) * weight;
                  mat_ss.accumulate(lfsv_s_v,i,lfsv_s_p,j, val);
                  mat_ss.accumulate(lfsv_s_p,j,lfsv_s_v,i, val * incomp_scaling);
                }

                for(unsigned int j=0; j<psize_n; ++j) {
                  RF val = 0.5*(phi_p_n[j]*normal[d]*phi_v_s[i]) * weight;
                  mat_sn.accumulate(lfsv_s_v,i,lfsv_n_p,j, val);
                  mat_ns.accumulate(lfsv_n_p,j,lfsv_s_v,i, val * incomp_scaling);
                }
              } // end i

              for(unsigned int i=0; i<vsize_n; ++i) {

                for (unsigned int j=0; j<psize_s;++j) {
                  // the normal vector flipped, thus the sign flips
                  RF val = -0.5*(phi_p_s[j]*normal[d]*phi_v_n[i]) * weight;
                  mat_ns.accumulate(lfsv_n_v,i,lfsv_s_p,j, val);
                  mat_sn.accumulate(lfsv_s_p,j,lfsv_n_v,i, val * incomp_scaling);
                }

                for (unsigned int j=0; j<psize_n;++j) {
                  // the normal vector flipped, thus the sign flips
                  RF val = -0.5*(phi_p_n[j]*normal[d]*phi_v_n[i]) * weight;
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
        const unsigned int dim = IG::Geometry::dimension;
        const unsigned int dimw = IG::Geometry::dimensionworld;

        // subspaces
        static_assert
          ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");

        typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
        const LFSV_PFS_V& lfsv_pfs_v = lfsv.template child<VBLOCK>();
        static_assert
          ((LFSV_PFS_V::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");

        // ... we assume all velocity components are the same
        typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
        const LFSV_V& lfsv_v = lfsv_pfs_v.template child<0>();
        const unsigned int vsize = lfsv_v.size();
        typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
        const LFSV_P& lfsv_p = lfsv.template child<PBLOCK>();
        const unsigned int psize = lfsv_p.size();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::Range RT;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;

        // make copy of inside cell w.r.t. the boundary
        auto inside_cell = ig.inside();

        // select quadrature rule
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        Dune::GeometryType gtface = ig.geometry().type();
        const int det_jac_order = gtface.isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        const int epsilon = prm.epsilonIPSymmetryFactor();
        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : rule)
          {
            // position of quadrature point in local coordinates of element
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(ip.position());

            const RF penalty_factor = prm.getFaceIP(ig,ip.position() );

            // value of velocity shape functions
            std::vector<RT> phi_v(vsize);
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

            // evaluate u
            Dune::FieldVector<RF,dim> u(0.0);
            for(unsigned int d=0; d<dim; ++d) {
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
            std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                    inside_cell.geometry(), local, grad_phi_v);

            // evaluate velocity jacobian
            Dune::FieldMatrix<RF,dim,dim> jacu(0.0);
            for(unsigned int d=0; d<dim; ++d) {
              const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);
              for(unsigned int i=0; i<vsize; i++)
                jacu[d].axpy(x(lfsv_v,i), grad_phi_v[i][0]);
            }

            const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(ip.position());
            const RF weight = ip.weight()*ig.geometry().integrationElement(ip.position());
            const RF mu = prm.mu(ig,ip.position());

            // evaluate boundary condition type
            typename PRM::Traits::BoundaryCondition::Type bctype(prm.bctype(ig,ip.position()));

            // Slip factor smoothly switching between slip and no slip conditions.
            RF slip_factor = 0.0;
            typedef NavierStokesDGImp::VariableBoundarySlipSwitch<PRM> BoundarySlipSwitch;
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
                const RF factor = weight * (1.0 - slip_factor);

                for(unsigned int d = 0; d < dim; ++d) {
                  const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);

                  //======================
                  // TODO
                  // Put loops over i together!
                  //======================

                  //================================================//
                  // - (\mu \int \nabla u. normal . v)
                  //================================================//
                  RF val = (jacu[d] * normal) * factor * mu;
                  for(unsigned int i=0; i<vsize; i++) {
                    r.accumulate(lfsv_v,i, -val * phi_v[i]);
                    r.accumulate(lfsv_v,i, epsilon * mu * (grad_phi_v[i][0] * normal) * u[d] * factor);

                    //============================================
                    // TODO
                    // add contribution from tull tensor
                    //============================================

                  } // end i

                  //================================================//
                  // \mu \int \sigma / |\gamma|^\beta v u
                  //================================================//
                  for(unsigned int i=0; i<vsize; i++) {
                    r.accumulate(lfsv_v,i, u[d] * phi_v[i] * penalty_factor * factor);
                  } // end i

                  //================================================//
                  // \int p v n
                  //================================================//
                  for(unsigned int i=0; i<vsize; i++) {
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

            //============================================
            // TODO
            // At the moment I don't care about slip velocity
            // boundary conditions.
            //============================================
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
        const unsigned int dim = IG::Geometry::dimension;
        const unsigned int dimw = IG::Geometry::dimensionworld;

        // subspaces
        static_assert
          ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");

        typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
        const LFSV_PFS_V& lfsv_pfs_v = lfsv.template child<VBLOCK>();
        static_assert
          ((LFSV_PFS_V::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");

        // ... we assume all velocity components are the same
        typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
        const LFSV_V& lfsv_v = lfsv_pfs_v.template child<0>();
        const unsigned int vsize = lfsv_v.size();
        typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
        const LFSV_P& lfsv_p = lfsv.template child<PBLOCK>();
        const unsigned int psize = lfsv_p.size();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::Range RT;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;

        // make copy of inside cell w.r.t. the boundary
        auto inside_cell = ig.inside();

        // select quadrature rule
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        Dune::GeometryType gtface = ig.geometry().type();
        const int det_jac_order = gtface.isSimplex() ? 0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        const int epsilon = prm.epsilonIPSymmetryFactor();
        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : rule)
          {
            // position of quadrature point in local coordinates of element
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(ip.position());

            const RF penalty_factor = prm.getFaceIP(ig,ip.position() );

            // value of velocity shape functions
            std::vector<RT> phi_v(vsize);
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);
            // and value of pressure shape functions
            std::vector<RT> phi_p(psize);
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                    inside_cell.geometry(), local, grad_phi_v);

            const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(ip.position());
            const RF weight = ip.weight()*ig.geometry().integrationElement(ip.position());
            const RF mu = prm.mu(ig,ip.position());

            // evaluate boundary condition type
            typename PRM::Traits::BoundaryCondition::Type bctype(prm.bctype(ig,ip.position()));

            // Slip factor smoothly switching between slip and no slip conditions.
            RF slip_factor = 0.0;
            typedef NavierStokesDGImp::VariableBoundarySlipSwitch<PRM> BoundarySlipSwitch;
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
                const RF factor = weight * (1.0 - slip_factor);

                for(unsigned int d = 0; d < dim; ++d) {
                  const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);

                  for(unsigned int i=0; i<vsize; i++) {

                    for(unsigned int j=0; j<vsize; j++) {
                      //================================================//
                      // - (\mu \int \nabla u. normal . v)
                      //================================================//
                      RF val = ((grad_phi_v[j][0]*normal)*phi_v[i]) * factor * mu;
                      mat.accumulate(lfsv_v,i,lfsv_v,j, - val);
                      mat.accumulate(lfsv_v,j,lfsv_v,i, epsilon*val);

                      // Assemble symmetric part for (grad u)^T
                      if(full_tensor) {
                        for(unsigned int dd = 0; dd < dim; ++dd) {
                          const LFSV_V& lfsv_v_dd = lfsv_pfs_v.child(dd);
                          RF Tval = ((grad_phi_v[j][0][d]*normal[dd])*phi_v[i]) * factor * mu;
                          mat.accumulate(lfsv_v,i,lfsv_v_dd,j, - Tval);
                          mat.accumulate(lfsv_v_dd,j,lfsv_v,i, epsilon*Tval);
                        }
                      }
                      //================================================//
                      // \mu \int \sigma / |\gamma|^\beta v u
                      //================================================//
                      mat.accumulate(lfsv_v,j,lfsv_v,i, phi_v[i] * phi_v[j] * penalty_factor * factor);
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

            //============================================
            // TODO
            // At the moment I don't care about slip velocity
            // boundary conditions.
            //============================================

            // if (bctype == BC::SlipVelocity)
            //   {
            //     const RF factor = weight * (slip_factor);

            //     //================================================//
            //     // - (\mu \int \nabla u. normal . v)
            //     //================================================//

            //     for (unsigned int i=0;i<vsize;++i) // ansatz
            //       {
            //         for (unsigned int j=0;j<vsize;++j) // test
            //           {
            //             RF ten_sum = 1.0;

            //             // Assemble symmetric part for (grad u)^T
            //             if(full_tensor)
            //               ten_sum = 2.0;

            //             RF val = ten_sum * ((grad_phi_v[j][0]*normal)*phi_v[i]) * factor * mu;
            //             for (unsigned int d=0;d<dim;++d)
            //               {
            //                 const LFSV_V& lfsv_v_d = lfsv_pfs_v.child(d);

            //                 for (unsigned int dd=0;dd<dim;++dd)
            //                   {
            //                     const LFSV_V& lfsv_v_dd = lfsv_pfs_v.child(dd);

            //                     mat.accumulate(lfsv_v_dd,i,lfsv_v_d,j, -val*normal[d]*normal[dd]);
            //                     mat.accumulate(lfsv_v_d,j,lfsv_v_dd,i, epsilon*val*normal[d]*normal[dd]);
            //                   }
            //               }
            //           }
            //       }

            //     //================================================//
            //     // \mu \int \sigma / |\gamma|^\beta v u
            //     //================================================//
            //     const RF p_factor = penalty_factor * factor;
            //     for (unsigned int i=0;i<vsize;++i)
            //       {
            //         for (unsigned int j=0;j<vsize;++j)
            //           {
            //             RF val = phi_v[i]*phi_v[j] * p_factor;
            //             for (unsigned int d=0;d<dim;++d)
            //               {
            //                 const LFSV_V& lfsv_v_d = lfsv_pfs_v.child(d);
            //                 for (unsigned int dd=0;dd<dim;++dd)
            //                   {
            //                     const LFSV_V& lfsv_v_dd = lfsv_pfs_v.child(dd);
            //                     mat.accumulate(lfsv_v_d,j,lfsv_v_dd,i, val*normal[d]*normal[dd]);
            //                   }
            //               }
            //           }
            //       }

            //   } // Slip Velocity
          } // end loop quadrature points
      } // end jacobian_boundary

      // volume integral depending only on test functions,
      // contains f on the right hand side
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // dimensions
        static const unsigned int dim = EG::Geometry::dimension;

        // subspaces
        static_assert
          ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");

        typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
        const LFSV_PFS_V& lfsv_pfs_v = lfsv.template child<VBLOCK>();

        static_assert
          ((LFSV_PFS_V::CHILDREN == dim),"You seem to use the wrong function space for StokesDG");

        // we assume all velocity components are the same type
        typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
        const LFSV_V& lfsv_v = lfsv_pfs_v.template child<0>();
        const unsigned int vsize = lfsv_v.size();
        typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
        const LFSV_P& lfsv_p = lfsv.template child<PBLOCK>();
        const unsigned int psize = lfsv_p.size();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::Range RT;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;
        typedef typename LFSV::Traits::SizeType size_type;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = gt.isSimplex() ?  0 : (dim-1);
        // quad order is velocity order + det_jac order + superintegration
        const int qorder = v_order + det_jac_order + superintegration_order;

        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            const Dune::FieldVector<DF,dim> local = ip.position();
            //const Dune::FieldVector<DF,dimw> global = eg.geometry().global(local);

            // values of velocity shape functions
            std::vector<RT> phi_v(vsize);
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

            // values of pressure shape functions
            std::vector<RT> phi_p(psize);
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            const RF weight = ip.weight() * eg.geometry().integrationElement(ip.position());

            // evaluate source term
            typename PRM::Traits::VelocityRange fval(prm.f(eg,local));

            //================================================//
            // \int (f*v)
            //================================================//
            const RF factor = weight;
            for (unsigned int d=0; d<dim; d++) {
              const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);

              // and store for each velocity component
              for (size_type i=0; i<vsize; i++) {
                RF val = phi_v[i]*factor;
                r.accumulate(lfsv_v,i, -fval[d] * val);
              }
            }

            const RF g2 = prm.g2(eg,ip.position());

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
        static const unsigned int dim = IG::Geometry::dimension;

        // subspaces
        static_assert
          ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");

        typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
        const LFSV_PFS_V& lfsv_pfs_v = lfsv.template child<VBLOCK>();

        static_assert
          ((LFSV_PFS_V::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");

        // ... we assume all velocity components are the same
        typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
        const LFSV_V& lfsv_v = lfsv_pfs_v.template child<0>();
        const unsigned int vsize = lfsv_v.size();
        typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
        const LFSV_P& lfsv_p = lfsv.template child<PBLOCK>();
        const unsigned int psize = lfsv_p.size();

        // domain and range field type
        typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
        typedef BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
        typedef typename BasisSwitch_V::DomainField DF;
        typedef typename BasisSwitch_V::Range RT;
        typedef typename BasisSwitch_V::RangeField RF;
        typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;

        // make copy of inside cell w.r.t. the boundary
        auto inside_cell = ig.inside();

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometry().type();
        const int v_order = FESwitch_V::basis(lfsv_v.finiteElement()).order();
        const int det_jac_order = gtface.isSimplex() ? 0 : (dim-2);
        const int jac_order = gtface.isSimplex() ? 0 : 1;
        const int qorder = 2*v_order + det_jac_order + jac_order + superintegration_order;
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        const int epsilon = prm.epsilonIPSymmetryFactor();
        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : rule)
          {
            // position of quadrature point in local coordinates of element
            Dune::FieldVector<DF,dim-1> flocal = ip.position();
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(flocal);
            //Dune::FieldVector<DF,dimw> global = ig.geometry().global(flocal);

            const RF penalty_factor = prm.getFaceIP(ig,flocal);

            // value of velocity shape functions
            std::vector<RT> phi_v(vsize);
            FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);
            // and value of pressure shape functions
            std::vector<RT> phi_p(psize);
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                    inside_cell.geometry(), local, grad_phi_v);

            const Dune::FieldVector<DF,dim> normal = ig.unitOuterNormal(ip.position());
            const RF weight = ip.weight()*ig.geometry().integrationElement(ip.position());
            const RF mu = prm.mu(ig,flocal);

            // evaluate boundary condition type
            typename PRM::Traits::BoundaryCondition::Type bctype(prm.bctype(ig,flocal));

            if (bctype == BC::VelocityDirichlet)
              {
                typename PRM::Traits::VelocityRange u0(prm.g(inside_cell,local));

                RF factor = mu * weight;
                for(unsigned int d = 0; d < dim; ++d) {
                  const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);

                  for(unsigned int i=0; i<vsize; i++) {
                    //================================================//
                    // \mu \int \nabla v \cdot u_0 \cdot n
                    //================================================//
                    r.accumulate(lfsv_v,i, -epsilon * (grad_phi_v[i][0] * normal) * factor * u0[d]);

                    // Assemble symmetric part for (grad u)^T
                    if(full_tensor) {
                      for(unsigned int dd = 0; dd < dim; ++dd) {
                        const LFSV_V& lfsv_v_dd = lfsv_pfs_v.child(dd);
                        RF Tval = (grad_phi_v[i][0][d]*normal[dd]) * factor;
                        r.accumulate(lfsv_v_dd,i, -epsilon * Tval * u0[d]);
                      }
                    }
                    //================================================//
                    // \int \sigma / |\gamma|^\beta v u_0
                    //================================================//
                    r.accumulate(lfsv_v,i, -phi_v[i] * penalty_factor * u0[d] * weight);

                  } // end i
                } // end d

                //================================================//
                // \int q u_0 n
                //================================================//
                for (unsigned int i=0;i<psize;++i) // test
                  {
                    RF val = phi_p[i]*(u0 * normal) * weight;
                    r.accumulate(lfsv_p,i, - val * incomp_scaling);
                  }
              } // end BC velocity
            if (bctype == BC::StressNeumann)
              {
                typename PRM::Traits::VelocityRange stress(prm.j(ig,flocal,normal));

                //std::cout << "Pdirichlet\n";
                //================================================//
                // \int p u n
                //================================================//
                for(unsigned int d = 0; d < dim; ++d) {
                  const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);

                  for(unsigned int i=0; i<vsize; i++)
                    r.accumulate(lfsv_v,i, stress[d] * phi_v[i] * weight);
                }
              }
          } // end loop quadrature points
      } // end lambda_boundary

    private :
      const PRM& prm;
      const int superintegration_order;
      Real current_dt;
    }; // end class DGNavierStokes

  } // end namespace PDELab
} // end namespace Dune
#endif
