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
        const int qorder = 2*v_order - 2 + 2*jac_order + det_jac_order + superintegration_order;
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        const RF incomp_scaling = prm.incompressibilityScaling(current_dt);

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            const Dune::FieldVector<DF,dim> local = ip.position();
            const RF mu = prm.mu(eg,local);

            // and value of pressure shape functions
            std::vector<RT> phi_p(psize);
            FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

            // compute gradients
            std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);
            BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                    eg.geometry(), local, grad_phi_v);

            const RF detj = eg.geometry().integrationElement(ip.position());
            const RF weight = ip.weight() * detj;

            //================================================//
            // \int (mu*grad_u*grad_v)
            //================================================//
            const RF factor = mu * weight;
            for (size_type j=0; j<vsize; j++)
              {
                for (size_type i=0; i<vsize; i++)
                  {
                    // grad_phi_j*grad_phi_i
                    RF val = (grad_phi_v[j][0]*grad_phi_v[i][0])*factor;

                    for (unsigned int d=0; d<dim; d++)
                      {
                        const LFSV_V& lfsv_v_d = lfsv_pfs_v.child(d);
                        mat.accumulate(lfsv_v_d,i,lfsv_v_d,j, val);

                        // Assemble symmetric part for (grad u)^T
                        if(full_tensor){
                          for (unsigned int dd=0; dd<dim; dd++){
                            RF Tval = (grad_phi_v[j][0][d]*grad_phi_v[i][0][dd])*factor;
                            const LFSV_V& lfsv_v_dd = lfsv_pfs_v.child(dd);
                            mat.accumulate(lfsv_v_d,i,lfsv_v_dd,j, Tval);
                          }
                        }

                      }
                  }
              }

            //================================================//
            // - q * div u
            // - p * div v
            //================================================//
            for (size_type j=0; j<psize; j++) // test (q)
              {
                RF val = -1.0 * phi_p[j]*weight;
                for (size_type i=0; i<vsize; i++) // ansatz (u)
                  {
                    for (unsigned int d=0; d<dim; d++)
                      {
                        const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);
                        mat.accumulate(lfsv_p,j,lfsv_v,i, val*grad_phi_v[i][0][d] * incomp_scaling);
                        mat.accumulate(lfsv_v,i,lfsv_p,j, val*grad_phi_v[i][0][d]);
                      }
                  }
              }
          }
      }

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

            //================================================//
            // - (\mu \int < \nabla u > . normal . [v])
            //================================================//
            assert(vsize_s == vsize_n);
            RF factor = mu * weight;
            for (unsigned int i=0;i<vsize_s;++i)
              {
                for (unsigned int j=0;j<vsize_s;++j)
                  {
                    RF val = (0.5*(grad_phi_v_s[i][0]*normal)*phi_v_s[j]) * factor;
                    for (unsigned int d=0;d<dim;++d)
                      {
                        const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.child(d);
                        mat_ss.accumulate(lfsv_s_v,j,lfsv_s_v,i, - val);
                        mat_ss.accumulate(lfsv_s_v,i,lfsv_s_v,j, epsilon*val );

                        // Assemble symmetric part for (grad u)^T
                        if(full_tensor){

                          for (unsigned int dd=0;dd<dim;++dd)
                            {
                              RF Tval = (0.5*(grad_phi_v_s[i][0][d]*normal[dd])*phi_v_s[j]) * factor;
                              const LFSV_V& lfsv_s_v_dd = lfsv_s_pfs_v.child(dd);
                              mat_ss.accumulate(lfsv_s_v,j,lfsv_s_v_dd,i, - Tval);
                              mat_ss.accumulate(lfsv_s_v_dd,i,lfsv_s_v,j, epsilon*Tval );
                            }
                        }
                      }
                  }
                for (unsigned int j=0;j<vsize_n;++j)
                  {
                    // the normal vector flipped, thus the sign flips
                    RF val = (-0.5*(grad_phi_v_s[i][0]*normal)*phi_v_n[j]) * factor;
                    for (unsigned int d=0;d<dim;++d)
                      {
                        const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.child(d);
                        const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.child(d);
                        mat_ns.accumulate(lfsv_n_v,j,lfsv_s_v,i,- val);
                        mat_sn.accumulate(lfsv_s_v,i,lfsv_n_v,j, epsilon*val);

                        // Assemble symmetric part for (grad u)^T
                        if(full_tensor){

                          for (unsigned int dd=0;dd<dim;++dd)
                            {
                              RF Tval = (-0.5*(grad_phi_v_s[i][0][d]*normal[dd])*phi_v_n[j]) * factor;
                              const LFSV_V& lfsv_s_v_dd = lfsv_s_pfs_v.child(dd);
                              mat_ns.accumulate(lfsv_n_v,j,lfsv_s_v_dd,i,- Tval);
                              mat_sn.accumulate(lfsv_s_v_dd,i,lfsv_n_v,j, epsilon*Tval);
                            }
                        }
                      }
                  }
              }
            for (unsigned int i=0;i<vsize_n;++i)
              {
                for (unsigned int j=0;j<vsize_s;++j)
                  {
                    RF val = (0.5*(grad_phi_v_n[i][0]*normal)*phi_v_s[j]) * factor;
                    for (unsigned int d=0;d<dim;++d)
                      {
                        const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.child(d);
                        const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.child(d);
                        mat_sn.accumulate(lfsv_s_v,j,lfsv_n_v,i, - val);
                        mat_ns.accumulate(lfsv_n_v,i,lfsv_s_v,j, epsilon*val );

                        // Assemble symmetric part for (grad u)^T
                        if(full_tensor){

                          for (unsigned int dd=0;dd<dim;++dd)
                            {
                              RF Tval = (0.5*(grad_phi_v_n[i][0][d]*normal[dd])*phi_v_s[j]) * factor;
                              const LFSV_V& lfsv_n_v_dd = lfsv_n_pfs_v.child(dd);
                              mat_sn.accumulate(lfsv_s_v,j,lfsv_n_v_dd,i, - Tval);
                              mat_ns.accumulate(lfsv_n_v_dd,i,lfsv_s_v,j, epsilon*Tval );
                            }
                        }
                      }
                  }
                for (unsigned int j=0;j<vsize_n;++j)
                  {
                    // the normal vector flipped, thus the sign flips
                    RF val = (-0.5*(grad_phi_v_n[i][0]*normal)*phi_v_n[j]) * factor;
                    for (unsigned int d=0;d<dim;++d)
                      {
                        const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.child(d);
                        mat_nn.accumulate(lfsv_n_v,j,lfsv_n_v,i,- val);
                        mat_nn.accumulate(lfsv_n_v,i,lfsv_n_v,j, epsilon*val);

                        // Assemble symmetric part for (grad u)^T
                        if(full_tensor){

                          for (unsigned int dd=0;dd<dim;++dd)
                            {
                              RF Tval = (-0.5*(grad_phi_v_n[i][0][d]*normal[dd])*phi_v_n[j]) * factor;
                              const LFSV_V& lfsv_n_v_dd = lfsv_n_pfs_v.child(dd);
                              mat_nn.accumulate(lfsv_n_v,j,lfsv_n_v_dd,i,- Tval);
                              mat_nn.accumulate(lfsv_n_v_dd,i,lfsv_n_v,j, epsilon*Tval);
                            }
                        }
                      }
                  }
              }
            //================================================//
            // \mu \int \sigma / |\gamma|^\beta v u
            //================================================//
            factor = penalty_factor * weight;
            for (unsigned int i=0;i<vsize_s;++i)
              {
                for (unsigned int j=0;j<vsize_s;++j)
                  {
                    RF val = phi_v_s[i]*phi_v_s[j] * factor;
                    for (unsigned int d=0;d<dim;++d)
                      {
                        const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.child(d);
                        mat_ss.accumulate(lfsv_s_v,j,lfsv_s_v,i, val);
                      }
                  }
                for (unsigned int j=0;j<vsize_n;++j)
                  {
                    RF val = phi_v_s[i]*phi_v_n[j] * factor;
                    for (unsigned int d=0;d<dim;++d)
                      {
                        const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.child(d);
                        const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.child(d);
                        mat_ns.accumulate(lfsv_n_v,j,lfsv_s_v,i, - val);
                      }
                  }
              }
            for (unsigned int i=0;i<vsize_n;++i)
              {
                for (unsigned int j=0;j<vsize_s;++j)
                  {
                    RF val = phi_v_n[i]*phi_v_s[j] * factor;
                    for (unsigned int d=0;d<dim;++d)
                      {
                        const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.child(d);
                        const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.child(d);
                        mat_sn.accumulate(lfsv_s_v,j,lfsv_n_v,i, - val);
                      }
                  }
                for (unsigned int j=0;j<vsize_n;++j)
                  {
                    RF val = phi_v_n[i]*phi_v_n[j] * factor;
                    for (unsigned int d=0;d<dim;++d)
                      {
                        const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.child(d);
                        mat_nn.accumulate(lfsv_n_v,j,lfsv_n_v,i, val);
                      }
                  }
              }
            //================================================//
            // \int <q> [u] n
            // \int <p> [v] n
            //================================================//
            for (unsigned int i=0;i<vsize_s;++i)
              {
                for (unsigned int j=0;j<psize_s;++j)
                  {
                    for (unsigned int d=0;d<dim;++d)
                      {
                        const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.child(d);
                        RF val = 0.5*(phi_p_s[j]*normal[d]*phi_v_s[i]) * weight;
                        mat_ss.accumulate(lfsv_s_v,i,lfsv_s_p,j, val);
                        mat_ss.accumulate(lfsv_s_p,j,lfsv_s_v,i, val * incomp_scaling);
                      }
                  }
                for (unsigned int j=0;j<psize_n;++j)
                  {
                    for (unsigned int d=0;d<dim;++d)
                      {
                        const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.child(d);
                        RF val = 0.5*(phi_p_n[j]*normal[d]*phi_v_s[i]) * weight;
                        mat_sn.accumulate(lfsv_s_v,i,lfsv_n_p,j, val);
                        mat_ns.accumulate(lfsv_n_p,j,lfsv_s_v,i, val * incomp_scaling);
                      }
                  }
              }
            for (unsigned int i=0;i<vsize_n;++i)
              {
                for (unsigned int j=0;j<psize_s;++j)
                  {
                    for (unsigned int d=0;d<dim;++d)
                      {
                        const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.child(d);

                        // the normal vector flipped, thus the sign flips
                        RF val = -0.5*(phi_p_s[j]*normal[d]*phi_v_n[i]) * weight;
                        mat_ns.accumulate(lfsv_n_v,i,lfsv_s_p,j, val);
                        mat_sn.accumulate(lfsv_s_p,j,lfsv_n_v,i, val * incomp_scaling);
                      }
                  }
                for (unsigned int j=0;j<psize_n;++j)
                  {
                    for (unsigned int d=0;d<dim;++d)
                      {
                        const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.child(d);

                        // the normal vector flipped, thus the sign flips
                        RF val = -0.5*(phi_p_n[j]*normal[d]*phi_v_n[i]) * weight;
                        mat_nn.accumulate(lfsv_n_v,i,lfsv_n_p,j, val);
                        mat_nn.accumulate(lfsv_n_p,j,lfsv_n_v,i, val * incomp_scaling);
                      }
                  }
              }
          }
      }

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

        // evaluate boundary condition type
        typename PRM::Traits::BoundaryCondition::Type bctype(prm.bctype(ig,rule.begin()->position()));

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
                const RF factor = weight * (1.0 - slip_factor);

                //================================================//
                // - (\mu \int \nabla u. normal . v)
                //================================================//
                for (unsigned int i=0;i<vsize;++i) // ansatz
                  {
                    for (unsigned int j=0;j<vsize;++j) // test
                      {
                        RF val = ((grad_phi_v[j][0]*normal)*phi_v[i]) * factor * mu;
                        for (unsigned int d=0;d<dim;++d)
                          {
                            const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);
                            mat.accumulate(lfsv_v,i,lfsv_v,j, - val);
                            mat.accumulate(lfsv_v,j,lfsv_v,i, epsilon*val);

                            // Assemble symmetric part for (grad u)^T
                            if(full_tensor){

                              for (unsigned int dd=0;dd<dim;++dd)
                                {
                                  RF Tval = ((grad_phi_v[j][0][d]*normal[dd])*phi_v[i]) * factor * mu;
                                  const LFSV_V& lfsv_v_dd = lfsv_pfs_v.child(dd);
                                  mat.accumulate(lfsv_v,i,lfsv_v_dd,j, - Tval);
                                  mat.accumulate(lfsv_v_dd,j,lfsv_v,i, epsilon*Tval);
                                }
                            }
                          }
                      }
                  }
                //================================================//
                // \int q u n
                // \int p v n
                //================================================//
                for (unsigned int i=0;i<vsize;++i) // ansatz
                  {
                    for (unsigned int j=0;j<psize;++j) // test
                      {
                        for (unsigned int d=0;d<dim;++d)
                          {
                            const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);
                            RF val = (phi_p[j]*normal[d]*phi_v[i]) * weight;
                            mat.accumulate(lfsv_p,j,lfsv_v,i, val * incomp_scaling); // q u n
                            mat.accumulate(lfsv_v,i,lfsv_p,j, val); // p v n
                          }
                      }
                  }
                //================================================//
                // \mu \int \sigma / |\gamma|^\beta v u
                //================================================//
                const RF p_factor = penalty_factor * factor;
                for (unsigned int i=0;i<vsize;++i)
                  {
                    for (unsigned int j=0;j<vsize;++j)
                      {
                        RF val = phi_v[i]*phi_v[j] * p_factor;
                        for (unsigned int d=0;d<dim;++d)
                          {
                            const LFSV_V& lfsv_v = lfsv_pfs_v.child(d);
                            mat.accumulate(lfsv_v,j,lfsv_v,i, val);
                          }
                      }
                  }

              } // Velocity Dirichlet
            if (bctype == BC::SlipVelocity)
              {
                const RF factor = weight * (slip_factor);

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

                        RF val = ten_sum * ((grad_phi_v[j][0]*normal)*phi_v[i]) * factor * mu;
                        for (unsigned int d=0;d<dim;++d)
                          {
                            const LFSV_V& lfsv_v_d = lfsv_pfs_v.child(d);

                            for (unsigned int dd=0;dd<dim;++dd)
                              {
                                const LFSV_V& lfsv_v_dd = lfsv_pfs_v.child(dd);

                                mat.accumulate(lfsv_v_dd,i,lfsv_v_d,j, -val*normal[d]*normal[dd]);
                                mat.accumulate(lfsv_v_d,j,lfsv_v_dd,i, epsilon*val*normal[d]*normal[dd]);
                              }
                          }
                      }
                  }

                //================================================//
                // \mu \int \sigma / |\gamma|^\beta v u
                //================================================//
                const RF p_factor = penalty_factor * factor;
                for (unsigned int i=0;i<vsize;++i)
                  {
                    for (unsigned int j=0;j<vsize;++j)
                      {
                        RF val = phi_v[i]*phi_v[j] * p_factor;
                        for (unsigned int d=0;d<dim;++d)
                          {
                            const LFSV_V& lfsv_v_d = lfsv_pfs_v.child(d);
                            for (unsigned int dd=0;dd<dim;++dd)
                              {
                                const LFSV_V& lfsv_v_dd = lfsv_pfs_v.child(dd);
                                mat.accumulate(lfsv_v_d,j,lfsv_v_dd,i, val*normal[d]*normal[dd]);
                              }
                          }
                      }
                  }

              } // Slip Velocity
          }
      }

    private :
      const PRM& prm;
      const int superintegration_order;
      Real current_dt;
    }; // end class DGNavierStokes

  } // end namespace PDELab
} // end namespace Dune
#endif
