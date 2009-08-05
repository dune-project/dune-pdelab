// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_DIFFUSIONDG_HH
#define DUNE_PDELAB_DIFFUSIONDG_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/grid/common/referenceelements.hh>

#include"../common/geometrywrapper.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"../gridoperatorspace/gridoperatorspaceutilities.hh"
#include"pattern.hh"
#include"flags.hh"

namespace Dune {
  namespace PDELab {

    // a local operator for solving the diffusion equation
    //     - div (K grad u) = 0 in \Omega,
    //                    u = g on \Gamma_D
    //    (- K grad u) * nu = j on \Gamma_N
    // discontinuous Galerkin method
    // G : grid function for Dirichlet boundary conditions
    template<typename K, typename F, typename B, typename G, typename J>
    class DiffusionDG :
      public LocalOperatorDefaultFlags,
      public FullSkeletonPattern, public FullVolumePattern
      //,public NumericalJacobianVolume<DiffusionDG<K, F, B, G, J> >
      //,public NumericalJacobianSkeleton<DiffusionDG<K, F, B, G, J> >
      //,public NumericalJacobianBoundary<DiffusionDG<K, F, B, G, J> >
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume    = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume   = true };
      enum { doLambdaSkeleton = false };
      enum { doLambdaBoundary = true };

      DiffusionDG (const K& k_, const F& f_, const B& b_, const G& g_, const J& j_, 
                   const int dg_method_, int qorder_=4) :
        k(k_), f(f_), b(b_), g(g_), j(j_), dg_method(dg_method_), qorder(qorder_)
      {
        // OBB
        if (dg_method == 0)
          {
            epsilon = 1.0;
            sigma = 0.0;
            beta = 0.0;
          }
        // NIPG
        else if (dg_method == 1)
          {
            epsilon = 1.0;
            sigma = 4.5;   // should be < 5 for stability reasons
            beta = 2.0 - 0.5*F::Traits::dimDomain;  // 2D => 1, 3D => 0.5
          }
        // SIPG
        else if (dg_method == 2)
          {
            epsilon = -1.0;
            sigma = 4.5;   // should be < 5 for stability reasons
            beta = 2.0 - 0.5*F::Traits::dimDomain;  // 2D => 1, 3D => 0.5
          }
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensionslocal
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        typename K::Traits::RangeType tensor(0.0);
        Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
        k.evaluate(eg.entity(),localcenter,tensor);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> js(lfsu.size());
            lfsu.localFiniteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradient to real element
            const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            for (size_t i=0; i<lfsu.size(); i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (size_t i=0; i<lfsu.size(); i++)
              gradu.axpy(x[i],gradphi[i]);

            // compute K * gradient of u
            Dune::FieldVector<RF,dim> Kgradu(0.0);
            tensor.umv(gradu,Kgradu);

            // integrate (K grad u)*grad phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_t i=0; i<lfsu.size(); i++)
              {
                r[i] += Kgradu*gradphi[i] * factor;
              }
          }
      }

      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;

        // dimensions
        const int dim = IG::dimension;
        const int dimw = IG::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        // normal of center in face's reference element
        const Dune::FieldVector<DF,IG::dimension-1>& face_center =
          Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(face_center);
            
        // evaluate diffusion tensor at elements' centers, assume they are constant over elements
        const Dune::FieldVector<DF,IG::dimension>& 
          inside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>& 
          outside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);
        typename K::Traits::RangeType permeability_s(0.0);
        typename K::Traits::RangeType permeability_n(0.0);
        k.evaluate(*(ig.inside()),inside_local,permeability_s);
        k.evaluate(*(ig.outside()),outside_local,permeability_n);

        /*for (unsigned int i = 0; i < K::Traits::GridViewType::dimension; ++i)
          {
          for (unsigned int j = 0; j < K::Traits::GridViewType::dimension; ++j)
          {
          if (permeability_s[i][j] * permeability_n[i][j] == 0.0)
          {
          permeability[i][j] = 0.0;
          }
          else
          {
          permeability[i][j] = 2.0 / (1.0/permeability_s[i][j] + 1.0/permeability_n[i][j]);
          }
          }
          }*/

        // penalty weight for NIPG / SIPG
        RF penalty_weight_s = sigma / pow(ig.inside()->geometry().volume(), beta);
        RF penalty_weight_n = sigma / pow(ig.outside()->geometry().volume(), beta);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // position of quadrature point in local coordinates of element
            Dune::FieldVector<DF,dim> local_s = ig.geometryInInside().global(it->position());
            Dune::FieldVector<DF,dim> local_n = ig.geometryInOutside().global(it->position());

            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> js_s(lfsv_s.size());
            lfsv_s.localFiniteElement().localBasis().evaluateJacobian(local_s,js_s);
            std::vector<JacobianType> js_n(lfsv_n.size());
            lfsv_n.localFiniteElement().localBasis().evaluateJacobian(local_n,js_n);

            // transform gradient to real element
            const Dune::FieldMatrix<DF,dimw,dim> jac_s = ig.inside()->geometry().jacobianInverseTransposed(local_s);
            std::vector<Dune::FieldVector<RF,dim> > gradphi_s(lfsv_s.size());
            for (size_t i=0; i<lfsv_s.size(); i++)
              {
                gradphi_s[i] = 0.0;
                jac_s.umv(js_s[i][0],gradphi_s[i]);
              }
            const Dune::FieldMatrix<DF,dimw,dim> jac_n = ig.outside()->geometry().jacobianInverseTransposed(local_n);
            std::vector<Dune::FieldVector<RF,dim> > gradphi_n(lfsv_n.size());
            for (size_t i=0; i<lfsv_n.size(); i++)
              {
                gradphi_n[i] = 0.0;
                jac_n.umv(js_n[i][0],gradphi_n[i]);
              }

            // evaluate test shape functions
            std::vector<RangeType> phi_s(lfsv_s.size());
            lfsv_s.localFiniteElement().localBasis().evaluateFunction(local_s,phi_s);
            std::vector<RangeType> phi_n(lfsv_n.size());
            lfsv_n.localFiniteElement().localBasis().evaluateFunction(local_n,phi_n);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu_s(0.0);
            for (size_t i=0; i<lfsu_s.size(); i++)
              {
                gradu_s.axpy(x_s[i],gradphi_s[i]);
              }
            Dune::FieldVector<RF,dim> gradu_n(0.0);
            for (size_t i=0; i<lfsu_n.size(); i++)
              {
                gradu_n.axpy(x_n[i],gradphi_n[i]);
              }

            // compute K * gradient of u
            Dune::FieldVector<RF,dim> kgradu_s(0.0);
            permeability_s.umv(gradu_s,kgradu_s);
            Dune::FieldVector<RF,dim> kgradu_n(0.0);
            permeability_n.umv(gradu_n,kgradu_n);

            // evaluate u in both elements self and neighbor
            RF u_s = 0.0;
            for (size_t i=0; i<lfsu_s.size(); i++)
              {
                u_s += x_s[i]*phi_s[i];
              }
            RF u_n = 0.0;
            for (size_t i=0; i<lfsu_n.size(); i++)
              {
                u_n += x_n[i]*phi_n[i];
              }

            // jump and average for u
            RF u_jump = u_s - u_n;

            // average on intersection of K * grad u * normal
            RF kgradunormal_average = (kgradu_s + kgradu_n)*normal * 0.5;

            // average on intersection of K * grad v * normal
            std::vector<Dune::FieldVector<RF,dim> > kgradphi_s(lfsu_s.size());
            std::vector<Dune::FieldVector<RF,dim> > kgradphi_n(lfsu_n.size());
            for (size_t i=0; i<lfsu_s.size(); i++)
              {
                permeability_s.mv(gradphi_s[i],kgradphi_s[i]);
              }
            for (size_t i=0; i<lfsu_n.size(); i++)
              {
                permeability_n.mv(gradphi_n[i],kgradphi_n[i]);
              }

            // integrate what needed
            RF factor = it->weight()*ig.geometry().integrationElement(it->position());
            for (unsigned int i=0; i<lfsv_s.size(); i++)
              {
                // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
                r_s[i] += penalty_weight_s * u_jump*phi_s[i]*factor;
                // epsilon * <Kgradv*my>[u] - <Kgradu*my>[v]
                r_s[i] += epsilon*(kgradphi_s[i]*normal)*0.5*u_jump*factor;
                r_s[i] -= phi_s[i]*kgradunormal_average*factor;
              }
            for (unsigned int i=0; i<lfsv_n.size(); i++)
              {
                // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
                r_n[i] += penalty_weight_n * u_jump*(-phi_n[i])*factor;
                // epsilon * <Kgradv*my>[u] - [v]<Kgradu*my>
                r_n[i] += epsilon*(kgradphi_n[i]*normal)*0.5*u_jump*factor;
                r_n[i] -= (-phi_n[i])*kgradunormal_average*factor;
              }
          }
      }

      // boundary integral depending on test and ansatz functions
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;

        // dimensions
        const int dim = IG::dimension;
        const int dimw = IG::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        // evaluate boundary condition type
        typename B::Traits::RangeType bctype;
        b.evaluate(ig,rule.begin()->position(),bctype);

        // Dirichlet boundary condition
        if (bctype > 0.0)
          {
            // center in face's reference element
            const Dune::FieldVector<DF,IG::dimension-1>& face_center =
              Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
            // outer normal, assuming it is constant over whole intersection
            const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(face_center);

            // evaluate diffusion tensor at cell center, assume it is constant over elements
            const Dune::FieldVector<DF,IG::dimension>
              localcenter = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
            typename K::Traits::RangeType tensor(0.0);
            k.evaluate(*ig.inside(),localcenter,tensor);

            // penalty weight for NIPG / SIPG
            RF penalty_weight = sigma / pow(ig.inside()->geometry().volume(), beta);

            // loop over quadrature points and integrate u * phi
            for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
              {
                // position of quadrature point in local coordinates of element
                Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

                // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
                std::vector<JacobianType> js(lfsv.size());
                lfsv.localFiniteElement().localBasis().evaluateJacobian(local,js);

                // transform gradient to real element
                const Dune::FieldMatrix<DF,dimw,dim> jac = ig.inside()->geometry().jacobianInverseTransposed(local);
                std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsv.size());
                for (size_t i=0; i<lfsv.size(); i++)
                  {
                    gradphi[i] = 0.0;
                    jac.umv(js[i][0],gradphi[i]);
                  }

                // evaluate test shape functions
                std::vector<RangeType> phi(lfsv.size());
                lfsv.localFiniteElement().localBasis().evaluateFunction(local,phi);

                // compute gradient of u
                Dune::FieldVector<RF,dim> gradu(0.0);
                for (size_t i=0; i<lfsu.size(); i++)
                  {
                    gradu.axpy(x[i],gradphi[i]);
                  }

                // compute K * gradient of u
                Dune::FieldVector<RF,dim> Kgradu(0.0);
                tensor.umv(gradu,Kgradu);

                // evaluate u
                RF u=0.0;
                for (size_t i=0; i<lfsu.size(); i++)
                  {
                    u += x[i]*phi[i];
                  }

                // integrate u
                RF factor = it->weight()*ig.geometry().integrationElement(it->position());
                for (size_t i=0; i<lfsv.size(); i++)
                  {
                    // Left hand side of NIPG / SIPG penalty term: sigma/|gamma|^beta*u*v
                    r[i] += penalty_weight*u*phi[i]*factor;
                    // epsilon*K*gradient of v*my*u - v*K*gradient of u*my
                    Dune::FieldVector<RF,dim> Kgradv(0.0);
                    tensor.umv(gradphi[i],Kgradv);
                    r[i] += epsilon * (Kgradv*normal)*u*factor;
                    r[i] -= phi[i]*(Kgradu*normal)*factor;
                  }
              }
          }
      }

      // volume integral depending only on test functions,
      // contains f on the right hand side
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate shape functions
            std::vector<RangeType> phi(lfsv.size());
            lfsv.localFiniteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate right hand side parameter function
            typename F::Traits::RangeType y;
            f.evaluate(eg.entity(),it->position(),y);

            // integrate f
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_t i=0; i<lfsv.size(); i++)
              {
                // f*v
                r[i] -= y*phi[i]*factor;
              }
          }
      }

      // boundary integral independent of ansatz functions,
      // Neumann and Dirichlet boundary conditions, DG penalty term's right hand side
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;

        // dimensions
        const int dim = IG::dimension;
        const int dimw = IG::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        // evaluate boundary condition type
        typename B::Traits::RangeType bctype;
        b.evaluate(ig,rule.begin()->position(),bctype);

        // Neumann boundary condition
        if (bctype == 0.0)
          {
            // loop over quadrature points and integrate normal flux
            for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
              {
                // position of quadrature point in local coordinates of element
                Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

                // evaluate test shape functions
                std::vector<RangeType> phi(lfsv.size());
                lfsv.localFiniteElement().localBasis().evaluateFunction(local,phi);

                // evaluate flux boundary condition
                typename J::Traits::RangeType y(0.0);
                j.evaluate(*(ig.inside()),local,y);

                // integrate J
                RF factor = it->weight()*ig.geometry().integrationElement(it->position());
                for (size_t i=0; i<lfsv.size(); i++)
                  {
                    // Neumann boundary condition: j*v
                    r[i] += y*phi[i]*factor;
                  }
              }
          }
        // Dirichlet boundary condition
        else if (bctype > 0.0)
          {
            // center in face's reference element
            const Dune::FieldVector<DF,IG::dimension-1>& face_center =
              Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
            // outer normal, assuming it is constant over whole intersection
            const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(face_center);
            // evaluate diffusion tensor at cell center, assume it is constant over elements
            const Dune::FieldVector<DF,IG::dimension>
              localcenter = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
            typename K::Traits::RangeType tensor(0.0);
            k.evaluate(*ig.inside(),localcenter,tensor);
            // penalty weight for NIPG / SIPG
            RF penalty_weight = sigma / pow(ig.inside()->geometry().volume(), beta);

            // loop over quadrature points and integrate g * phi
            for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
              {
                // position of quadrature point in local coordinates of element
                Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

                // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
                std::vector<JacobianType> js(lfsv.size());
                lfsv.localFiniteElement().localBasis().evaluateJacobian(local,js);

                // transform gradient to real element
                const Dune::FieldMatrix<DF,dimw,dim> jac = ig.inside()->geometry().jacobianInverseTransposed(local);
                std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsv.size());
                for (size_t i=0; i<lfsv.size(); i++)
                  {
                    gradphi[i] = 0.0;
                    jac.umv(js[i][0],gradphi[i]);
                  }

                // evaluate test shape functions
                std::vector<RangeType> phi(lfsv.size());
                lfsv.localFiniteElement().localBasis().evaluateFunction(local,phi);

                // evaluate Dirichlet boundary condition
                typename G::Traits::RangeType y;
                g.evaluate(*(ig.inside()),local,y);

                // integrate G
                RF factor = it->weight()*ig.geometry().integrationElement(it->position());
                for (size_t i=0; i<lfsv.size(); i++)
                  {
                    // Right hand side of NIPG / SIPG penalty term: -sigma / |gamma|^beta * g
                    r[i] -= penalty_weight*y*phi[i]*factor;
                    // Dirichlet boundary: -epsilon*K*gradient of v*my*g
                    Dune::FieldVector<RF,dim> Kgradv(0.0);
                    tensor.umv(gradphi[i],Kgradv);
                    r[i] -= epsilon * (Kgradv*normal)*y*factor;
                  }
              }
          }
        else
          {
            DUNE_THROW(Dune::Exception, "Unrecognized or unsupported boundary condition. (" << bctype << ")");
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            LocalMatrix<R>& mat) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        typename K::Traits::RangeType tensor;
        Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
        k.evaluate(eg.entity(),localcenter,tensor);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> js(lfsu.size());
            lfsu.localFiniteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradient to real element
            const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            for (typename LFSU::Traits::SizeType i=0; i<lfsu.size(); i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // compute K * gradient of shape functions
            std::vector<Dune::FieldVector<RF,dim> > Kgradphi(lfsu.size());
            for (typename LFSU::Traits::SizeType i=0; i<lfsu.size(); i++)
              {
                tensor.mv(gradphi[i],Kgradphi[i]);
              }

            // integrate (K grad phi_j)*grad phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (typename LFSU::Traits::SizeType j=0; j<lfsu.size(); j++)
              {
                for (typename LFSU::Traits::SizeType i=0; i<lfsu.size(); i++)
                  {
                    // K*gradient of phi_j*K*gradient of phi_i
                    mat(i,j) += ( Kgradphi[j]*gradphi[i])*factor;
                  }
              }
          }
      }

      // jacobian of skeleton term
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                              LocalMatrix<R>& mat_ss, LocalMatrix<R>& mat_sn,
                              LocalMatrix<R>& mat_ns, LocalMatrix<R>& mat_nn) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;

        // dimensions
        const int dim = IG::dimension;
        const int dimw = IG::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        // center in face's reference element
        const Dune::FieldVector<DF,IG::dimension-1>& face_center =
          Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(face_center);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        const Dune::FieldVector<DF,IG::dimension>& 
          inside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>& 
          outside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);
        typename K::Traits::RangeType permeability_s(0.0);
        typename K::Traits::RangeType permeability_n(0.0);
        k.evaluate(*(ig.inside()),inside_local,permeability_s);
        k.evaluate(*(ig.outside()),outside_local,permeability_n);
        /*for (unsigned int i = 0; i < K::Traits::GridViewType::dimension; ++i)
          {
          for (unsigned int j = 0; j < K::Traits::GridViewType::dimension; ++j)
          {
          if (permeability_s[i][j] * permeability_n[i][j] == 0.0)
          {
          permeability[i][j] = 0.0;
          }
          else
          {
          permeability[i][j] = 2.0 / (1.0/permeability_s[i][j] + 1.0/permeability_n[i][j]);
          }
          }
          }*/

        // penalty weight for NIPG / SIPG
        RF penalty_weight_s = sigma / pow(ig.inside()->geometry().volume(), beta);
        RF penalty_weight_n = sigma / pow(ig.outside()->geometry().volume(), beta);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // position of quadrature point in local coordinates of element
            Dune::FieldVector<DF,dim> local_s = ig.geometryInInside().global(it->position());
            Dune::FieldVector<DF,dim> local_n = ig.geometryInOutside().global(it->position());

            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> js_s(lfsv_s.size());
            lfsv_s.localFiniteElement().localBasis().evaluateJacobian(local_s,js_s);
            std::vector<JacobianType> js_n(lfsv_n.size());
            lfsv_n.localFiniteElement().localBasis().evaluateJacobian(local_n,js_n);

            // transform gradient to real element
            const Dune::FieldMatrix<DF,dimw,dim> jac_s = ig.inside()->geometry().jacobianInverseTransposed(local_s);
            std::vector<Dune::FieldVector<RF,dim> > gradphi_s(lfsv_s.size());
            for (size_t i=0; i<lfsv_s.size(); i++)
              {
                gradphi_s[i] = 0.0;
                jac_s.umv(js_s[i][0],gradphi_s[i]);
              }
            const Dune::FieldMatrix<DF,dimw,dim> jac_n = ig.outside()->geometry().jacobianInverseTransposed(local_n);
            std::vector<Dune::FieldVector<RF,dim> > gradphi_n(lfsv_n.size());
            for (size_t i=0; i<lfsv_n.size(); i++)
              {
                gradphi_n[i] = 0.0;
                jac_n.umv(js_n[i][0],gradphi_n[i]);
              }

            // evaluate test shape functions
            std::vector<RangeType> phi_s(lfsv_s.size());
            lfsv_s.localFiniteElement().localBasis().evaluateFunction(local_s,phi_s);
            std::vector<RangeType> phi_n(lfsv_n.size());
            lfsv_n.localFiniteElement().localBasis().evaluateFunction(local_n,phi_n);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu_s(0.0);
            for (size_t i=0; i<lfsu_s.size(); i++)
              {
                gradu_s.axpy(x_s[i],gradphi_s[i]);
              }
            Dune::FieldVector<RF,dim> gradu_n(0.0);
            for (size_t i=0; i<lfsu_n.size(); i++)
              {
                gradu_n.axpy(x_n[i],gradphi_n[i]);
              }

            // compute K * gradient of u
            Dune::FieldVector<RF,dim> kgradu_s(0.0);
            permeability_s.umv(gradu_s,kgradu_s);
            Dune::FieldVector<RF,dim> kgradu_n(0.0);
            permeability_n.umv(gradu_n,kgradu_n);

            // evaluate u in both elements self and neighbor
            RF u_s = 0.0;
            for (size_t i=0; i<lfsu_s.size(); i++)
              {
                u_s += x_s[i]*phi_s[i];
              }
            RF u_n = 0.0;
            for (size_t i=0; i<lfsu_n.size(); i++)
              {
                u_n += x_n[i]*phi_n[i];
              }

            // average on intersection of K * grad v * normal
            std::vector<Dune::FieldVector<RF,dim> > kgradphi_s(lfsu_s.size());
            std::vector<Dune::FieldVector<RF,dim> > kgradphi_n(lfsu_n.size());
            for (size_t i=0; i<lfsu_s.size(); i++)
              {
                permeability_s.mv(gradphi_s[i],kgradphi_s[i]);
              }
            for (size_t i=0; i<lfsu_n.size(); i++)
              {
                permeability_n.mv(gradphi_n[i],kgradphi_n[i]);
              }

            // integrate what needed
            RF factor = it->weight()*ig.geometry().integrationElement(it->position());
            for (typename LFSU::Traits::SizeType j=0; j<lfsu_s.size(); j++)
              {
                for (typename LFSU::Traits::SizeType i=0; i<lfsu_s.size(); i++)
                  {
                    // epsilon*(K*gradient of phi_s_j + K*gradient of phi_n_j)/2*(phi_s_i - phi_n_i) - (phi_s_j - phi_n_j)*(K*gradient of phi_s_i + K*gradient of phi_n_i)/2
                    mat_ss(i,j) += (epsilon*0.5*(kgradphi_s[i]*normal)*(phi_s[j]) - (phi_s[i])*0.5*(kgradphi_s[j]*normal) )*factor;
                    // NIPG / SIPG term: (phi_n_j - phi_s_j)*(phi_n_i - phi_s_i)
                    mat_ss(i,j) += penalty_weight_s*(phi_s[j])*(phi_s[i]) *factor;
                  }
                for (typename LFSU::Traits::SizeType i=0; i<lfsu_n.size(); i++)
                  {
                    // epsilon*(K*gradient of phi_s_j + K*gradient of phi_n_j)/2*(phi_s_i - phi_n_i) - (phi_s_j - phi_n_j)*(K*gradient of phi_s_i + K*gradient of phi_n_i)/2
                    mat_ns(i,j) += (epsilon*0.5*(kgradphi_n[i]*normal)*(phi_s[j]) - (-phi_n[i])*0.5*(kgradphi_s[j]*normal) )*factor;
                    // NIPG / SIPG term: (phi_n_j - phi_s_j)*(phi_n_i - phi_s_i)
                    mat_ns(i,j) += penalty_weight_n*(phi_s[j])*(-phi_n[i]) *factor;
                  }
              }
            for (typename LFSU::Traits::SizeType j=0; j<lfsu_n.size(); j++)
              {
                for (typename LFSU::Traits::SizeType i=0; i<lfsu_s.size(); i++)
                  {
                    // epsilon*(K*gradient of phi_s_j + K*gradient of phi_n_j)/2*(phi_s_i - phi_n_i) - (phi_s_j - phi_n_j)*(K*gradient of phi_s_i + K*gradient of phi_n_i)/2
                    mat_sn(i,j) += (epsilon*0.5*(kgradphi_s[i]*normal)*(-phi_n[j]) - (phi_s[i])*0.5*(kgradphi_n[j]*normal) )*factor;
                    // NIPG / SIPG term: (phi_n_j - phi_s_j)*(phi_n_i - phi_s_i)
                    mat_sn(i,j) += penalty_weight_s*(-phi_n[j])*(phi_s[i]) *factor;
                  }
                for (typename LFSU::Traits::SizeType i=0; i<lfsu_n.size(); i++)
                  {
                    // epsilon*(K*gradient of phi_s_j + K*gradient of phi_n_j)/2*(phi_s_i - phi_n_i) - (phi_s_j - phi_n_j)*(K*gradient of phi_s_i + K*gradient of phi_n_i)/2
                    mat_nn(i,j) += (epsilon*0.5*(kgradphi_n[i]*normal)*(-phi_n[j]) - (-phi_n[i])*0.5*(kgradphi_n[j]*normal) )*factor;
                    // NIPG / SIPG term: (phi_n_j - phi_s_j)*(phi_n_i - phi_s_i)
                    mat_nn(i,j) += penalty_weight_n*(-phi_n[j])*(-phi_n[i]) *factor;
                  }
              }
          }
      }

      // jacobian of volume term
      // TODO: Ist die uebergebene Matrix initialisiert? Bastian nachfragen. ////////////////////////////////////////////////////////////////
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu, const X& x, const LFSV& lfsv,
                              LocalMatrix<R>& mat) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::LocalFiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;

        // dimensions
        const int dim = IG::dimension;
        const int dimw = IG::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        // evaluate boundary condition type
        typename B::Traits::RangeType bctype;
        b.evaluate(ig,rule.begin()->position(),bctype);

        // Dirichlet boundary condition
        if (bctype > 0.0)
          {
            // center in face's reference element
            const Dune::FieldVector<DF,IG::dimension-1>& face_center =
              Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
            // outer normal, assuming it is constant over whole intersection
            const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(face_center);

            // evaluate diffusion tensor at cell center, assume it is constant over elements
            const Dune::FieldVector<DF,IG::dimension>
              localcenter = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
            typename K::Traits::RangeType tensor(0.0);
            k.evaluate(*ig.inside(),localcenter,tensor);

            // penalty weight for NIPG / SIPG
            RF penalty_weight = sigma / pow(ig.inside()->geometry().volume(), beta);

            // loop over quadrature points and integrate u * phi
            for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
              {
                // position of quadrature point in local coordinates of element
                Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

                // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
                std::vector<JacobianType> js(lfsv.size());
                lfsv.localFiniteElement().localBasis().evaluateJacobian(local,js);

                // transform gradient to real element
                const Dune::FieldMatrix<DF,dimw,dim> jac = ig.inside()->geometry().jacobianInverseTransposed(local);
                std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsv.size());
                for (size_t i=0; i<lfsv.size(); i++)
                  {
                    gradphi[i] = 0.0;
                    jac.umv(js[i][0],gradphi[i]);
                  }

                // evaluate test shape functions
                std::vector<RangeType> phi(lfsv.size());
                lfsv.localFiniteElement().localBasis().evaluateFunction(local,phi);

                // compute gradient of u
                Dune::FieldVector<RF,dim> gradu(0.0);
                for (size_t i=0; i<lfsu.size(); i++)
                  {
                    gradu.axpy(x[i],gradphi[i]);
                  }

                // compute K * gradient of v
                std::vector<Dune::FieldVector<RF,dim> > kgradphi(lfsu.size());
                for (size_t i=0; i<lfsu.size(); i++)
                  {
                    tensor.mv(gradphi[i],kgradphi[i]);
                  }

                // integrate
                RF factor = it->weight()*ig.geometry().integrationElement(it->position());
                for (typename LFSU::Traits::SizeType j=0; j<lfsu.size(); j++)
                  {
                    for (typename LFSU::Traits::SizeType i=0; i<lfsu.size(); i++)
                      {
                        // epsilon*K*gradient of phi_i*my*phi_j - phi_i*K*gradient of phi_j*my
                        mat(i,j) += (epsilon*(kgradphi[i]*normal)*phi[j] - phi[i]*(kgradphi[j]*normal))*factor;
                        // NIPG / SIPG penalty term: sigma / |gamma|^beta *phi_j*phi_i
                        mat(i,j) += (penalty_weight*phi[j]*phi[i])*factor;
                      }
                  }
              }
          }
      }

    private:
      const K& k;
      const F& f;
      const B& b;
      const G& g;
      const J& j;
      // select dg method
      double dg_method;
      // values for NIPG / NIPG
      double epsilon;
      double sigma;
      double beta;
      int qorder;
    };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
