// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CONVECTIONDIFFUSIONDG_HH
#define DUNE_PDELAB_CONVECTIONDIFFUSIONDG_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

#include"convectiondiffusionparameter.hh"

#define USECACHE 1

namespace Dune {
  namespace PDELab {

    struct ConvectionDiffusionDGMethod
    {
      enum Type { NIPG, SIPG };
    };

    struct ConvectionDiffusionDGWeights
    {
      enum Type { weightsOn, weightsOff };
    };

    /** a local operator for solving the convection-diffusion equation with discontinuous Galerkin
     *  
     * \f{align*}{
     *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \\
     *                                              u &=& g \mbox{ on } \partial\Omega_D \\
     *                (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
     *                        -(A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_O
     * \f}
     * Note:
     *  - This formulation is valid for velocity fields which are non-divergence free.
     *  - Outflow boundary conditions should only be set on the outflow boundary
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T, typename FiniteElementMap>
    class ConvectionDiffusionDG 
      : public Dune::PDELab::NumericalJacobianApplyVolume<ConvectionDiffusionDG<T,FiniteElementMap> >,
        public Dune::PDELab::NumericalJacobianApplySkeleton<ConvectionDiffusionDG<T,FiniteElementMap> >,
        public Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionDG<T,FiniteElementMap> >,
        public Dune::PDELab::FullSkeletonPattern, 
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::LocalOperatorDefaultFlags,
        public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>
    {
      enum { dim = T::Traits::GridViewType::dimension };
 
      typedef typename T::Traits::RangeFieldType Real;
      typedef typename ConvectionDiffusionBoundaryConditions::Type BCType;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume  = true };

      //! constructor: pass parameter object
      ConvectionDiffusionDG (T& param_, 
                             ConvectionDiffusionDGMethod::Type method_=ConvectionDiffusionDGMethod::NIPG, 
                             ConvectionDiffusionDGWeights::Type weights_=ConvectionDiffusionDGWeights::weightsOff,
                             Real alpha_=0.0, int intorderadd_=0) 
        : param(param_), method(method_), weights(weights_),
          alpha(alpha_), intorderadd(intorderadd_), quadrature_factor(2),
          Dune::PDELab::NumericalJacobianApplyVolume<ConvectionDiffusionDG<T,FiniteElementMap> >(1.0e-7),
          Dune::PDELab::NumericalJacobianApplySkeleton<ConvectionDiffusionDG<T,FiniteElementMap> >(1.0e-7),
          Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionDG<T,FiniteElementMap> >(1.0e-7)
      {
        theta = 1.0;
        if (method==ConvectionDiffusionDGMethod::SIPG) theta = -1.0;
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;
        
        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;
        const int intorder = intorderadd+quadrature_factor*lfsu.finiteElement().localBasis().order();

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        typename T::Traits::PermTensorType A;
        Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
        A = param.A(eg.entity(),localcenter);

        // transformation
        Dune::FieldMatrix<DF,dimw,dim> jac;

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
#if USECACHE==0
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);
#else
            const std::vector<RangeType>& phi = cache.evaluateFunction(it->position(),lfsu.finiteElement().localBasis());
#endif

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++) u += x[lfsu.localIndex(i)]*phi[i];

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
#if USECACHE==0
            std::vector<JacobianType> js(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);
#else
            const std::vector<JacobianType>& js = cache.evaluateJacobian(it->position(),lfsu.finiteElement().localBasis());
#endif

            // transform gradients of shape functions to real element
            jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            for (size_type i=0; i<lfsu.size(); i++)
              jac.mv(js[i][0],gradphi[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (size_type i=0; i<lfsu.size(); i++)
              gradu.axpy(x[lfsu.localIndex(i)],gradphi[i]);

            // compute K * gradient of u
            Dune::FieldVector<RF,dim> Agradu(0.0);
            A.umv(gradu,Agradu);

            // evaluate velocity field
            typename T::Traits::RangeType b = param.b(eg.entity(),it->position());

            // evaluate reaction term
            typename T::Traits::RangeFieldType c = param.c(eg.entity(),it->position());

            // integrate (K grad u - bu)*grad phi_i + a*u*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsv.size(); i++)
              r[lfsv.localIndex(i)] += ( Agradu*gradphi[i] - u*(b*gradphi[i]) + c*u*phi[i] )*factor;
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
                            Dune::PDELab::LocalMatrix<R>& mat) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;
        
        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;
        const int intorder = intorderadd+quadrature_factor*lfsu.finiteElement().localBasis().order();

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        typename T::Traits::PermTensorType A;
        Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
        A = param.A(eg.entity(),localcenter);

        // transformation
        Dune::FieldMatrix<DF,dimw,dim> jac;

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
#if USECACHE==0
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);
#else
            const std::vector<RangeType>& phi = cache.evaluateFunction(it->position(),lfsu.finiteElement().localBasis());
#endif

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
#if USECACHE==0
            std::vector<JacobianType> js(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);
#else
            const std::vector<JacobianType>& js = cache.evaluateJacobian(it->position(),lfsu.finiteElement().localBasis());
#endif

            // transform gradients of shape functions to real element
            jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            std::vector<Dune::FieldVector<RF,dim> > Agradphi(lfsu.size());
            for (size_type i=0; i<lfsu.size(); i++)
              {
                jac.mv(js[i][0],gradphi[i]);
                A.mv(gradphi[i],Agradphi[i]);
              }

            // evaluate velocity field
            typename T::Traits::RangeType b = param.b(eg.entity(),it->position());

            // evaluate reaction term
            typename T::Traits::RangeFieldType c = param.c(eg.entity(),it->position());

            // integrate (K grad u - bu)*grad phi_i + a*u*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat(lfsu.localIndex(i),lfsu.localIndex(j)) +=  ( Agradphi[j]*gradphi[i] - phi[j]*(b*gradphi[i]) + c*phi[j]*phi[i] )*factor;
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
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        const int intorder = intorderadd+quadrature_factor*std::max(lfsu_s.finiteElement().localBasis().order(),
                                                                    lfsu_n.finiteElement().localBasis().order());
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>& 
          outside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);
        typename T::Traits::PermTensorType A_s, A_n;
        A_s = param.A(*(ig.inside()),inside_local);
        A_n = param.A(*(ig.outside()),outside_local);

        // face diameter; this should be revised for anisotropic meshes?
        DF h_s, h_n;
        DF hmax_s, hmax_n;
        element_size(ig.inside()->geometry(),h_s,hmax_s);
        element_size(ig.outside()->geometry(),h_n,hmax_n);
        RF h_F = std::min(h_s,h_n);
        h_F = std::min(ig.inside()->geometry().volume(),ig.outside()->geometry().volume())/ig.geometry().volume(); // Houston!

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

        // tensor times normal
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        Dune::FieldVector<RF,dim> An_F_n;
        A_n.mv(n_F,An_F_n);

        // compute weights
        RF omega_s;
        RF omega_n;
        RF harmonic_average(0.0);
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          {
            RF delta_s = (An_F_s*n_F);
            RF delta_n = (An_F_n*n_F);
            omega_s = delta_n/(delta_s+delta_n+1e-20);
            omega_n = delta_s/(delta_s+delta_n+1e-20);
            harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
          }
        else
          {
            omega_s = omega_n = 0.5;
            harmonic_average = 1.0;
          }

        // get polynomial degree
        int degree = std::max(lfsu_s.finiteElement().localBasis().order(),
                              lfsu_n.finiteElement().localBasis().order());

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // exact normal
            const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(it->position());

            // position of quadrature point in local coordinates of elements 
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());
            Dune::FieldVector<DF,dim> iplocal_n = ig.geometryInOutside().global(it->position());

            // evaluate basis functions
#if USECACHE==0
            std::vector<RangeType> phi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);
            std::vector<RangeType> phi_n(lfsu_n.size());
            lfsu_n.finiteElement().localBasis().evaluateFunction(iplocal_n,phi_n);
#else
            const std::vector<RangeType>& phi_s = cache.evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());
            const std::vector<RangeType>& phi_n = cache.evaluateFunction(iplocal_n,lfsu_n.finiteElement().localBasis());
#endif

            // evaluate u
            RF u_s=0.0; for (size_type i=0; i<lfsu_s.size(); i++) u_s += x_s[lfsu_s.localIndex(i)]*phi_s[i];
            RF u_n=0.0; for (size_type i=0; i<lfsu_n.size(); i++) u_n += x_n[lfsu_n.localIndex(i)]*phi_n[i];

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
#if USECACHE==0
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
            std::vector<JacobianType> gradphi_n(lfsu_n.size());
            lfsu_n.finiteElement().localBasis().evaluateJacobian(iplocal_n,gradphi_n);
#else
            const std::vector<JacobianType>& gradphi_s = cache.evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
            const std::vector<JacobianType>& gradphi_n = cache.evaluateJacobian(iplocal_n,lfsu_n.finiteElement().localBasis());
#endif

            // transform gradients of shape functions to real element
            jac = ig.inside()->geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
            jac = ig.outside()->geometry().jacobianInverseTransposed(iplocal_n);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_n(lfsu_n.size());
            for (size_type i=0; i<lfsu_n.size(); i++) jac.mv(gradphi_n[i][0],tgradphi_n[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu_s(0.0);
            for (size_type i=0; i<lfsu_s.size(); i++) gradu_s.axpy(x_s[lfsu_s.localIndex(i)],tgradphi_s[i]);
            Dune::FieldVector<RF,dim> gradu_n(0.0);
            for (size_type i=0; i<lfsu_n.size(); i++) gradu_n.axpy(x_n[lfsu_n.localIndex(i)],tgradphi_n[i]);

            // evaluate velocity field and upwinding, assume H(div) velocity field => may choose any side
            typename T::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
            RF normalflux = b*n_F_local;
            RF omegaup_s, omegaup_n;
            if (normalflux>=0.0) 
              {
                omegaup_s = 1.0;
                omegaup_n = 0.0;
              }
            else
              {
                omegaup_s = 0.0;
                omegaup_n = 1.0;
              }

            // integration factor
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());

            // convection term
            RF term1 = (omegaup_s*u_s + omegaup_n*u_n) * normalflux *factor;
            for (size_type i=0; i<lfsu_s.size(); i++) r_s[lfsu_s.localIndex(i)] += term1 * phi_s[i];
            for (size_type i=0; i<lfsu_n.size(); i++) r_n[lfsu_n.localIndex(i)] -= term1 * phi_n[i];

            // diffusion term
            RF term2 =  -(omega_s*(An_F_s*gradu_s) + omega_n*(An_F_n*gradu_n)) * factor;
            for (size_type i=0; i<lfsu_s.size(); i++) r_s[lfsu_s.localIndex(i)] += term2 * phi_s[i];
            for (size_type i=0; i<lfsu_n.size(); i++) r_n[lfsu_n.localIndex(i)] -= term2 * phi_n[i];

            // (non-)symmetric IP term
            RF term3 = (u_s-u_n) * factor;
            for (size_type i=0; i<lfsu_s.size(); i++) 
              r_s[lfsu_s.localIndex(i)] += term3 * theta * omega_s * (An_F_s*tgradphi_s[i]);
            for (size_type i=0; i<lfsu_n.size(); i++) 
              r_n[lfsu_n.localIndex(i)] += term3 * theta * omega_n * (An_F_n*tgradphi_n[i]);

            // standard IP term integral
            RF term4 = penalty_factor * (u_s-u_n) * factor;
            for (size_type i=0; i<lfsu_s.size(); i++) r_s[lfsu_s.localIndex(i)] += term4 * phi_s[i];
            for (size_type i=0; i<lfsu_n.size(); i++) r_n[lfsu_n.localIndex(i)] -= term4 * phi_n[i];
          }
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void jacobian_skeleton (const IG& ig, 
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                              Dune::PDELab::LocalMatrix<R>& mat_ss, Dune::PDELab::LocalMatrix<R>& mat_sn, 
                              Dune::PDELab::LocalMatrix<R>& mat_ns, Dune::PDELab::LocalMatrix<R>& mat_nn) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        const int intorder = intorderadd+quadrature_factor*std::max(lfsu_s.finiteElement().localBasis().order(),
                                                                    lfsu_n.finiteElement().localBasis().order());
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>& 
          outside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);
        typename T::Traits::PermTensorType A_s, A_n;
        A_s = param.A(*(ig.inside()),inside_local);
        A_n = param.A(*(ig.outside()),outside_local);

        // face diameter; this should be revised for anisotropic meshes?
        DF h_s, h_n;
        DF hmax_s, hmax_n;
        element_size(ig.inside()->geometry(),h_s,hmax_s);
        element_size(ig.outside()->geometry(),h_n,hmax_n);
        RF h_F = std::min(h_s,h_n);
        h_F = std::min(ig.inside()->geometry().volume(),ig.outside()->geometry().volume())/ig.geometry().volume(); // Houston!

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

        // tensor times normal
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        Dune::FieldVector<RF,dim> An_F_n;
        A_n.mv(n_F,An_F_n);

        // compute weights
        RF omega_s;
        RF omega_n;
        RF harmonic_average(0.0);
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          {
            RF delta_s = (An_F_s*n_F);
            RF delta_n = (An_F_n*n_F);
            omega_s = delta_n/(delta_s+delta_n+1e-20);
            omega_n = delta_s/(delta_s+delta_n+1e-20);
            harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
          }
        else
          {
            omega_s = omega_n = 0.5;
            harmonic_average = 1.0;
          }

        // get polynomial degree
        int degree = std::max(lfsu_s.finiteElement().localBasis().order(),
                              lfsu_n.finiteElement().localBasis().order());

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // exact normal
            const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(it->position());

            // position of quadrature point in local coordinates of elements 
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());
            Dune::FieldVector<DF,dim> iplocal_n = ig.geometryInOutside().global(it->position());

            // evaluate basis functions
#if USECACHE==0
            std::vector<RangeType> phi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);
            std::vector<RangeType> phi_n(lfsu_n.size());
            lfsu_n.finiteElement().localBasis().evaluateFunction(iplocal_n,phi_n);
#else
            const std::vector<RangeType>& phi_s = cache.evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());
            const std::vector<RangeType>& phi_n = cache.evaluateFunction(iplocal_n,lfsu_n.finiteElement().localBasis());
#endif

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
#if USECACHE==0
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
            std::vector<JacobianType> gradphi_n(lfsu_n.size());
            lfsu_n.finiteElement().localBasis().evaluateJacobian(iplocal_n,gradphi_n);
#else
            const std::vector<JacobianType>& gradphi_s = cache.evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
            const std::vector<JacobianType>& gradphi_n = cache.evaluateJacobian(iplocal_n,lfsu_n.finiteElement().localBasis());
#endif

            // transform gradients of shape functions to real element
            jac = ig.inside()->geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
            jac = ig.outside()->geometry().jacobianInverseTransposed(iplocal_n);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_n(lfsu_n.size());
            for (size_type i=0; i<lfsu_n.size(); i++) jac.mv(gradphi_n[i][0],tgradphi_n[i]);

            // evaluate velocity field and upwinding, assume H(div) velocity field => may choose any side
            typename T::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
            RF normalflux = b*n_F_local;
            RF omegaup_s, omegaup_n;
            if (normalflux>=0.0) 
              {
                omegaup_s = 1.0;
                omegaup_n = 0.0;
              }
            else
              {
                omegaup_s = 0.0;
                omegaup_n = 1.0;
              }

            // integration factor
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());
            RF ipfactor = penalty_factor * factor;

            // do all terms in the order: I convection, II diffusion, III consistency, IV ip
            for (size_type j=0; j<lfsu_s.size(); j++) {
              RF temp1 = -(An_F_s*tgradphi_s[j])*omega_s*factor;
              for (size_type i=0; i<lfsu_s.size(); i++) {
                mat_ss(lfsu_s.localIndex(i),lfsu_s.localIndex(j)) += omegaup_s * phi_s[j] * normalflux *factor * phi_s[i];
                mat_ss(lfsu_s.localIndex(i),lfsu_s.localIndex(j)) += temp1 * phi_s[i];
                mat_ss(lfsu_s.localIndex(i),lfsu_s.localIndex(j)) += phi_s[j] * factor * theta * omega_s * (An_F_s*tgradphi_s[i]);
                mat_ss(lfsu_s.localIndex(i),lfsu_s.localIndex(j)) += phi_s[j] * ipfactor * phi_s[i];
              }
            }
            for (size_type j=0; j<lfsu_n.size(); j++) {
              RF temp1 = -(An_F_n*tgradphi_n[j])*omega_n*factor;
              for (size_type i=0; i<lfsu_s.size(); i++) {
                mat_sn(lfsu_s.localIndex(i),lfsu_n.localIndex(j)) += omegaup_n * phi_n[j] * normalflux *factor * phi_s[i];
                mat_sn(lfsu_s.localIndex(i),lfsu_n.localIndex(j)) += temp1 * phi_s[i];
                mat_sn(lfsu_s.localIndex(i),lfsu_n.localIndex(j)) += -phi_n[j] * factor * theta * omega_s * (An_F_s*tgradphi_s[i]);
                mat_sn(lfsu_s.localIndex(i),lfsu_n.localIndex(j)) += -phi_n[j] * ipfactor * phi_s[i];
              }
            }
            for (size_type j=0; j<lfsu_s.size(); j++) {
              RF temp1 = -(An_F_s*tgradphi_s[j])*omega_s*factor;
              for (size_type i=0; i<lfsu_n.size(); i++) {
                mat_ns(lfsu_n.localIndex(i),lfsu_s.localIndex(j)) -= omegaup_s * phi_s[j] * normalflux *factor * phi_n[i];
                mat_ns(lfsu_n.localIndex(i),lfsu_s.localIndex(j)) -= temp1 * phi_n[i];
                mat_ns(lfsu_n.localIndex(i),lfsu_s.localIndex(j)) += phi_s[j] * factor * theta * omega_n * (An_F_n*tgradphi_n[i]);
                mat_ns(lfsu_n.localIndex(i),lfsu_s.localIndex(j)) -= phi_s[j] * ipfactor * phi_n[i];
              }
            }
            for (size_type j=0; j<lfsu_n.size(); j++) {
              RF temp1 = -(An_F_n*tgradphi_n[j])*omega_n*factor;
              for (size_type i=0; i<lfsu_n.size(); i++) {
                mat_nn(lfsu_n.localIndex(i),lfsu_n.localIndex(j)) -= omegaup_n * phi_n[j] * normalflux *factor * phi_n[i];
                mat_nn(lfsu_n.localIndex(i),lfsu_n.localIndex(j)) -= temp1 * phi_n[i];
                mat_nn(lfsu_n.localIndex(i),lfsu_n.localIndex(j)) += -phi_n[j] * factor * theta * omega_n * (An_F_n*tgradphi_n[i]);
                mat_nn(lfsu_n.localIndex(i),lfsu_n.localIndex(j)) -= -phi_n[j] * ipfactor * phi_n[i];
              }
            }
          }
      }

      // boundary integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        const int intorder = intorderadd+quadrature_factor*lfsu_s.finiteElement().localBasis().order();
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        typename T::Traits::PermTensorType A_s;
        A_s = param.A(*(ig.inside()),inside_local);

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // face diameter
        DF h_s;
        DF hmax_s;
        element_size(ig.inside()->geometry(),h_s,hmax_s);
        RF h_F = h_s;
        h_F = ig.inside()->geometry().volume()/ig.geometry().volume(); // Houston!

        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

        // evaluate boundary condition
        const Dune::FieldVector<DF,dim-1> 
          face_local = Dune::GenericReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        BCType bctype = param.bctype(ig.intersection(),face_local);

        // compute weights
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        RF harmonic_average;
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          harmonic_average = An_F_s*n_F;
        else
          harmonic_average = 1.0;

        // get polynomial degree
        int degree = lfsu_s.finiteElement().localBasis().order();

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // position of quadrature point in local coordinates of elements 
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());

            // local normal
            const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(it->position());

            // evaluate basis functions
#if USECACHE==0
            std::vector<RangeType> phi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);
#else
            const std::vector<RangeType>& phi_s = cache.evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());
#endif

            // integration factor
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());

            if (bctype == ConvectionDiffusionBoundaryConditions::Neumann)
              {
                // evaluate flux boundary condition
                RF j = param.j(ig.intersection(),it->position());
                
                // integrate
                for (size_type i=0; i<lfsv_s.size(); i++) r_s[lfsu_s.localIndex(i)] += j * phi_s[i] * factor;

                continue;
              }

            // evaluate u
            RF u_s=0.0; for (size_type i=0; i<lfsu_s.size(); i++) u_s += x_s[lfsu_s.localIndex(i)]*phi_s[i];

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            typename T::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
            RF normalflux = b*n_F_local;

            if (bctype == ConvectionDiffusionBoundaryConditions::Outflow)
              {
                if (normalflux<-1e-30)
                  DUNE_THROW(Dune::Exception,"Outflow boundary condition on inflow!");

                // convection term
                RF term1 = u_s * normalflux *factor;
                for (size_type i=0; i<lfsu_s.size(); i++) r_s[lfsu_s.localIndex(i)] += term1 * phi_s[i];

                // evaluate flux boundary condition
                RF o = param.o(ig.intersection(),it->position());

                // integrate
                for (size_type i=0; i<lfsv_s.size(); i++) r_s[lfsu_s.localIndex(i)] += o * phi_s[i] * factor;

                continue;
              }

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
#if USECACHE==0
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
#else
            const std::vector<JacobianType>& gradphi_s = cache.evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
#endif

            // transform gradients of shape functions to real element
            jac = ig.inside()->geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu_s(0.0);
            for (size_type i=0; i<lfsu_s.size(); i++) gradu_s.axpy(x_s[lfsu_s.localIndex(i)],tgradphi_s[i]);

            // evaluate Dirichlet boundary condition
            RF g = param.g(*(ig.inside()),iplocal_s);

            // upwind
            RF omegaup_s, omegaup_n;
            if (normalflux>=0.0) 
              {
                omegaup_s = 1.0;
                omegaup_n = 0.0;
              }
            else
              {
                omegaup_s = 0.0;
                omegaup_n = 1.0;
              }

            // convection term
            RF term1 = (omegaup_s*u_s + omegaup_n*g) * normalflux *factor;
            for (size_type i=0; i<lfsu_s.size(); i++) r_s[lfsu_s.localIndex(i)] += term1 * phi_s[i];

            // diffusion term
            RF term2 =  (An_F_s*gradu_s) * factor;
            for (size_type i=0; i<lfsu_s.size(); i++) r_s[lfsu_s.localIndex(i)] -= term2 * phi_s[i];

            // (non-)symmetric IP term
            RF term3 = (u_s-g) * factor;
            for (size_type i=0; i<lfsu_s.size(); i++) r_s[lfsu_s.localIndex(i)] += term3 * theta * (An_F_s*tgradphi_s[i]);

            // standard IP term
            RF term4 = penalty_factor * (u_s-g) * factor;
            for (size_type i=0; i<lfsu_s.size(); i++) r_s[lfsu_s.localIndex(i)] += term4 * phi_s[i];
          }
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              Dune::PDELab::LocalMatrix<R>& mat_ss) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        const int intorder = intorderadd+quadrature_factor*lfsu_s.finiteElement().localBasis().order();
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        typename T::Traits::PermTensorType A_s;
        A_s = param.A(*(ig.inside()),inside_local);

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // face diameter
        DF h_s;
        DF hmax_s;
        element_size(ig.inside()->geometry(),h_s,hmax_s);
        RF h_F = h_s;
        h_F = ig.inside()->geometry().volume()/ig.geometry().volume(); // Houston!

        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

        // evaluate boundary condition
        const Dune::FieldVector<DF,dim-1> 
          face_local = Dune::GenericReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        BCType bctype = param.bctype(ig.intersection(),face_local);

        // compute weights
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        RF harmonic_average;
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          harmonic_average = An_F_s*n_F;
        else
          harmonic_average = 1.0;

        // get polynomial degree
        int degree = lfsu_s.finiteElement().localBasis().order();

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // Neumann boundary makes no contribution to boundary
        if (bctype == ConvectionDiffusionBoundaryConditions::Neumann) return;

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // position of quadrature point in local coordinates of elements 
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());

            // local normal
            const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(it->position());

            // evaluate basis functions
#if USECACHE==0
            std::vector<RangeType> phi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);
#else
            const std::vector<RangeType>& phi_s = cache.evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());
#endif

            // integration factor
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            typename T::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
            RF normalflux = b*n_F_local;

            if (bctype == ConvectionDiffusionBoundaryConditions::Outflow)
              {
                if (normalflux<-1e-30)
                  DUNE_THROW(Dune::Exception,"Outflow boundary condition on inflow!");

                // convection term
                for (size_type j=0; j<lfsu_s.size(); j++) 
                  for (size_type i=0; i<lfsu_s.size(); i++) 
                    mat_ss(lfsu_s.localIndex(i),lfsu_s.localIndex(j)) += phi_s[j] * normalflux * factor * phi_s[i];

                continue;
              }

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
#if USECACHE==0
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
#else
            const std::vector<JacobianType>& gradphi_s = cache.evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
#endif

            // transform gradients of shape functions to real element
            jac = ig.inside()->geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);

            // upwind
            RF omegaup_s, omegaup_n;
            if (normalflux>=0.0) 
              {
                omegaup_s = 1.0;
                omegaup_n = 0.0;
              }
            else
              {
                omegaup_s = 0.0;
                omegaup_n = 1.0;
              }

            // convection term
            for (size_type j=0; j<lfsu_s.size(); j++) 
              for (size_type i=0; i<lfsu_s.size(); i++) 
                mat_ss(lfsu_s.localIndex(i),lfsu_s.localIndex(j)) += omegaup_s * phi_s[j] * normalflux * factor * phi_s[i];

            // diffusion term
            for (size_type j=0; j<lfsu_s.size(); j++) 
              for (size_type i=0; i<lfsu_s.size(); i++) 
                mat_ss(lfsu_s.localIndex(i),lfsu_s.localIndex(j)) -= (An_F_s*tgradphi_s[j]) * factor * phi_s[i];

            // (non-)symmetric IP term
            for (size_type j=0; j<lfsu_s.size(); j++) 
              for (size_type i=0; i<lfsu_s.size(); i++) 
                mat_ss(lfsu_s.localIndex(i),lfsu_s.localIndex(j)) += phi_s[j] * factor * theta * (An_F_s*tgradphi_s[i]);

            // standard IP term
            for (size_type j=0; j<lfsu_s.size(); j++) 
              for (size_type i=0; i<lfsu_s.size(); i++) 
                mat_ss(lfsu_s.localIndex(i),lfsu_s.localIndex(j)) +=  penalty_factor * phi_s[j] * phi_s[i] * factor;
          }
      }

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = EG::Geometry::dimension;
        const int intorder = intorderadd+2*lfsv.finiteElement().localBasis().order();
        
        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate shape functions 
#if USECACHE==0
            std::vector<RangeType> phi(lfsv.size());
            lfsv.finiteElement().localBasis().evaluateFunction(it->position(),phi);
#else
            const std::vector<RangeType>& phi = cache.evaluateFunction(it->position(),lfsv.finiteElement().localBasis());
#endif

            // evaluate right hand side parameter function
            Real f;
            f = param.f(eg.entity(),it->position());

            // integrate f
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsv.size(); i++)
              r[lfsv.localIndex(i)] -= f*phi[i]*factor;
          }
      }

    private:
      T& param;  // two phase parameter class
      ConvectionDiffusionDGMethod::Type method;
      ConvectionDiffusionDGWeights::Type weights;
      Real alpha, beta;
      int quadrature_factor;
      int intorderadd;
      Real theta;
      typedef typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;
      Dune::PDELab::LocalBasisCache<LocalBasisType> cache;

      template<class GEO>
      void element_size (const GEO& geo, typename GEO::ctype& hmin, typename GEO::ctype hmax) const
      {
        typedef typename GEO::ctype DF;
        hmin = 1.0E100;
        hmax = -1.0E00;
        const int dim = GEO::coorddimension;
        if (dim==1)
          {
            Dune::FieldVector<DF,dim> x = geo.corner(0);
            x -= geo.corner(1);
            hmin = hmax = x.two_norm();
            return;
          }
        else
          {
            Dune::GeometryType gt = geo.type();
            for (int i=0; i<Dune::GenericReferenceElements<DF,dim>::general(gt).size(dim-1); i++)
              {
                Dune::FieldVector<DF,dim> x = geo.corner(Dune::GenericReferenceElements<DF,dim>::general(gt).subEntity(i,dim-1,0,dim));
                x -= geo.corner(Dune::GenericReferenceElements<DF,dim>::general(gt).subEntity(i,dim-1,1,dim));
                hmin = std::min(hmin,x.two_norm());
                hmax = std::max(hmax,x.two_norm());
              }
            return;
          }
      }
    };
  }
}
#endif
