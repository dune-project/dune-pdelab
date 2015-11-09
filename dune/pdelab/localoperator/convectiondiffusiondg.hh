// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CONVECTIONDIFFUSIONDG_HH
#define DUNE_PDELAB_CONVECTIONDIFFUSIONDG_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

#include"convectiondiffusionparameter.hh"

#ifndef USECACHE
#define USECACHE 1
//#define USECACHE 0
#endif

/**
   \todo update quadrature order to work with lfsv != lfsu
   \todo update alpha_* to work with lfsv != lfsu (./)
   \todo update jacobian_* to work with lfsv != lfsu
   \todo update caches to work with lfsv != lfsu
 */

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
     *                (b(x) u - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
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
                             Real alpha_=0.0,
                             int intorderadd_=0
                             )
        : Dune::PDELab::NumericalJacobianApplyVolume<ConvectionDiffusionDG<T,FiniteElementMap> >(1.0e-7),
          Dune::PDELab::NumericalJacobianApplySkeleton<ConvectionDiffusionDG<T,FiniteElementMap> >(1.0e-7),
          Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionDG<T,FiniteElementMap> >(1.0e-7),
          param(param_), method(method_), weights(weights_),
          alpha(alpha_), intorderadd(intorderadd_), quadrature_factor(2),
          cache(20)
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
        const int dim = EG::Entity::dimension;
        const int order = std::max(lfsu.finiteElement().localBasis().order(),
            lfsv.finiteElement().localBasis().order());
        const int intorder = intorderadd + quadrature_factor * order;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        typename T::Traits::PermTensorType A;
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
        A = param.A(eg.entity(),localcenter);

        // transformation
        typename EG::Geometry::JacobianInverseTransposed jac;

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            // evaluate basis functions
#if USECACHE==0
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(ip.position(),phi);
            std::vector<RangeType> psi(lfsv.size());
            lfsv.finiteElement().localBasis().evaluateFunction(ip.position(),psi);
#else
            const std::vector<RangeType>& phi = cache[order].evaluateFunction(ip.position(),lfsu.finiteElement().localBasis());
            const std::vector<RangeType>& psi = cache[order].evaluateFunction(ip.position(),lfsv.finiteElement().localBasis());
#endif

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // evaluate gradient of basis functions
#if USECACHE==0
            std::vector<JacobianType> js(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateJacobian(ip.position(),js);
            std::vector<JacobianType> js_v(lfsv.size());
            lfsv.finiteElement().localBasis().evaluateJacobian(ip.position(),js_v);
#else
            const std::vector<JacobianType>& js = cache[order].evaluateJacobian(ip.position(),lfsu.finiteElement().localBasis());
            const std::vector<JacobianType>& js_v = cache[order].evaluateJacobian(ip.position(),lfsv.finiteElement().localBasis());
#endif

            // transform gradients of shape functions to real element
            jac = eg.geometry().jacobianInverseTransposed(ip.position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            for (size_type i=0; i<lfsu.size(); i++)
              jac.mv(js[i][0],gradphi[i]);

            std::vector<Dune::FieldVector<RF,dim> > gradpsi(lfsv.size());
            for (size_type i=0; i<lfsv.size(); i++)
              jac.mv(js_v[i][0],gradpsi[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (size_type i=0; i<lfsu.size(); i++)
              gradu.axpy(x(lfsu,i),gradphi[i]);

            // compute A * gradient of u
            Dune::FieldVector<RF,dim> Agradu(0.0);
            A.umv(gradu,Agradu);

            // evaluate velocity field
            typename T::Traits::RangeType b = param.b(eg.entity(),ip.position());

            // evaluate reaction term
            typename T::Traits::RangeFieldType c = param.c(eg.entity(),ip.position());

            // integrate (A grad u - bu)*grad phi_i + a*u*phi_i
            RF factor = ip.weight() * eg.geometry().integrationElement(ip.position());
            for (size_type i=0; i<lfsv.size(); i++)
              r.accumulate(lfsv,i,( Agradu*gradpsi[i] - u*(b*gradpsi[i]) + c*u*psi[i] )*factor);
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
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
        const int dim = EG::Entity::dimension;
        const int order = std::max(lfsu.finiteElement().localBasis().order(),
            lfsv.finiteElement().localBasis().order());
        const int intorder = intorderadd + quadrature_factor * order;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        typename T::Traits::PermTensorType A;
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
        A = param.A(eg.entity(),localcenter);

        // transformation
        typename EG::Geometry::JacobianInverseTransposed jac;

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            // evaluate basis functions
#if USECACHE==0
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(ip.position(),phi);
#else
            const std::vector<RangeType>& phi = cache[order].evaluateFunction(ip.position(),lfsu.finiteElement().localBasis());
#endif

            // evaluate gradient of basis functions
#if USECACHE==0
            std::vector<JacobianType> js(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateJacobian(ip.position(),js);
#else
            const std::vector<JacobianType>& js = cache[order].evaluateJacobian(ip.position(),lfsu.finiteElement().localBasis());
#endif

            // transform gradients of shape functions to real element
            jac = eg.geometry().jacobianInverseTransposed(ip.position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            std::vector<Dune::FieldVector<RF,dim> > Agradphi(lfsu.size());
            for (size_type i=0; i<lfsu.size(); i++)
              {
                jac.mv(js[i][0],gradphi[i]);
                A.mv(gradphi[i],Agradphi[i]);
              }

            // evaluate velocity field
            typename T::Traits::RangeType b = param.b(eg.entity(),ip.position());

            // evaluate reaction term
            typename T::Traits::RangeFieldType c = param.c(eg.entity(),ip.position());

            // integrate (A grad u - bu)*grad phi_i + a*u*phi_i
            RF factor = ip.weight() * eg.geometry().integrationElement(ip.position());
            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat.accumulate(lfsu,i,lfsu,j,( Agradphi[j]*gradphi[i] - phi[j]*(b*gradphi[i]) + c*phi[j]*phi[i] )*factor);
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
        const int order = std::max(
            std::max(lfsu_s.finiteElement().localBasis().order(),
                lfsu_n.finiteElement().localBasis().order()),
            std::max(lfsv_s.finiteElement().localBasis().order(),
                lfsv_n.finiteElement().localBasis().order())
            );
        const int intorder = intorderadd+quadrature_factor*order;

        // make copy of inside and outside cell w.r.t. the intersection
        auto inside_cell = ig.inside();
        auto outside_cell = ig.outside();

        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(inside_cell.type()).position(0,0);
        const Dune::FieldVector<DF,dim>&
          outside_local = Dune::ReferenceElements<DF,dim>::general(outside_cell.type()).position(0,0);
        typename T::Traits::PermTensorType A_s, A_n;
        A_s = param.A(inside_cell,inside_local);
        A_n = param.A(outside_cell,outside_local);

        // face diameter; this should be revised for anisotropic meshes?
        DF h_s, h_n;
        DF hmax_s = 0.;
        DF hmax_n = 0.;
        element_size(inside_cell.geometry(),h_s,hmax_s);
        element_size(outside_cell.geometry(),h_n,hmax_n);
        RF h_F = std::min(h_s,h_n);
        h_F = std::min(inside_cell.geometry().volume(),outside_cell.geometry().volume())/ig.geometry().volume(); // Houston!

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // transformation
        typename IG::Entity::Geometry::JacobianInverseTransposed jac;

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
        const int order_s = lfsu_s.finiteElement().localBasis().order();
        const int order_n = lfsu_n.finiteElement().localBasis().order();
        int degree = std::max( order_s, order_n );

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            // exact normal
            const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(ip.position());

            // position of quadrature point in local coordinates of elements
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(ip.position());
            Dune::FieldVector<DF,dim> iplocal_n = ig.geometryInOutside().global(ip.position());

            // evaluate basis functions
#if USECACHE==0
            std::vector<RangeType> phi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);
            std::vector<RangeType> phi_n(lfsu_n.size());
            lfsu_n.finiteElement().localBasis().evaluateFunction(iplocal_n,phi_n);
            std::vector<RangeType> psi_s(lfsv_s.size());
            lfsv_s.finiteElement().localBasis().evaluateFunction(iplocal_s,psi_s);
            std::vector<RangeType> psi_n(lfsv_n.size());
            lfsv_n.finiteElement().localBasis().evaluateFunction(iplocal_n,psi_n);
#else
            const std::vector<RangeType>& phi_s = cache[order_s].evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());
            const std::vector<RangeType>& phi_n = cache[order_n].evaluateFunction(iplocal_n,lfsu_n.finiteElement().localBasis());
            const std::vector<RangeType>& psi_s = cache[order_s].evaluateFunction(iplocal_s,lfsv_s.finiteElement().localBasis());
            const std::vector<RangeType>& psi_n = cache[order_n].evaluateFunction(iplocal_n,lfsv_n.finiteElement().localBasis());
#endif

            // evaluate u
            RF u_s=0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              u_s += x_s(lfsu_s,i)*phi_s[i];
            RF u_n=0.0;
            for (size_type i=0; i<lfsu_n.size(); i++)
              u_n += x_n(lfsu_n,i)*phi_n[i];

            // evaluate gradient of basis functions
#if USECACHE==0
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
            std::vector<JacobianType> gradphi_n(lfsu_n.size());
            lfsu_n.finiteElement().localBasis().evaluateJacobian(iplocal_n,gradphi_n);
            std::vector<JacobianType> gradpsi_s(lfsv_s.size());
            lfsv_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradpsi_s);
            std::vector<JacobianType> gradpsi_n(lfsv_n.size());
            lfsv_n.finiteElement().localBasis().evaluateJacobian(iplocal_n,gradpsi_n);
#else
            const std::vector<JacobianType>& gradphi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
            const std::vector<JacobianType>& gradphi_n = cache[order_n].evaluateJacobian(iplocal_n,lfsu_n.finiteElement().localBasis());
            const std::vector<JacobianType>& gradpsi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsv_s.finiteElement().localBasis());
            const std::vector<JacobianType>& gradpsi_n = cache[order_n].evaluateJacobian(iplocal_n,lfsv_n.finiteElement().localBasis());
#endif

            // transform gradients of shape functions to real element
            jac = inside_cell.geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
            std::vector<Dune::FieldVector<RF,dim> > tgradpsi_s(lfsv_s.size());
            for (size_type i=0; i<lfsv_s.size(); i++) jac.mv(gradpsi_s[i][0],tgradpsi_s[i]);
            jac = outside_cell.geometry().jacobianInverseTransposed(iplocal_n);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_n(lfsu_n.size());
            for (size_type i=0; i<lfsu_n.size(); i++) jac.mv(gradphi_n[i][0],tgradphi_n[i]);
            std::vector<Dune::FieldVector<RF,dim> > tgradpsi_n(lfsv_n.size());
            for (size_type i=0; i<lfsv_n.size(); i++) jac.mv(gradpsi_n[i][0],tgradpsi_n[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu_s(0.0);
            for (size_type i=0; i<lfsu_s.size(); i++)
              gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);
            Dune::FieldVector<RF,dim> gradu_n(0.0);
            for (size_type i=0; i<lfsu_n.size(); i++)
              gradu_n.axpy(x_n(lfsu_n,i),tgradphi_n[i]);

            // evaluate velocity field and upwinding, assume H(div) velocity field => may choose any side
            typename T::Traits::RangeType b = param.b(inside_cell,iplocal_s);
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
            RF factor = ip.weight() * ig.geometry().integrationElement(ip.position());

            // convection term
            RF term1 = (omegaup_s*u_s + omegaup_n*u_n) * normalflux *factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsv_s,i,term1 * psi_s[i]);
            for (size_type i=0; i<lfsv_n.size(); i++)
              r_n.accumulate(lfsv_n,i,-term1 * psi_n[i]);

            // diffusion term
            RF term2 =  -(omega_s*(An_F_s*gradu_s) + omega_n*(An_F_n*gradu_n)) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsv_s,i,term2 * psi_s[i]);
            for (size_type i=0; i<lfsv_n.size(); i++)
              r_n.accumulate(lfsv_n,i,-term2 * psi_n[i]);

            // (non-)symmetric IP term
            RF term3 = (u_s-u_n) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsv_s,i,term3 * theta * omega_s * (An_F_s*tgradpsi_s[i]));
            for (size_type i=0; i<lfsv_n.size(); i++)
              r_n.accumulate(lfsv_n,i,term3 * theta * omega_n * (An_F_n*tgradpsi_n[i]));

            // standard IP term integral
            RF term4 = penalty_factor * (u_s-u_n) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsv_s,i,term4 * psi_s[i]);
            for (size_type i=0; i<lfsv_n.size(); i++)
              r_n.accumulate(lfsv_n,i,-term4 * psi_n[i]);
          }
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                              M& mat_ss, M& mat_sn,
                              M& mat_ns, M& mat_nn) const
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
        const int order = std::max(
            std::max(lfsu_s.finiteElement().localBasis().order(),
                lfsu_n.finiteElement().localBasis().order()),
            std::max(lfsv_s.finiteElement().localBasis().order(),
                lfsv_n.finiteElement().localBasis().order())
            );
        const int intorder = intorderadd+quadrature_factor*order;

        // make copy of inside and outside cell w.r.t. the intersection
        auto inside_cell = ig.inside();
        auto outside_cell = ig.outside();

        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(inside_cell.type()).position(0,0);
        const Dune::FieldVector<DF,dim>&
          outside_local = Dune::ReferenceElements<DF,dim>::general(outside_cell.type()).position(0,0);
        typename T::Traits::PermTensorType A_s, A_n;
        A_s = param.A(inside_cell,inside_local);
        A_n = param.A(outside_cell,outside_local);

        // face diameter; this should be revised for anisotropic meshes?
        DF h_s, h_n;
        DF hmax_s = 0., hmax_n = 0.;
        element_size(inside_cell.geometry(),h_s,hmax_s);
        element_size(outside_cell.geometry(),h_n,hmax_n);
        RF h_F = std::min(h_s,h_n);
        h_F = std::min(inside_cell.geometry().volume(),outside_cell.geometry().volume())/ig.geometry().volume(); // Houston!

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // transformation
        typename IG::Entity::Geometry::JacobianInverseTransposed jac;

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
        const int order_s = lfsu_s.finiteElement().localBasis().order();
        const int order_n = lfsu_n.finiteElement().localBasis().order();
        int degree = std::max( order_s, order_n );

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            // exact normal
            const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(ip.position());

            // position of quadrature point in local coordinates of elements
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(ip.position());
            Dune::FieldVector<DF,dim> iplocal_n = ig.geometryInOutside().global(ip.position());

            // evaluate basis functions
#if USECACHE==0
            std::vector<RangeType> phi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);
            std::vector<RangeType> phi_n(lfsu_n.size());
            lfsu_n.finiteElement().localBasis().evaluateFunction(iplocal_n,phi_n);
#else
            const std::vector<RangeType>& phi_s = cache[order_s].evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());
            const std::vector<RangeType>& phi_n = cache[order_n].evaluateFunction(iplocal_n,lfsu_n.finiteElement().localBasis());
#endif

            // evaluate gradient of basis functions
#if USECACHE==0
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
            std::vector<JacobianType> gradphi_n(lfsu_n.size());
            lfsu_n.finiteElement().localBasis().evaluateJacobian(iplocal_n,gradphi_n);
#else
            const std::vector<JacobianType>& gradphi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
            const std::vector<JacobianType>& gradphi_n = cache[order_n].evaluateJacobian(iplocal_n,lfsu_n.finiteElement().localBasis());
#endif

            // transform gradients of shape functions to real element
            jac = inside_cell.geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
            jac = outside_cell.geometry().jacobianInverseTransposed(iplocal_n);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_n(lfsu_n.size());
            for (size_type i=0; i<lfsu_n.size(); i++) jac.mv(gradphi_n[i][0],tgradphi_n[i]);

            // evaluate velocity field and upwinding, assume H(div) velocity field => may choose any side
            typename T::Traits::RangeType b = param.b(inside_cell,iplocal_s);
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
            RF factor = ip.weight() * ig.geometry().integrationElement(ip.position());
            RF ipfactor = penalty_factor * factor;

            // do all terms in the order: I convection, II diffusion, III consistency, IV ip
            for (size_type j=0; j<lfsu_s.size(); j++) {
              RF temp1 = -(An_F_s*tgradphi_s[j])*omega_s*factor;
              for (size_type i=0; i<lfsu_s.size(); i++) {
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,omegaup_s * phi_s[j] * normalflux *factor * phi_s[i]);
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,temp1 * phi_s[i]);
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * factor * theta * omega_s * (An_F_s*tgradphi_s[i]));
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * ipfactor * phi_s[i]);
              }
            }
            for (size_type j=0; j<lfsu_n.size(); j++) {
              RF temp1 = -(An_F_n*tgradphi_n[j])*omega_n*factor;
              for (size_type i=0; i<lfsu_s.size(); i++) {
                mat_sn.accumulate(lfsu_s,i,lfsu_n,j,omegaup_n * phi_n[j] * normalflux *factor * phi_s[i]);
                mat_sn.accumulate(lfsu_s,i,lfsu_n,j,temp1 * phi_s[i]);
                mat_sn.accumulate(lfsu_s,i,lfsu_n,j,-phi_n[j] * factor * theta * omega_s * (An_F_s*tgradphi_s[i]));
                mat_sn.accumulate(lfsu_s,i,lfsu_n,j,-phi_n[j] * ipfactor * phi_s[i]);
              }
            }
            for (size_type j=0; j<lfsu_s.size(); j++) {
              RF temp1 = -(An_F_s*tgradphi_s[j])*omega_s*factor;
              for (size_type i=0; i<lfsu_n.size(); i++) {
                mat_ns.accumulate(lfsu_n,i,lfsu_s,j,-omegaup_s * phi_s[j] * normalflux *factor * phi_n[i]);
                mat_ns.accumulate(lfsu_n,i,lfsu_s,j,-temp1 * phi_n[i]);
                mat_ns.accumulate(lfsu_n,i,lfsu_s,j,phi_s[j] * factor * theta * omega_n * (An_F_n*tgradphi_n[i]));
                mat_ns.accumulate(lfsu_n,i,lfsu_s,j,-phi_s[j] * ipfactor * phi_n[i]);
              }
            }
            for (size_type j=0; j<lfsu_n.size(); j++) {
              RF temp1 = -(An_F_n*tgradphi_n[j])*omega_n*factor;
              for (size_type i=0; i<lfsu_n.size(); i++) {
                mat_nn.accumulate(lfsu_n,i,lfsu_n,j,-omegaup_n * phi_n[j] * normalflux *factor * phi_n[i]);
                mat_nn.accumulate(lfsu_n,i,lfsu_n,j,-temp1 * phi_n[i]);
                mat_nn.accumulate(lfsu_n,i,lfsu_n,j,-phi_n[j] * factor * theta * omega_n * (An_F_n*tgradphi_n[i]));
                mat_nn.accumulate(lfsu_n,i,lfsu_n,j,phi_n[j] * ipfactor * phi_n[i]);
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
        const int order = std::max(
            lfsu_s.finiteElement().localBasis().order(),
            lfsv_s.finiteElement().localBasis().order()
            );
        const int intorder = intorderadd+quadrature_factor*order;

        // make copy of inside cell w.r.t. the boundary
        auto inside_cell = ig.inside();

        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(inside_cell.type()).position(0,0);
        typename T::Traits::PermTensorType A_s;
        A_s = param.A(inside_cell,inside_local);

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // face diameter
        DF h_s;
        DF hmax_s = 0.;
        element_size(inside_cell.geometry(),h_s,hmax_s);
        RF h_F = h_s;
        h_F = inside_cell.geometry().volume()/ig.geometry().volume(); // Houston!

        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

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
        const int order_s = lfsu_s.finiteElement().localBasis().order();
        int degree = order_s;

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            BCType bctype = param.bctype(ig.intersection(),ip.position());

            if (bctype == ConvectionDiffusionBoundaryConditions::None)
              continue;

            // position of quadrature point in local coordinates of elements
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(ip.position());

            // local normal
            const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(ip.position());

            // evaluate basis functions
#if USECACHE==0
            std::vector<RangeType> phi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);
            std::vector<RangeType> psi_s(lfsv_s.size());
            lfsv_s.finiteElement().localBasis().evaluateFunction(iplocal_s,psi_s);
#else
            const std::vector<RangeType>& phi_s = cache[order_s].evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());
            const std::vector<RangeType>& psi_s = cache[order_s].evaluateFunction(iplocal_s,lfsv_s.finiteElement().localBasis());
#endif

            // integration factor
            RF factor = ip.weight() * ig.geometry().integrationElement(ip.position());

            if (bctype == ConvectionDiffusionBoundaryConditions::Neumann)
              {
                // evaluate flux boundary condition
                RF j = param.j(ig.intersection(),ip.position());

                // integrate
                for (size_type i=0; i<lfsv_s.size(); i++)
                  r_s.accumulate(lfsv_s,i,j * psi_s[i] * factor);

                continue;
              }

            // evaluate u
            RF u_s=0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              u_s += x_s(lfsu_s,i)*phi_s[i];

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            typename T::Traits::RangeType b = param.b(inside_cell,iplocal_s);
            RF normalflux = b*n_F_local;

            if (bctype == ConvectionDiffusionBoundaryConditions::Outflow)
              {
                if (normalflux<-1e-30)
                  DUNE_THROW(Dune::Exception,
                    "Outflow boundary condition on inflow! [b("
                    << ig.geometry().global(ip.position()) << ") = "
                    << b << ")");

                // convection term
                RF term1 = u_s * normalflux *factor;
                for (size_type i=0; i<lfsv_s.size(); i++)
                  r_s.accumulate(lfsv_s,i,term1 * psi_s[i]);

                // evaluate flux boundary condition
                RF o = param.o(ig.intersection(),ip.position());

                // integrate
                for (size_type i=0; i<lfsv_s.size(); i++)
                  r_s.accumulate(lfsv_s,i,o * psi_s[i] * factor);

                continue;
              }

            // evaluate gradient of basis functions
            assert (bctype == ConvectionDiffusionBoundaryConditions::Dirichlet);
#if USECACHE==0
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
            std::vector<JacobianType> gradpsi_s(lfsv_s.size());
            lfsv_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradpsi_s);
#else
            const std::vector<JacobianType>& gradphi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
            const std::vector<JacobianType>& gradpsi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsv_s.finiteElement().localBasis());
#endif

            // transform gradients of shape functions to real element
            jac = inside_cell.geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
            std::vector<Dune::FieldVector<RF,dim> > tgradpsi_s(lfsv_s.size());
            for (size_type i=0; i<lfsv_s.size(); i++) jac.mv(gradpsi_s[i][0],tgradpsi_s[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu_s(0.0);
            for (size_type i=0; i<lfsu_s.size(); i++)
              gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);

            // evaluate Dirichlet boundary condition
            RF g = param.g(inside_cell,iplocal_s);

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
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsv_s,i,term1 * psi_s[i]);

            // diffusion term
            RF term2 =  (An_F_s*gradu_s) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsv_s,i,-term2 * psi_s[i]);

            // (non-)symmetric IP term
            RF term3 = (u_s-g) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsv_s,i,term3 * theta * (An_F_s*tgradpsi_s[i]));

            // standard IP term
            RF term4 = penalty_factor * (u_s-g) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsv_s,i,term4 * psi_s[i]);
          }
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              M& mat_ss) const
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
        const int order = std::max(
            lfsu_s.finiteElement().localBasis().order(),
            lfsv_s.finiteElement().localBasis().order()
            );
        const int intorder = intorderadd+quadrature_factor*order;

        // make copy of inside cell w.r.t. the boundary
        auto inside_cell = ig.inside();

        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(inside_cell.type()).position(0,0);
        typename T::Traits::PermTensorType A_s;
        A_s = param.A(inside_cell,inside_local);

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // face diameter
        DF h_s;
        DF hmax_s = 0.;
        element_size(inside_cell.geometry(),h_s,hmax_s);
        RF h_F = h_s;
        h_F = inside_cell.geometry().volume()/ig.geometry().volume(); // Houston!

        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

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
        const int order_s = lfsu_s.finiteElement().localBasis().order();
        int degree = order_s;

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // Neumann boundary makes no contribution to boundary
        //if (bctype == ConvectionDiffusionBoundaryConditions::Neumann) return;

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            BCType bctype = param.bctype(ig.intersection(),ip.position());

            if (bctype == ConvectionDiffusionBoundaryConditions::None ||
                bctype == ConvectionDiffusionBoundaryConditions::Neumann)
              continue;

            // position of quadrature point in local coordinates of elements
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(ip.position());

            // local normal
            const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(ip.position());

            // evaluate basis functions
#if USECACHE==0
            std::vector<RangeType> phi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);
#else
            const std::vector<RangeType>& phi_s = cache[order_s].evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());
#endif

            // integration factor
            RF factor = ip.weight() * ig.geometry().integrationElement(ip.position());

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            typename T::Traits::RangeType b = param.b(inside_cell,iplocal_s);
            RF normalflux = b*n_F_local;

            if (bctype == ConvectionDiffusionBoundaryConditions::Outflow)
              {
                if (normalflux<-1e-30)
                  DUNE_THROW(Dune::Exception,
                    "Outflow boundary condition on inflow! [b("
                    << ig.geometry().global(ip.position()) << ") = "
                    << b << ")" << n_F_local << " " << normalflux);

                // convection term
                for (size_type j=0; j<lfsu_s.size(); j++)
                  for (size_type i=0; i<lfsu_s.size(); i++)
                    mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * normalflux * factor * phi_s[i]);

                continue;
              }

            // evaluate gradient of basis functions
#if USECACHE==0
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
#else
            const std::vector<JacobianType>& gradphi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
#endif

            // transform gradients of shape functions to real element
            jac = inside_cell.geometry().jacobianInverseTransposed(iplocal_s);
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
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,omegaup_s * phi_s[j] * normalflux * factor * phi_s[i]);

            // diffusion term
            for (size_type j=0; j<lfsu_s.size(); j++)
              for (size_type i=0; i<lfsu_s.size(); i++)
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,-(An_F_s*tgradphi_s[j]) * factor * phi_s[i]);

            // (non-)symmetric IP term
            for (size_type j=0; j<lfsu_s.size(); j++)
              for (size_type i=0; i<lfsu_s.size(); i++)
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * factor * theta * (An_F_s*tgradphi_s[i]));

            // standard IP term
            for (size_type j=0; j<lfsu_s.size(); j++)
              for (size_type i=0; i<lfsu_s.size(); i++)
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,penalty_factor * phi_s[j] * phi_s[i] * factor);
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
        const int dim = EG::Entity::dimension;
        const int order = lfsv.finiteElement().localBasis().order();
        const int intorder = intorderadd + 2 * order;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (const auto& ip : rule)
          {
            // evaluate shape functions
#if USECACHE==0
            std::vector<RangeType> phi(lfsv.size());
            lfsv.finiteElement().localBasis().evaluateFunction(ip.position(),phi);
#else
            const std::vector<RangeType>& phi = cache[order].evaluateFunction(ip.position(),lfsv.finiteElement().localBasis());
#endif

            // evaluate right hand side parameter function
            Real f;
            f = param.f(eg.entity(),ip.position());

            // integrate f
            RF factor = ip.weight() * eg.geometry().integrationElement(ip.position());
            for (size_type i=0; i<lfsv.size(); i++)
              r.accumulate(lfsv,i,-f*phi[i]*factor);
          }
      }

      //! set time in parameter class
      void setTime (double t)
      {
        Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>::setTime(t);
        param.setTime(t);
      }

    private:
      T& param;  // two phase parameter class
      ConvectionDiffusionDGMethod::Type method;
      ConvectionDiffusionDGWeights::Type weights;
      Real alpha, beta;
      int intorderadd;
      int quadrature_factor;
      Real theta;
      typedef typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;

      typedef Dune::PDELab::LocalBasisCache<LocalBasisType> Cache;

      // In theory it is possible that one and the same local operator is
      // called first with a finite element of one type and later with a
      // finite element of another type.  Since finite elements of different
      // type will usually produce different results for the same local
      // coordinate they cannot share a cache.  Here we use a vector of caches
      // to allow for different orders of the shape functions, which should be
      // enough to support p-adaptivity.  (Another likely candidate would be
      // differing geometry types, i.e. hybrid meshes.)

      std::vector<Cache> cache;

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
            for (int i=0; i<Dune::ReferenceElements<DF,dim>::general(gt).size(dim-1); i++)
              {
                Dune::FieldVector<DF,dim> x = geo.corner(Dune::ReferenceElements<DF,dim>::general(gt).subEntity(i,dim-1,0,dim));
                x -= geo.corner(Dune::ReferenceElements<DF,dim>::general(gt).subEntity(i,dim-1,1,dim));
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
