// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONFASTDG_HH
#define DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONFASTDG_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

/**
   \todo update quadrature order to work with lfsv != lfsu
   \todo update alpha_* to work with lfsv != lfsu (./)
   \todo update jacobian_* to work with lfsv != lfsu
   \todo update caches to work with lfsv != lfsu
 */
namespace Dune {
  namespace PDELab {

    /** a local operator for solving the convection-diffusion equation with discontinuous Galerkin
     *
     * \f{align*}{
     *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \\
     *                                              u &=& g \mbox{ on } \partial\Omega_D \\
     *                (b(x) u - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
     *                        -(A(x)\nabla u) \cdot n &=& o \mbox{ on } \partial\Omega_O
     * \f}
     * Note:
     *  - This formulation is valid for velocity fields which are non-divergence free.
     *  - Outflow boundary conditions should only be set on the outflow boundary
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T, typename FiniteElementMap>
    class ConvectionDiffusionFastDG
      : public Dune::PDELab::NumericalJacobianApplyVolume<ConvectionDiffusionFastDG<T,FiniteElementMap> >,
        public Dune::PDELab::NumericalJacobianApplySkeleton<ConvectionDiffusionFastDG<T,FiniteElementMap> >,
        public Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionFastDG<T,FiniteElementMap> >,
        public Dune::PDELab::FullSkeletonPattern,
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::LocalOperatorDefaultFlags,
        public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>
    {
      enum { dim = T::Traits::GridViewType::dimension };

      using Real = typename T::Traits::RangeFieldType;
      using BCType = typename ConvectionDiffusionBoundaryConditions::Type;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume  = true };

      /** \brief constructor: pass parameter object and define DG-method
       * \param[in] param_   Reference to parameter object.
       * \param[in] method_  Interior penalty Galerkin method. Default is skew-symmetric.
       * \param[in] weights_ Weighted averages for diffusion tensor. Default is no weighting.
       * \param[in] alpha_   Penalization constant. Default is zero.
       *
       * Collecting the input parameters above, the default is the OBB-method.
       */
      ConvectionDiffusionFastDG (T& param_,
                                 ConvectionDiffusionDGMethod::Type method_=ConvectionDiffusionDGMethod::NIPG,
                                 ConvectionDiffusionDGWeights::Type weights_=ConvectionDiffusionDGWeights::weightsOff,
                                 Real alpha_=0.0,
                                 int intorderadd_=0)
        : Dune::PDELab::NumericalJacobianApplyVolume<ConvectionDiffusionFastDG<T,FiniteElementMap> >(1.0e-7),
          Dune::PDELab::NumericalJacobianApplySkeleton<ConvectionDiffusionFastDG<T,FiniteElementMap> >(1.0e-7),
          Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionFastDG<T,FiniteElementMap> >(1.0e-7),
          param(param_), method(method_), weights(weights_),
          alpha(alpha_), intorderadd(intorderadd_), quadrature_factor(2),
          cache(20)
      {
        theta = 1.0;
        if (method==ConvectionDiffusionDGMethod::SIPG) theta = -1.0;
        if (method==ConvectionDiffusionDGMethod::IIPG) theta = 0.0;
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSU::Traits::SizeType;

        // dimensions
        const int dim = EG::Entity::dimension;
        const int order = std::max(lfsu.finiteElement().localBasis().order(),
            lfsv.finiteElement().localBasis().order());

        // Get cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        auto ref_el = referenceElement(geo);
        auto localcenter = ref_el.position(0,0);
        auto A = param.A(cell,localcenter);

        // Initialize vectors outside for loop
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
        std::vector<Dune::FieldVector<RF,dim> > gradpsi(lfsv.size());
        Dune::FieldVector<RF,dim> gradu(0.0);
        Dune::FieldVector<RF,dim> Agradu(0.0);

        // Transformation matrix
        typename EG::Geometry::JacobianInverseTransposed jac;

        // Weight of time stepping
        auto t_weight = r.weight();

        // loop over quadrature points
        auto intorder = intorderadd + quadrature_factor * order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            auto& phi = cache[order].evaluateFunction(ip.position(),lfsu.finiteElement().localBasis());
            auto& psi = cache[order].evaluateFunction(ip.position(),lfsv.finiteElement().localBasis());

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x[i]*phi[i];

            // evaluate gradient of basis functions
            auto& js = cache[order].evaluateJacobian(ip.position(),lfsu.finiteElement().localBasis());
            auto& js_v = cache[order].evaluateJacobian(ip.position(),lfsv.finiteElement().localBasis());

            // transform gradients of shape functions to real element
            jac = geo.jacobianInverseTransposed(ip.position());
            for (size_type i=0; i<lfsu.size(); i++)
              jac.mv(js[i][0],gradphi[i]);

            for (size_type i=0; i<lfsv.size(); i++)
              jac.mv(js_v[i][0],gradpsi[i]);

            // compute gradient of u
            gradu = 0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              gradu.axpy(x[i],gradphi[i]);

            // compute A * gradient of u
            A.mv(gradu,Agradu);

            // evaluate velocity field
            auto b = param.b(cell,ip.position());

            // evaluate reaction term
            auto c = param.c(cell,ip.position());

            // integrate (A grad u - bu)*grad phi_i + a*u*phi_i
            RF factor = ip.weight() * geo.integrationElement(ip.position()) * t_weight;
            for (size_type i=0; i<lfsv.size(); i++)
              r.data()[i] += ( Agradu*gradpsi[i] - u*(b*gradpsi[i]) + c*u*psi[i] )*factor;
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSU::Traits::SizeType;

        // dimensions
        const int dim = EG::Entity::dimension;
        const int order = std::max(lfsu.finiteElement().localBasis().order(),
            lfsv.finiteElement().localBasis().order());

        // Get cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        auto ref_el = referenceElement(geo);
        auto localcenter = ref_el.position(0,0);
        auto A = param.A(cell,localcenter);

        // Initialize vectors outside for loop
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
        std::vector<Dune::FieldVector<RF,dim> > Agradphi(lfsu.size());

        // Transformation matrix
        typename EG::Geometry::JacobianInverseTransposed jac;

        // Weight of time stepping
        auto t_weight = mat.weight();

        // loop over quadrature points
        auto intorder = intorderadd + quadrature_factor * order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            auto& phi = cache[order].evaluateFunction(ip.position(),lfsu.finiteElement().localBasis());

            // evaluate gradient of basis functions
            auto& js = cache[order].evaluateJacobian(ip.position(),lfsu.finiteElement().localBasis());

            // transform gradients of shape functions to real element
            jac = geo.jacobianInverseTransposed(ip.position());
            for (size_type i=0; i<lfsu.size(); i++)
              {
                jac.mv(js[i][0],gradphi[i]);
                A.mv(gradphi[i],Agradphi[i]);
              }

            // evaluate velocity field
            auto b = param.b(cell,ip.position());

            // evaluate reaction term
            auto c = param.c(cell,ip.position());

            // integrate (A grad u - bu)*grad phi_i + a*u*phi_i
            auto factor = ip.weight() * geo.integrationElement(ip.position()) * t_weight;
            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat.data()[lfsu.size()*i+j] += (( Agradphi[j]*gradphi[i] - phi[j]*(b*gradphi[i]) + c*phi[j]*phi[i] )*factor);
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
        // define types
        using RF = typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSV::Traits::SizeType;

        // dimensions
        const int dim = IG::dimension;
        const int order = std::max(
            std::max(lfsu_s.finiteElement().localBasis().order(),
                lfsu_n.finiteElement().localBasis().order()),
            std::max(lfsv_s.finiteElement().localBasis().order(),
                lfsv_n.finiteElement().localBasis().order())
            );

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

        // evaluate permeability tensors
        auto ref_el_inside = referenceElement(geo_inside);
        auto ref_el_outside = referenceElement(geo_outside);
        auto local_inside = ref_el_inside.position(0,0);
        auto local_outside = ref_el_outside.position(0,0);
        auto A_s = param.A(cell_inside,local_inside);
        auto A_n = param.A(cell_outside,local_outside);

        // face diameter for anisotropic meshes taken from Paul Houston et al.
        // this formula ensures coercivity of the bilinear form
        auto h_F = std::min(geo_inside.volume(),geo_outside.volume())/geo.volume();

        // tensor times normal
        auto n_F = ig.centerUnitOuterNormal();
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
        auto order_s = lfsu_s.finiteElement().localBasis().order();
        auto order_n = lfsu_n.finiteElement().localBasis().order();
        auto degree = std::max( order_s, order_n );

        // penalty factor
        auto penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // Initialize vectors outside for loop
        std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
        std::vector<Dune::FieldVector<RF,dim> > tgradpsi_s(lfsv_s.size());
        std::vector<Dune::FieldVector<RF,dim> > tgradphi_n(lfsu_n.size());
        std::vector<Dune::FieldVector<RF,dim> > tgradpsi_n(lfsv_n.size());
        Dune::FieldVector<RF,dim> gradu_s(0.0);
        Dune::FieldVector<RF,dim> gradu_n(0.0);

        // Transformation matrix
        typename IG::Entity::Geometry::JacobianInverseTransposed jac;

        // Weight of time stepping
        auto t_weight = r_s.weight();
        assert(r_s.weight()==r_n.weight());

        // loop over quadrature points
        auto intorder = intorderadd+quadrature_factor*order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // exact normal
            auto n_F_local = ig.unitOuterNormal(ip.position());

            // position of quadrature point in local coordinates of elements
            auto iplocal_s = geo_in_inside.global(ip.position());
            auto iplocal_n = geo_in_outside.global(ip.position());

            // evaluate basis functions
            auto& phi_s = cache[order_s].evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());
            auto& phi_n = cache[order_n].evaluateFunction(iplocal_n,lfsu_n.finiteElement().localBasis());
            auto& psi_s = cache[order_s].evaluateFunction(iplocal_s,lfsv_s.finiteElement().localBasis());
            auto& psi_n = cache[order_n].evaluateFunction(iplocal_n,lfsv_n.finiteElement().localBasis());

            // evaluate u
            RF u_s=0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              u_s += x_s[i]*phi_s[i];
            RF u_n=0.0;
            for (size_type i=0; i<lfsu_n.size(); i++)
              u_n += x_n[i]*phi_n[i];

            // evaluate gradient of basis functions
            auto& gradphi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
            auto& gradphi_n = cache[order_n].evaluateJacobian(iplocal_n,lfsu_n.finiteElement().localBasis());
            auto& gradpsi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsv_s.finiteElement().localBasis());
            auto& gradpsi_n = cache[order_n].evaluateJacobian(iplocal_n,lfsv_n.finiteElement().localBasis());

            // transform gradients of shape functions to real element
            jac = geo_inside.jacobianInverseTransposed(iplocal_s);
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
            for (size_type i=0; i<lfsv_s.size(); i++) jac.mv(gradpsi_s[i][0],tgradpsi_s[i]);
            jac = geo_outside.jacobianInverseTransposed(iplocal_n);
            for (size_type i=0; i<lfsu_n.size(); i++) jac.mv(gradphi_n[i][0],tgradphi_n[i]);
            for (size_type i=0; i<lfsv_n.size(); i++) jac.mv(gradpsi_n[i][0],tgradpsi_n[i]);

            // compute gradient of u
            gradu_s = 0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              gradu_s.axpy(x_s[i],tgradphi_s[i]);
            gradu_n = 0.0;
            for (size_type i=0; i<lfsu_n.size(); i++)
              gradu_n.axpy(x_n[i],tgradphi_n[i]);

            // evaluate velocity field and upwinding, assume H(div) velocity field => may choose any side
            auto b = param.b(cell_inside,iplocal_s);
            auto normalflux = b*n_F_local;
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
            auto factor = ip.weight()*geo.integrationElement(ip.position()) * t_weight;

            // convection term
            auto term1 = (omegaup_s*u_s + omegaup_n*u_n) * normalflux *factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.data()[i] += (term1 * psi_s[i]);
            for (size_type i=0; i<lfsv_n.size(); i++)
              r_n.data()[i] += (-term1 * psi_n[i]);

            // diffusion term
            auto term2 =  -(omega_s*(An_F_s*gradu_s) + omega_n*(An_F_n*gradu_n)) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.data()[i] += (term2 * psi_s[i]);
            for (size_type i=0; i<lfsv_n.size(); i++)
              r_n.data()[i] += (-term2 * psi_n[i]);

            // (non-)symmetric IP term
            auto term3 = (u_s-u_n) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.data()[i] += (term3 * theta * omega_s * (An_F_s*tgradpsi_s[i]));
            for (size_type i=0; i<lfsv_n.size(); i++)
              r_n.data()[i] += (term3 * theta * omega_n * (An_F_n*tgradpsi_n[i]));

            // standard IP term integral
            auto term4 = penalty_factor * (u_s-u_n) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.data()[i] += (term4 * psi_s[i]);
            for (size_type i=0; i<lfsv_n.size(); i++)
              r_n.data()[i] += (-term4 * psi_n[i]);
          }
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                              M& mat_ss, M& mat_sn,
                              M& mat_ns, M& mat_nn) const
      {
        // define types
        using RF = typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSV::Traits::SizeType;

        // dimensions
        const int dim = IG::dimension;
        const int order = std::max(
            std::max(lfsu_s.finiteElement().localBasis().order(),
                lfsu_n.finiteElement().localBasis().order()),
            std::max(lfsv_s.finiteElement().localBasis().order(),
                lfsv_n.finiteElement().localBasis().order())
            );

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

        // evaluate permeability tensors
        auto ref_el_inside = referenceElement(geo_inside);
        auto ref_el_outside = referenceElement(geo_outside);
        auto local_inside = ref_el_inside.position(0,0);
        auto local_outside = ref_el_outside.position(0,0);
        auto A_s = param.A(cell_inside,local_inside);
        auto A_n = param.A(cell_outside,local_outside);

        // face diameter for anisotropic meshes taken from Paul Houston et al.
        // this formula ensures coercivity of the bilinear form
        auto h_F = std::min(geo_inside.volume(),geo_outside.volume())/geo.volume();

        // tensor times normal
        auto n_F = ig.centerUnitOuterNormal();
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
        auto order_s = lfsu_s.finiteElement().localBasis().order();
        auto order_n = lfsu_n.finiteElement().localBasis().order();
        auto degree = std::max( order_s, order_n );

        // penalty factor
        auto penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // Initialize vectors outside for loop
        std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
        std::vector<Dune::FieldVector<RF,dim> > tgradphi_n(lfsu_n.size());

        // Transformation matrix
        typename IG::Entity::Geometry::JacobianInverseTransposed jac;

        // Weight of time stepping
        auto t_weight = mat_ss.weight();
        assert(mat_ss.weight()==mat_sn.weight());
        assert(mat_ss.weight()==mat_ns.weight());
        assert(mat_ss.weight()==mat_nn.weight());

        // loop over quadrature points
        const int intorder = intorderadd+quadrature_factor*order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // exact normal
            auto n_F_local = ig.unitOuterNormal(ip.position());

            // position of quadrature point in local coordinates of elements
            auto iplocal_s = geo_in_inside.global(ip.position());
            auto iplocal_n = geo_in_outside.global(ip.position());

            // evaluate basis functions
            auto& phi_s = cache[order_s].evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());
            auto& phi_n = cache[order_n].evaluateFunction(iplocal_n,lfsu_n.finiteElement().localBasis());

            // evaluate gradient of basis functions
            auto& gradphi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
            auto& gradphi_n = cache[order_n].evaluateJacobian(iplocal_n,lfsu_n.finiteElement().localBasis());

            // transform gradients of shape functions to real element
            jac = geo_inside.jacobianInverseTransposed(iplocal_s);
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
            jac = geo_outside.jacobianInverseTransposed(iplocal_n);
            for (size_type i=0; i<lfsu_n.size(); i++) jac.mv(gradphi_n[i][0],tgradphi_n[i]);

            // evaluate velocity field and upwinding, assume H(div) velocity field => may choose any side
            auto b = param.b(cell_inside,iplocal_s);
            auto normalflux = b*n_F_local;
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
            auto factor = ip.weight() * geo.integrationElement(ip.position()) * t_weight;
            auto ipfactor = penalty_factor * factor;

            // do all terms in the order: I convection, II diffusion, III consistency, IV ip
            for (size_type j=0; j<lfsu_s.size(); j++) {
              auto temp1 = -(An_F_s*tgradphi_s[j])*omega_s*factor;
              for (size_type i=0; i<lfsu_s.size(); i++) {
                mat_ss.data()[lfsu_s.size()*i+j] += (omegaup_s * phi_s[j] * normalflux *factor * phi_s[i]);
                mat_ss.data()[lfsu_s.size()*i+j] += (temp1 * phi_s[i]);
                mat_ss.data()[lfsu_s.size()*i+j] += (phi_s[j] * factor * theta * omega_s * (An_F_s*tgradphi_s[i]));
                mat_ss.data()[lfsu_s.size()*i+j] += (phi_s[j] * ipfactor * phi_s[i]);
              }
            }
            for (size_type j=0; j<lfsu_n.size(); j++) {
              auto temp1 = -(An_F_n*tgradphi_n[j])*omega_n*factor;
              for (size_type i=0; i<lfsu_s.size(); i++) {
                mat_sn.data()[lfsu_s.size()*i+j] += (omegaup_n * phi_n[j] * normalflux *factor * phi_s[i]);
                mat_sn.data()[lfsu_s.size()*i+j] += (temp1 * phi_s[i]);
                mat_sn.data()[lfsu_s.size()*i+j] += (-phi_n[j] * factor * theta * omega_s * (An_F_s*tgradphi_s[i]));
                mat_sn.data()[lfsu_s.size()*i+j] += (-phi_n[j] * ipfactor * phi_s[i]);
              }
            }
            for (size_type j=0; j<lfsu_s.size(); j++) {
              auto temp1 = -(An_F_s*tgradphi_s[j])*omega_s*factor;
              for (size_type i=0; i<lfsu_n.size(); i++) {
                mat_ns.data()[lfsu_n.size()*i+j] += (-omegaup_s * phi_s[j] * normalflux *factor * phi_n[i]);
                mat_ns.data()[lfsu_n.size()*i+j] += (-temp1 * phi_n[i]);
                mat_ns.data()[lfsu_n.size()*i+j] += (phi_s[j] * factor * theta * omega_n * (An_F_n*tgradphi_n[i]));
                mat_ns.data()[lfsu_n.size()*i+j] += (-phi_s[j] * ipfactor * phi_n[i]);
              }
            }
            for (size_type j=0; j<lfsu_n.size(); j++) {
              auto temp1 = -(An_F_n*tgradphi_n[j])*omega_n*factor;
              for (size_type i=0; i<lfsu_n.size(); i++) {
                mat_nn.data()[lfsu_n.size()*i+j] += (-omegaup_n * phi_n[j] * normalflux *factor * phi_n[i]);
                mat_nn.data()[lfsu_n.size()*i+j] += (-temp1 * phi_n[i]);
                mat_nn.data()[lfsu_n.size()*i+j] += (-phi_n[j] * factor * theta * omega_n * (An_F_n*tgradphi_n[i]));
                mat_nn.data()[lfsu_n.size()*i+j] += (phi_n[j] * ipfactor * phi_n[i]);
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
        // define types
        using RF = typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSV::Traits::SizeType;

        // dimensions
        const int dim = IG::dimension;
        const int order = std::max(
            lfsu_s.finiteElement().localBasis().order(),
            lfsv_s.finiteElement().localBasis().order()
            );

        // References to the inside cell
        const auto& cell_inside = ig.inside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // Get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // evaluate permeability tensors
        auto ref_el_inside = referenceElement(geo_inside);
        auto local_inside = ref_el_inside.position(0,0);
        auto A_s = param.A(cell_inside,local_inside);

        // face diameter for anisotropic meshes taken from Paul Houston et al.
        // this formula ensures coercivity of the bilinear form
        auto h_F = geo_inside.volume()/geo.volume();

        // compute weights
        auto n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        RF harmonic_average;
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          harmonic_average = An_F_s*n_F;
        else
          harmonic_average = 1.0;

        // get polynomial degree
        auto order_s = lfsu_s.finiteElement().localBasis().order();
        auto degree = order_s;

        // penalty factor
        auto penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // Initialize vectors outside for loop
        std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
        std::vector<Dune::FieldVector<RF,dim> > tgradpsi_s(lfsv_s.size());
        Dune::FieldVector<RF,dim> gradu_s(0.0);

        // Transformation matrix
        typename IG::Entity::Geometry::JacobianInverseTransposed jac;

        // Weight of time stepping
        auto t_weight = r_s.weight();

        // loop over quadrature points
        auto intorder = intorderadd+quadrature_factor*order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            auto bctype = param.bctype(ig.intersection(),ip.position());

            if (bctype == ConvectionDiffusionBoundaryConditions::None)
              continue;

            // position of quadrature point in local coordinates of elements
            auto iplocal_s = geo_in_inside.global(ip.position());

            // local normal
            auto n_F_local = ig.unitOuterNormal(ip.position());

            // evaluate basis functions
            auto& phi_s = cache[order_s].evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());
            auto& psi_s = cache[order_s].evaluateFunction(iplocal_s,lfsv_s.finiteElement().localBasis());

            // integration factor
            RF factor = ip.weight() * geo.integrationElement(ip.position()) * t_weight;

            if (bctype == ConvectionDiffusionBoundaryConditions::Neumann)
              {
                // evaluate flux boundary condition
                auto j = param.j(ig.intersection(),ip.position());

                // integrate
                for (size_type i=0; i<lfsv_s.size(); i++)
                  r_s.accumulate(lfsv_s,i,j * psi_s[i] * factor);

                continue;
              }

            // evaluate u
            RF u_s=0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              u_s += x_s[i]*phi_s[i];

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            auto b = param.b(cell_inside,iplocal_s);
            auto normalflux = b*n_F_local;

            if (bctype == ConvectionDiffusionBoundaryConditions::Outflow)
              {
                if (normalflux<-1e-30)
                  DUNE_THROW(Dune::Exception,
                    "Outflow boundary condition on inflow! [b("
                    << geo.global(ip.position()) << ") = "
                    << b << ")");

                // convection term
                auto term1 = u_s * normalflux *factor;
                for (size_type i=0; i<lfsv_s.size(); i++)
                  r_s.data()[i] += (term1 * psi_s[i]);

                // evaluate flux boundary condition
                auto o = param.o(ig.intersection(),ip.position());

                // integrate
                for (size_type i=0; i<lfsv_s.size(); i++)
                  r_s.data()[i] += (o * psi_s[i] * factor);

                continue;
              }

            // evaluate gradient of basis functions
            assert (bctype == ConvectionDiffusionBoundaryConditions::Dirichlet);
            auto& gradphi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
            auto& gradpsi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsv_s.finiteElement().localBasis());

            // transform gradients of shape functions to real element
            jac = geo_inside.jacobianInverseTransposed(iplocal_s);
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
            for (size_type i=0; i<lfsv_s.size(); i++) jac.mv(gradpsi_s[i][0],tgradpsi_s[i]);

            // compute gradient of u
            gradu_s = 0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              gradu_s.axpy(x_s[i],tgradphi_s[i]);

            // evaluate Dirichlet boundary condition
            auto g = param.g(cell_inside,iplocal_s);

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
            auto term1 = (omegaup_s*u_s + omegaup_n*g) * normalflux *factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.data()[i] += (term1 * psi_s[i]);

            // diffusion term
            auto term2 =  (An_F_s*gradu_s) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.data()[i] += (-term2 * psi_s[i]);

            // (non-)symmetric IP term
            auto term3 = (u_s-g) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.data()[i] += (term3 * theta * (An_F_s*tgradpsi_s[i]));

            // standard IP term
            auto term4 = penalty_factor * (u_s-g) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.data()[i] += (term4 * psi_s[i]);
          }
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              M& mat_ss) const
      {
        // define types
        using RF = typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSV::Traits::SizeType;

        // dimensions
        const int dim = IG::dimension;
        const int order = std::max(
            lfsu_s.finiteElement().localBasis().order(),
            lfsv_s.finiteElement().localBasis().order()
            );

        // References to the inside cell
        const auto& cell_inside = ig.inside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // Get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // evaluate permeability tensors
        auto ref_el_inside = referenceElement(geo_inside);
        auto local_inside = ref_el_inside.position(0,0);
        auto A_s = param.A(cell_inside,local_inside);

        // face diameter for anisotropic meshes taken from Paul Houston et al.
        // this formula ensures coercivity of the bilinear form
        auto h_F = geo_inside.volume()/geo.volume();

        // compute weights
        auto n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        RF harmonic_average;
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          harmonic_average = An_F_s*n_F;
        else
          harmonic_average = 1.0;

        // get polynomial degree
        auto order_s = lfsu_s.finiteElement().localBasis().order();
        auto degree = order_s;

        // penalty factor
        auto penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // Initialize vectors outside for loop
        std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());

        // Transformation matrix
        typename IG::Entity::Geometry::JacobianInverseTransposed jac;

        // Weight of time stepping
        auto t_weight = mat_ss.weight();

        // loop over quadrature points
        auto intorder = intorderadd+quadrature_factor*order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            auto bctype = param.bctype(ig.intersection(),ip.position());

            if (bctype == ConvectionDiffusionBoundaryConditions::None ||
                bctype == ConvectionDiffusionBoundaryConditions::Neumann)
              continue;

            // position of quadrature point in local coordinates of elements
            auto iplocal_s = geo_in_inside.global(ip.position());

            // local normal
            auto n_F_local = ig.unitOuterNormal(ip.position());

            // evaluate basis functions
            auto& phi_s = cache[order_s].evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());

            // integration factor
            auto factor = ip.weight() * geo.integrationElement(ip.position()) * t_weight;

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            auto b = param.b(cell_inside,iplocal_s);
            auto normalflux = b*n_F_local;

            if (bctype == ConvectionDiffusionBoundaryConditions::Outflow)
              {
                if (normalflux<-1e-30)
                  DUNE_THROW(Dune::Exception,
                    "Outflow boundary condition on inflow! [b("
                    << geo.global(ip.position()) << ") = "
                    << b << ")" << n_F_local << " " << normalflux);

                // convection term
                for (size_type j=0; j<lfsu_s.size(); j++)
                  for (size_type i=0; i<lfsu_s.size(); i++)
                    mat_ss.data()[lfsu_s.size()*i+j] += (phi_s[j] * normalflux * factor * phi_s[i]);

                continue;
              }

            // evaluate gradient of basis functions
            auto& gradphi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());

            // transform gradients of shape functions to real element
            jac = geo_inside.jacobianInverseTransposed(iplocal_s);
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
                mat_ss.data()[lfsu_s.size()*i+j] += (omegaup_s * phi_s[j] * normalflux * factor * phi_s[i]);

            // diffusion term
            for (size_type j=0; j<lfsu_s.size(); j++)
              for (size_type i=0; i<lfsu_s.size(); i++)
                mat_ss.data()[lfsu_s.size()*i+j] += (-(An_F_s*tgradphi_s[j]) * factor * phi_s[i]);

            // (non-)symmetric IP term
            for (size_type j=0; j<lfsu_s.size(); j++)
              for (size_type i=0; i<lfsu_s.size(); i++)
                mat_ss.data()[lfsu_s.size()*i+j] += (phi_s[j] * factor * theta * (An_F_s*tgradphi_s[i]));

            // standard IP term
            for (size_type j=0; j<lfsu_s.size(); j++)
              for (size_type i=0; i<lfsu_s.size(); i++)
                mat_ss.data()[lfsu_s.size()*i+j] += (penalty_factor * phi_s[j] * phi_s[i] * factor);
          }
      }

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // define types
        using size_type = typename LFSV::Traits::SizeType;

        // Get cell
        const auto& cell = eg.entity();

        // get geometries
        auto geo = eg.geometry();

        // Weight of time stepping
        auto t_weight = r.weight();

        // loop over quadrature points
        auto order = lfsv.finiteElement().localBasis().order();
        auto intorder = intorderadd + 2 * order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate shape functions
            auto& phi = cache[order].evaluateFunction(ip.position(),lfsv.finiteElement().localBasis());

            // evaluate right hand side parameter function
            auto f = param.f(cell,ip.position());

            // integrate f
            auto factor = ip.weight() * geo.integrationElement(ip.position()) * t_weight;
            for (size_type i=0; i<lfsv.size(); i++)
              r.data()[i] += (-f*phi[i]*factor);
          }
      }

      //! set time in parameter class
      void setTime (Real t)
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

      using LocalBasisType = typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType;
      using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;

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
        using DF = typename GEO::ctype;
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
#endif // DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONFASTDG_HH
