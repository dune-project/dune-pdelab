// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_TEST_NONLINEARCONVECTIONDIFFUSIONDG_HH
#define DUNE_PDELAB_TEST_NONLINEARCONVECTIONDIFFUSIONDG_HH

#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>

namespace Dune {
  namespace PDELab {

    /** A local operator that adds a nonlinearity q(u) to CovectionDiffusionDG
     *
     * The main purpose of this operator is for testing nonlinear problem
     * related code. For that reason it was not put into the localoperator
     * folder but directly into the test folder.
     *
     * \f{align*}{
     *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u + q(u) &=& f \mbox{ in } \Omega,  \\
     *                                                     u &=& g \mbox{ on } \partial\Omega_D \\
     *                       (b(x) u - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
     *                               -(A(x)\nabla u) \cdot n &=& o \mbox{ on } \partial\Omega_O
     * \f}
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T, typename FiniteElementMap>
    class NonlinearConvectionDiffusionDG :
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

      enum { isLinear = false };

      /** \brief constructor: pass parameter object and define DG-method
       * \param[in] param_   Reference to parameter object.
       * \param[in] method_  Interior penalty Galerkin method. Default is skew-symmetric.
       * \param[in] weights_ Weighted averages for diffusion tensor. Default is no weighting.
       * \param[in] alpha_   Penalization constant. Default is zero.
       *
       * Collecting the input parameters above, the default is the OBB-method.
       */
      NonlinearConvectionDiffusionDG (T& param_,
                                      ConvectionDiffusionDGMethod::Type method_=ConvectionDiffusionDGMethod::NIPG,
                                      ConvectionDiffusionDGWeights::Type weights_=ConvectionDiffusionDGWeights::weightsOff,
                                      Real alpha_=0.0,
                                      int intorderadd_=0
                                      )
        : param(param_)
        , method(method_)
        , weights(weights_)
        , alpha(alpha_)
        , intorderadd(intorderadd_)
        , quadrature_factor(2)
        , cache(20)
        , linearLocalOperator(param_, method_, weights_, alpha_, intorderadd_)
      {
        theta = 1.0;
        if (method==ConvectionDiffusionDGMethod::SIPG) theta = -1.0;
        if (method==ConvectionDiffusionDGMethod::IIPG) theta = 0.0;
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // Integrate linear part here
        linearLocalOperator.alpha_volume(eg, lfsu, x, lfsv, r);

        // Define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSU::Traits::SizeType;

        // Get geometry
        auto geo = eg.geometry();

        // Loop over quadrature points
        const int order = std::max(lfsu.finiteElement().localBasis().order(),
            lfsv.finiteElement().localBasis().order());
        auto intorder = intorderadd + quadrature_factor * order;
        for (const auto& ip : quadratureRule(geo, intorder))
        {
          // Evaluate basis functions
          auto& phi = cache[order].evaluateFunction(ip.position(), lfsu.finiteElement().localBasis());
          auto& psi = cache[order].evaluateFunction(ip.position(), lfsv.finiteElement().localBasis());

          // Evaluate u
          RF u = 0.0;
          for (size_type i=0; i<lfsu.size(); i++)
            u += x(lfsu,i) * phi[i];

          // Evaluate nonlinearity
          auto q = param.q(u);

          // Add integral q(u)*phi_i
          RF factor = ip.weight() * geo.integrationElement(ip.position());
          for (size_type i=0; i<lfsv.size(); i++)
            r.accumulate(lfsv, i, q * psi[i] * factor);
        }
      }

      // apply jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
      void jacobian_apply_volume (const EG& eg, const LFSU& lfsu, const X& x, const Z& z, const LFSV& lfsv, Y& y) const
      {
        // Integrate linear part here
        linearLocalOperator.jacobian_apply_volume(eg, lfsu, z, lfsv, y);

        // Define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSU::Traits::SizeType;

        // Get geometry
        auto geo = eg.geometry();

        // Loop over quadrature points
        const int order = std::max(lfsu.finiteElement().localBasis().order(),
            lfsv.finiteElement().localBasis().order());
        auto intorder = intorderadd + quadrature_factor * order;
        for (const auto& ip : quadratureRule(geo, intorder))
        {
          // Evaluate basis functions
          auto& phi = cache[order].evaluateFunction(ip.position(), lfsu.finiteElement().localBasis());
          auto& psi = cache[order].evaluateFunction(ip.position(), lfsv.finiteElement().localBasis());

          // Evaluate u
          RF u = 0.0;
          for (size_type i=0; i<lfsu.size(); i++)
            u += x(lfsu,i) * phi[i];

          // Evaluate nonlinearity
          auto qprime = param.qprime(u);

          // Evaluate z
          RF z_eval = 0.0;
          for (size_type i=0; i<lfsu.size(); i++)
            z_eval += z(lfsu,i) * phi[i];

          // Add integral qprime * z * psi_i
          RF factor = ip.weight() * geo.integrationElement(ip.position());
          for (size_type i=0; i<lfsv.size(); i++)
            y.accumulate(lfsv, i, qprime * z_eval * psi[i] * factor);
        }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // Integrate linear part here
        linearLocalOperator.jacobian_volume(eg, lfsu, x, lfsv, mat);

        // Define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSU::Traits::SizeType;

        // Get geometry
        auto geo = eg.geometry();

        // Loop over quadrature points
        const int order = std::max(lfsu.finiteElement().localBasis().order(),
            lfsv.finiteElement().localBasis().order());
        auto intorder = intorderadd + quadrature_factor * order;
        for (const auto& ip : quadratureRule(geo, intorder))
        {
          // Evaluate basis functions
          auto& phi = cache[order].evaluateFunction(ip.position(), lfsu.finiteElement().localBasis());
          auto& psi = cache[order].evaluateFunction(ip.position(), lfsv.finiteElement().localBasis());

          // Evaluate u
          RF u=0.0;
          for (size_type i=0; i<lfsu.size(); i++)
            u += x(lfsu,i)*phi[i];

          // Evaluate nonlinearity
          auto qprime = param.qprime(u);

          // Add integral q(u)*phi_j*psi_i
          RF factor = ip.weight() * geo.integrationElement(ip.position());
          for (size_type j=0; j<lfsu.size(); j++)
            for (size_type i=0; i<lfsu.size(); i++)
              mat.accumulate(lfsv, i, lfsu, j, qprime * phi[j]* psi[i] * factor);
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
        linearLocalOperator.alpha_skeleton(ig, lfsu_s, x_s, lfsv_s, lfsu_n, x_n, lfsv_n, r_s, r_n);
      }

      // apply jacobian of skeleton term
      template<typename IG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
      void jacobian_apply_skeleton (const IG& ig,
                                    const LFSU& lfsu_s, const X& x_s, const Z& z_s, const LFSV& lfsv_s,
                                    const LFSU& lfsu_n, const X& x_n, const Z& z_n, const LFSV& lfsv_n,
                                    Y& y_s, Y& y_n) const
      {
        linearLocalOperator.jacobian_apply_skeleton(ig, lfsu_s, z_s, lfsv_s, lfsu_n, z_n, lfsv_n, y_s, y_n);
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                              M& mat_ss, M& mat_sn,
                              M& mat_ns, M& mat_nn) const
      {
        linearLocalOperator.jacobian_skeleton(ig, lfsu_s, x_s, lfsv_s, lfsu_n, x_n, lfsv_n, mat_ss, mat_sn, mat_ns, mat_nn);
      }

      // boundary integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                         const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                         R& r_s) const
      {
        linearLocalOperator.alpha_boundary(ig, lfsu_s, x_s, lfsv_s, r_s);
      }

      // apply jacobian of boundary term
      template<typename IG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
      void jacobian_apply_boundary (const IG& ig,
                                    const LFSU& lfsu_s, const X& x_s, const Z& z_s, const LFSV& lfsv_s,
                                    Y& y_s) const
      {
        linearLocalOperator.jacobian_apply_boundary(ig, lfsu_s, z_s, lfsv_s, y_s);
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              M& mat_ss) const
      {
        linearLocalOperator.jacobian_boundary(ig, lfsu_s, x_s, lfsv_s, mat_ss);
      }

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        linearLocalOperator.lambda_volume(eg, lfsv, r);
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

      ConvectionDiffusionDG<T, FiniteElementMap> linearLocalOperator;
    };

  } // namespace PDELab
} // namespace Dune
#endif
