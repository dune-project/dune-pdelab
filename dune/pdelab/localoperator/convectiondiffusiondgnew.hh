// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONDGNEW_HH
#define DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONDGNEW_HH

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

#include"convectiondiffusionparameter.hh"

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
    template<typename Problem>
    class ConvectionDiffusionDGNew :
      //public Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionDG<T,FiniteElementMap> >,
      public Dune::PDELab::FullSkeletonPattern,
      public Dune::PDELab::FullVolumePattern
    //public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>
    {
      //enum { dim = T::Traits::GridViewType::dimension };

      using Real = typename Problem::Real;
      using BCType = typename ConvectionDiffusionBoundaryConditions::Type;

    public:

      ConvectionDiffusionDGNew (Problem& param_,
                             ConvectionDiffusionDGMethod::Type method_=ConvectionDiffusionDGMethod::NIPG,
                             ConvectionDiffusionDGWeights::Type weights_=ConvectionDiffusionDGWeights::weightsOff,
                             Real alpha_=0.0,
                             int intorderadd_=0
                             )
        :
        param(param_), method(method_), weights(weights_),
        alpha(alpha_), intorderadd(intorderadd_), quadrature_factor(2)
      {
        theta = 1.0;
        if (method==ConvectionDiffusionDGMethod::SIPG) theta = -1.0;
        if (method==ConvectionDiffusionDGMethod::IIPG) theta = 0.0;
      }


      template<typename Context>
      void volumePattern(Context& ctx) const
      {
        for (auto i : ctx.test().space())
          for (auto j : ctx.trial().space())
            ctx.pattern().addLink(ctx.test().space(),i,ctx.trial().space(),j);
      }

      template<typename Context>
      void skeletonPattern(Context& ctx) const
      {
        auto& inside  = ctx.inside();
        auto& outside = ctx.outside();

        auto& pattern_io = ctx.pattern(inside,outside);
        for (auto i : inside.test().space())
          for (auto j : outside.trial().space())
            pattern_io.addLink(inside.test().space(),i,outside.trial().space(),j);

        auto& pattern_oi = ctx.pattern(outside,inside);
        for (auto i : outside.test().space())
          for (auto j : inside.trial().space())
            pattern_oi.addLink(outside.test().space(),i,inside.trial().space(),j);
      }


      template<typename Context>
      void volumeIntegral(Context& ctx) const
      {
        // extract some useful types
        using RF        = LocalOperator::RangeField<Context>;
        using Range     = LocalOperator::Range<Context>;
        using Gradient  = LocalOperator::Gradient<Context>;
        using size_type = std::size_t;

        // short cuts for cell, bases and function spaces
        auto& cell = ctx.cell();
        auto& trial_space = ctx.trial().functionSpace();
        auto& test_space = ctx.test().functionSpace();
        auto& trial_basis = ctx.trial().basis();
        auto& test_basis = ctx.test().basis();

        // determine integration order
        auto intorder = intorderadd + quadrature_factor * std::max(trial_basis.order(),test_basis.order());

        auto A = cell.A(cell.centroid());

        // variables for quantities used inside the quadrature loop
        Range u;
        Gradient gradu, Agradu;
        typename Context::Traits::Velocity b;
        RF c, f;

        // loop over quadrature points
        for (auto ip : ctx.domain().quadratureRule(intorder))
          {

            // evaluate test basis functions
            auto psi = test_basis(ip);

            if (not ctx.skipVariablePart())
              {
                // update all variables dependent on A if A is not cell-wise constant
                if (!cell.permeabilityIsConstantPerCell())
                  {
                    A = cell.A(ip);
                  }

                // evaluate trial basis functions
                auto phi = trial_basis(ip);

                // evaluate u
                u = 0.0;
                for (size_type i=0 ; i < trial_space.size() ; ++i)
                  u += cell.argument()(trial_space,i) * phi[i];

                // evaluate gradient of basis functions
                auto gradphi = trial_basis.gradients(ip);

                // compute gradient of u
                gradu = 0.0;
                for (size_type i = 0 ; i < trial_space.size() ; ++i)
                  gradu.axpy(cell.argument()(trial_space,i),gradphi[i]);

                // compute A * gradient of u
                A.mv(gradu,Agradu);

                // evaluate velocity field
                b = cell.b(ip);

                // evaluate reaction term
                c = cell.c(ip);
              }

            if (not ctx.skipConstantPart())
              {
                // evaluate right hand side parameter function
                f = cell.f(ip);
              }

            for (auto [dof, i] : ctx.residual(test_space))
              {
                RF r = 0.0;
                if (not ctx.skipVariablePart())
                  {
                    auto gradpsi = test_basis.gradients(ip);

                    // integrate (A grad u - bu)*grad psi_i + a*u*psi_i
                    r += Agradu*gradpsi[i] - u*(b*gradpsi[i]) + c*u*psi[i];
                  }
                if (not ctx.skipConstantPart())
                  {
                    // integrate -f*psi_i
                    r -= f * psi[i];
                  }
                dof += r*ip.weight();
              }
          }
      }


      template<typename Context>
      void volumeJacobian(Context& ctx) const
      {
        // extract some useful types
        using size_type = std::size_t;

        // short cuts for cell, bases and function spaces
        auto& cell = ctx.cell();
        auto domain = ctx.domain();
        auto& trial_space = ctx.trial().functionSpace();
        auto& test_space = ctx.test().functionSpace();
        auto& trial_basis = ctx.trial().basis();
        auto& test_basis = ctx.test().basis();

        // determine integration order
        auto intorder = intorderadd + quadrature_factor * std::max(trial_basis.order(),test_basis.order());

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        auto A = cell.A(cell.centroid());

        // allocate storage for A * gradphi
        cell.Agradphi.resize(trial_space.size());

        // loop over quadrature points
        for (auto ip : domain.quadratureRule(intorder))
          {

            // update all variables dependent on A if A is not cell-wise constant
            if (not cell.permeabilityIsConstantPerCell())
            {
              A = cell.A(ip);
            }

            // evaluate basis functions and gradients
            auto phi = trial_basis.values(ip);
            auto psi = test_basis.values(ip);

            auto gradphi = trial_basis.gradients(ip);
            auto gradpsi = test_basis.gradients(ip);

            for (size_type i = 0 ; i < trial_basis.size() ; ++i)
              A.mv(gradphi[i],cell.Agradphi[i]);

            // evaluate velocity field
            auto b = cell.b(ip);

            // evaluate reaction term
            auto c = cell.c(ip);

            // integrate (A grad u - bu)*grad phi_i + a*u*phi_i
            for (auto [trialFunctions,i] : cell.jacobian(test_space,trial_space))
              for (auto [dof,j] : trialFunctions)
                dof += (cell.Agradphi[j]*gradpsi[i] - phi[j]*(b*gradpsi[i]) + c*phi[j]*psi[i]) * ip.weight();

          }
      }

      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename Context>
      void skeletonIntegral(Context& ctx) const
      {
        if (ctx.skipVariablePart())
          return;

        // extract some useful types
        using RF        = LocalOperator::RangeField<Context>;
        using Range     = LocalOperator::Range<Context>;
        using Gradient  = LocalOperator::Gradient<Context>;
        using size_type = std::size_t;

        auto& inside  = ctx.inside();
        auto& outside = ctx.outside();
        auto  domain  = ctx.domain();

        // determine integration order
        auto order = std::max(
          std::max(inside.trial().basis().order(),outside.test().basis().order()),
          std::max(inside.test().basis().order(),outside.test().basis().order())
          );
        auto intorder = intorderadd + quadrature_factor * order;

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        auto A_s = inside.A(inside.centroid());
        auto A_n = outside.A(outside.centroid());

        // variables for quantities used inside the quadrature loop
        Range u_s, u_n;
        Gradient gradu_s, gradu_n;
        typename Context::Traits::Velocity b;

        // face diameter for anisotropic meshes taken from Paul Houston et al.
        // this formula ensures coercivity of the bilinear form
        auto h_F = std::min(inside.volume(),outside.volume())/domain.volume();

        // tensor times normal
        auto n_F = domain.centerUnitOuterNormal();
        auto An_F_s = n_F;
        auto An_F_n = n_F;
        A_s.mv(n_F,An_F_s);
        A_n.mv(n_F,An_F_n);

        // compute weights
        RF omega_s;
        RF omega_n;
        RF harmonic_average(0.0);
        if (weights == ConvectionDiffusionDGWeights::weightsOn)
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
        auto degree = std::max(inside.trial().basis().order(), outside.trial().basis().order());

        // penalty factor
        auto penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+Context::dimWorld-1);

        // loop over quadrature points
        for (const auto& ip : domain.quadratureRule(intorder))
          {

            // exact normal
            auto n_F_local = ip.unitOuterNormal();

            // update all variables dependent on A if A is not cell-wise constant
            if (not inside.permeabilityIsConstantPerCell())
            {
              A_s = inside.A(ip);
              A_n = outside.A(ip);
              A_s.mv(n_F_local,An_F_s);
              A_n.mv(n_F_local,An_F_n);
              if (weights==ConvectionDiffusionDGWeights::weightsOn)
                {
                  RF delta_s = (An_F_s*n_F);
                  RF delta_n = (An_F_n*n_F);
                  omega_s = delta_n/(delta_s+delta_n+1e-20);
                  omega_n = delta_s/(delta_s+delta_n+1e-20);
                  harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
                  penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+Context::dimWorld-1);
                }
            }

            // evaluate basis functions
            auto phi_s = inside.trial().basis()(ip);
            auto phi_n = outside.trial().basis()(ip);
            auto psi_s = inside.test().basis()(ip);
            auto psi_n = outside.test().basis()(ip);

            // evaluate u
            u_s=0.0;
            for (size_type i = 0 ; i < inside.trial().basis().size() ; ++i)
              u_s += inside.argument()(inside.trial().space(),i) * phi_s[i];
            u_n=0.0;
            for (size_type i = 0 ; i < outside.trial().basis().size() ; ++i)
              u_n += outside.argument()(outside.trial().space(),i) * phi_n[i];

            // evaluate gradient of basis functions
            auto gradphi_s = inside.trial().basis().gradients(ip);
            auto gradphi_n = outside.trial().basis().gradients(ip);
            auto gradpsi_s = inside.test().basis().gradients(ip);
            auto gradpsi_n = outside.test().basis().gradients(ip);

            // compute gradient of u
            gradu_s = 0.0;
            for (size_type i = 0 ; i < inside.trial().space().size() ; ++i)
              gradu_s.axpy(inside.argument()(inside.trial().space(),i),gradphi_s[i]);
            gradu_n = 0.0;
            for (size_type i = 0 ; i < outside.trial().space().size() ; ++i)
              gradu_n.axpy(outside.argument()(outside.trial().space(),i),gradphi_n[i]);

            // evaluate velocity field and upwinding, assume H(div) velocity field => may choose any side
            b = inside.b(ip);
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

            auto inside_residual  = inside.residual(inside.test().space());
            auto outside_residual = outside.residual(outside.test().space());

            // convection term
            auto term1 = (omegaup_s*u_s + omegaup_n*u_n) * normalflux * ip.weight();
            for (auto [dof, i] : inside_residual)
              dof += term1 * psi_s[i];
            for (auto [dof, i] : outside_residual)
              dof -= term1 * psi_n[i];

            // diffusion term
            auto term2 =  -(omega_s*(An_F_s*gradu_s) + omega_n*(An_F_n*gradu_n)) * ip.weight();
            for (size_type i = 0 ; i < inside_residual.size() ; ++i)
              inside_residual.accumulate(i,term2 * psi_s[i]);
            for (size_type i = 0 ; i < outside_residual.size() ; ++i)
              outside_residual.accumulate(i,-term2 * psi_n[i]);

            // (non-)symmetric IP term
            auto term3 = (u_s-u_n) * ip.weight();
            for (size_type i = 0 ; i < inside.trial().space().size() ; ++i)
              inside_residual.accumulate(inside.trial().space(),i,term3 * theta * omega_s * (An_F_s * gradpsi_s[i]));
            for (size_type i = 0 ; i < outside.trial().space().size() ; ++i)
              outside_residual.accumulate(outside.trial().space(),i,term3 * theta * omega_n * (An_F_n * gradpsi_n[i]));

            // standard IP term integral
            auto term4 = penalty_factor * (u_s-u_n) * ip.weight();
            for (auto [dof, i] : inside_residual)
              dof += term4 * psi_s[i];
            for (auto [dof, i] : outside_residual)
              dof -= term4 * psi_n[i];
          }
      }


      template<typename Context>
      void skeletonJacobian(Context& ctx) const
      {
        // extract some useful types
        using RF        = LocalOperator::RangeField<Context>;
        using size_type = std::size_t;

        auto& inside  = ctx.inside();
        auto& outside = ctx.outside();
        auto  domain  = ctx.domain();

        // determine integration order
        auto order = std::max(
          std::max(inside.trial().basis().order(),outside.test().basis().order()),
          std::max(inside.test().basis().order(),outside.test().basis().order())
          );
        auto intorder = intorderadd + quadrature_factor * order;

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        auto A_s = inside.A(inside.centroid());
        auto A_n = outside.A(outside.centroid());

        // face diameter for anisotropic meshes taken from Paul Houston et al.
        // this formula ensures coercivity of the bilinear form
        auto h_F = std::min(inside.volume(),outside.volume())/domain.volume();

        // tensor times normal
        auto n_F = domain.centerUnitOuterNormal();
        auto An_F_s = n_F;
        auto An_F_n = n_F;
        A_s.mv(n_F,An_F_s);
        A_n.mv(n_F,An_F_n);

        // compute weights
        RF omega_s;
        RF omega_n;
        RF harmonic_average(0.0);
        if (weights == ConvectionDiffusionDGWeights::weightsOn)
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
        auto degree = std::max(inside.trial().basis().order(), outside.trial().basis().order());

        // penalty factor
        auto penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+Context::dimWorld-1);

        // loop over quadrature points
        for (auto ip : domain.quadratureRule(intorder))
          {

            // exact normal
            auto n_F_local = ip.unitOuterNormal();

            // update all variables dependent on A if A is not cell-wise constant
            if (not inside.permeabilityIsConstantPerCell())
            {
              A_s = inside.A(ip);
              A_n = outside.A(ip);
              A_s.mv(n_F_local,An_F_s);
              A_n.mv(n_F_local,An_F_n);
              if (weights==ConvectionDiffusionDGWeights::weightsOn)
                {
                  RF delta_s = (An_F_s*n_F);
                  RF delta_n = (An_F_n*n_F);
                  omega_s = delta_n/(delta_s+delta_n+1e-20);
                  omega_n = delta_s/(delta_s+delta_n+1e-20);
                  harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
                  penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+Context::dimWorld-1);
                }
            }

            // evaluate basis functions
            auto phi_s = inside.trial().basis()(ip);
            auto phi_n = outside.trial().basis()(ip);
            auto psi_s = inside.test().basis()(ip);
            auto psi_n = outside.test().basis()(ip);

            // evaluate gradients of basis functions
            auto gradphi_s = inside.trial().basis().gradients(ip);
            auto gradphi_n = outside.trial().basis().gradients(ip);
            auto gradpsi_s = inside.test().basis().gradients(ip);
            auto gradpsi_n = outside.test().basis().gradients(ip);

            // evaluate velocity field and upwinding, assume H(div) velocity field => may choose any side
            auto b = inside.b(ip);
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
            auto ipfactor = penalty_factor * ip.weight();

            // do all terms in the order: I convection, II diffusion, III consistency, IV ip
            for (auto [trialFunctions,i] : ctx.jacobian(inside.test().space(),inside.trial().space()))
              for (auto [dof,j] : trialFunctions)
                {
                  dof += omegaup_s * phi_s[j] * normalflux * ip.weight() * psi_s[i];
                  dof -= (An_F_s*gradphi_s[j])*omega_s*ip.weight() * psi_s[i];
                  dof += phi_s[j] * ip.weight() * theta * omega_s * (An_F_s*gradpsi_s[i]);
                  dof += phi_s[j] * ipfactor * psi_s[i];
                }
            for (auto [trialFunctions,i] : ctx.jacobian(inside.test().space(),outside.trial().space()))
              for (auto [dof,j] : trialFunctions)
                {
                  dof += omegaup_n * phi_n[j] * normalflux * ip.weight() * psi_s[i];
                  dof -= (An_F_n*gradphi_n[j])*omega_n*ip.weight() * psi_s[i];
                  dof -= phi_n[j] * ip.weight() * theta * omega_s * (An_F_s*gradpsi_s[i]);
                  dof -= phi_n[j] * ipfactor * psi_s[i];
                }
            for (auto [trialFunctions,i] : ctx.jacobian(outside.test().space(),inside.trial().space()))
              for (auto [dof,j] : trialFunctions)
                {
                  dof -= omegaup_s * phi_s[j] * normalflux * ip.weight() * psi_n[i];
                  dof += (An_F_s*gradphi_s[j])*omega_s*ip.weight() * psi_n[i];
                  dof += phi_s[j] * ip.weight() * theta * omega_n * (An_F_n*gradpsi_n[i]);
                  dof -= phi_s[j] * ipfactor * psi_n[i];
                }
            for (auto [trialFunctions,i] : ctx.jacobian(outside.test().space(),outside.trial().space()))
              for (auto [dof,j] : trialFunctions)
                {
                  dof -= omegaup_n * phi_n[j] * normalflux * ip.weight() * psi_n[i];
                  dof += (An_F_n*gradphi_n[j])*omega_n*ip.weight() * psi_n[i];
                  dof -= phi_n[j] * ip.weight() * theta * omega_n * (An_F_n*gradpsi_n[i]);
                  dof += phi_n[j] * ipfactor * psi_n[i];
                }
          }
      }

      // boundary integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename Context>
      void boundaryIntegral(Context& ctx) const
      {
        // extract some useful types
        using RF        = LocalOperator::RangeField<Context>;
        using Gradient  = LocalOperator::Gradient<Context>;
        using size_type = std::size_t;

        auto domain = ctx.domain();
        auto& cell   = ctx.cell();

        // dimensions
        auto order    = std::max(ctx.trial().basis().order(),ctx.test().basis().order());
        auto intorder = intorderadd + quadrature_factor * order;

        auto A_s = cell.A(cell.centroid());

        // face diameter for anisotropic meshes taken from Paul Houston et al.
        // this formula ensures coercivity of the bilinear form
        auto h_F = cell.volume()/domain.volume();

        // compute weights
        auto n_F = domain.centerUnitOuterNormal();
        auto An_F_s = n_F;
        A_s.mv(n_F,An_F_s);

        RF harmonic_average;
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          harmonic_average = An_F_s*n_F;
        else
          harmonic_average = 1.0;

        // get polynomial degree
        auto degree = ctx.trial().basis().order();

        // penalty factor
        auto penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+Context::dimWorld-1);

        // loop over quadrature points
        for (auto ip : domain.quadratureRule(intorder))
          {

            auto bctype = ctx.bctype(ip);

            if (bctype == ConvectionDiffusionBoundaryConditions::None)
              continue;

            // local normal
            auto n_F_local = ip.unitOuterNormal();

            // update all variables dependent on A if A is not cell-wise constant
            if (not cell.permeabilityIsConstantPerCell())
            {
              A_s = cell.A(ip);
              A_s.mv(n_F_local,An_F_s);
              if (weights==ConvectionDiffusionDGWeights::weightsOn)
                {
                  harmonic_average = An_F_s*n_F;
                  penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+Context::dimWorld-1);
                }
            }

            // evaluate basis functions
            auto phi_s = cell.trial().basis()(ip);
            auto psi_s = cell.test().basis()(ip);

            if (bctype == ConvectionDiffusionBoundaryConditions::Neumann)
              {
                // evaluate flux boundary condition
                auto j = ctx.j(ip);

                // integrate
                for (size_type i = 0 ; i < ctx.test().space().size() ; ++i)
                  cell.residual().accumulate(ctx.test().space(),i,j * psi_s[i] * ip.weight());

                continue;
              }

            // evaluate u
            RF u_s = 0.0;
            for (size_type i = 0 ; i < ctx.trial().space().size() ; ++i)
              u_s += cell.argument()(ctx.trial().space(),i) * phi_s[i];

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            auto b = cell.b(ip);
            auto normalflux = b * n_F_local;

            if (bctype == ConvectionDiffusionBoundaryConditions::Outflow)
              {
                if (normalflux < -1e-30)
                  DUNE_THROW(Dune::Exception,
                    "Outflow boundary condition on inflow! [b("
                    << ip.global() << ") = "
                    << b << ")");

                // convection term
                auto term1 = u_s * normalflux * ip.weight();
                for (size_type i = 0 ; i < ctx.test().space().size() ; ++i)
                  cell.residual().accumulate(ctx.test().space(),i,term1 * psi_s[i]);

                // evaluate flux boundary condition
                auto o = ctx.o(ip) * ip.weight();

                // integrate
                for (size_type i = 0 ; i < ctx.test().space().size() ; ++i)
                  cell.residual().accumulate(ctx.test().space(),i,o * psi_s[i]);

                continue;
              }

            // evaluate gradient of basis functions
            assert (bctype == ConvectionDiffusionBoundaryConditions::Dirichlet);
            auto gradphi_s = cell.trial().basis().gradients(ip);
            auto gradpsi_s = cell.test().basis().gradients(ip);

            // compute gradient of u
            Gradient gradu_s(0.0);
            for (size_type i = 0 ; i < cell.trial().space().size(); ++i)
              gradu_s.axpy(cell.argument()(cell.trial().space(),i),gradphi_s[i]);

            // evaluate Dirichlet boundary condition
            auto g = cell.g(ip);

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
            auto term1 = (omegaup_s*u_s + omegaup_n*g) * normalflux * ip.weight();
            for (size_type i = 0 ; i < ctx.test().space().size() ; ++i)
              cell.residual().accumulate(ctx.test().space(),i,term1 * psi_s[i]);

            // diffusion term
            auto term2 = (An_F_s*gradu_s) * ip.weight();
            for (size_type i = 0 ; i < ctx.test().space().size() ; ++i)
              cell.residual().accumulate(ctx.test().space(),i,-term2 * psi_s[i]);

            // (non-)symmetric IP term
            auto term3 = (u_s-g) * ip.weight();
            for (size_type i = 0 ; i < ctx.test().space().size() ; ++i)
              cell.residual().accumulate(ctx.test().space(),i,term3 * theta * (An_F_s*gradpsi_s[i]));

            // standard IP term
            auto term4 = penalty_factor * (u_s-g) * ip.weight();
            for (size_type i = 0 ; i < ctx.test().space().size() ; ++i)
              cell.residual().accumulate(ctx.test().space(),i,term4 * psi_s[i]);
          }
      }


      template<typename Context>
      void boundaryJacobian(Context& ctx) const
      {
        // extract some useful types
        using RF = LocalOperator::RangeField<Context>;

        auto domain = ctx.domain();
        auto& cell  = ctx.cell();

        // dimensions
        auto order    = std::max(ctx.trial().basis().order(),ctx.test().basis().order());
        auto intorder = intorderadd + quadrature_factor * order;

        auto A_s = cell.A(cell.centroid());

        // face diameter for anisotropic meshes taken from Paul Houston et al.
        // this formula ensures coercivity of the bilinear form
        auto h_F = cell.volume()/domain.volume();

        // compute weights
        auto n_F = domain.centerUnitOuterNormal();
        auto An_F_s = n_F;
        A_s.mv(n_F,An_F_s);

        RF harmonic_average;
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          harmonic_average = An_F_s*n_F;
        else
          harmonic_average = 1.0;

        // get polynomial degree
        auto degree = ctx.trial().basis().order();

        // penalty factor
        auto penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+Context::dimWorld-1);

        // loop over quadrature points
        for (auto ip : domain.quadratureRule(intorder))
          {
            auto bctype = ctx.bctype(ip);

            if (bctype == ConvectionDiffusionBoundaryConditions::None ||
                bctype == ConvectionDiffusionBoundaryConditions::Neumann)
              continue;

            // local normal
            auto n_F_local = domain.unitOuterNormal(ip);

            // update all variables dependent on A if A is not cell-wise constant
            if (not cell.permeabilityIsConstantPerCell())
            {
              A_s = cell.A(ip);
              A_s.mv(n_F_local,An_F_s);
              if (weights==ConvectionDiffusionDGWeights::weightsOn)
                {
                  harmonic_average = An_F_s*n_F;
                  penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+Context::dimWorld-1);
                }
            }

            // evaluate basis functions
            auto phi_s = cell.trial().basis()(ip);
            auto psi_s = cell.test().basis()(ip);

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            auto b = cell.b(ip);
            auto normalflux = b*n_F_local;

            auto mat_ss = ctx.jacobian(cell.test().space(),cell.trial().space());

            if (bctype == ConvectionDiffusionBoundaryConditions::Outflow)
              {
                if (normalflux<-1e-30)
                  DUNE_THROW(Dune::Exception,
                    "Outflow boundary condition on inflow! [b("
                    << ip.global() << ") = "
                    << b << ")" << n_F_local << " " << normalflux);

                for (auto [trialFunctions,i] : mat_ss)
                  for (auto [dof,j] : trialFunctions)
                    dof += phi_s[j] * normalflux * ip.weight() * psi_s[i];

                continue;
              }

            // evaluate gradient of basis functions
            auto gradphi_s = cell.trial().basis().gradients(ip);
            auto gradpsi_s = cell.test().basis().gradients(ip);

            // upwind
            RF omegaup_s = normalflux>=0.0 ? 1.0 : 0.0;

            for (auto [trialFunctions,i] : mat_ss)
              for (auto [dof,j] : trialFunctions)
                {
                  // convection term
                  dof += omegaup_s * phi_s[j] * normalflux * ip.weight() * psi_s[i];

                  // diffusion term
                  dof -= (An_F_s*gradphi_s[j]) * ip.weight() * psi_s[i];

                  // (non-)symmetric IP term
                  dof += phi_s[j] * ip.weight() * theta * (An_F_s*gradpsi_s[i]);

                  // standard IP term
                  dof += penalty_factor * phi_s[j] * psi_s[i] * ip.weight();
                }
          }
      }

      //! set time in parameter class
      /*
      void setTime (Real t)
      {
        Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>::setTime(t);
        param.setTime(t);
      }
      */

      template<typename Ctx>
      struct CellCache
        : public Ctx
      {

        CellCache(Ctx&& ctx)
          : Ctx(std::move(ctx))
        {}

        std::vector<LocalOperator::Gradient<Ctx>> Agradphi;
      };

      template<typename Ctx>
      using CellContext = typename Problem::template ProblemContext<CellCache<Ctx>>;

      template<typename Ctx>
      using Context = typename Problem::template ProblemContext<Ctx>;

      const Problem& problem() const
      {
        return param;
      }

    private:
      const Problem& param;  // two phase parameter class
      ConvectionDiffusionDGMethod::Type method;
      ConvectionDiffusionDGWeights::Type weights;
      Real alpha, beta;
      int intorderadd;
      int quadrature_factor;
      Real theta;

      //using LocalBasisType = typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType;
      //using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;

      // In theory it is possible that one and the same local operator is
      // called first with a finite element of one type and later with a
      // finite element of another type.  Since finite elements of different
      // type will usually produce different results for the same local
      // coordinate they cannot share a cache.  Here we use a vector of caches
      // to allow for different orders of the shape functions, which should be
      // enough to support p-adaptivity.  (Another likely candidate would be
      // differing geometry types, i.e. hybrid meshes.)

      // std::vector<Cache> cache;

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
#endif // DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONDG_HH
