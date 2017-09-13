// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_DIFFUSIONMIXED_HH
#define DUNE_PDELAB_LOCALOPERATOR_DIFFUSIONMIXED_HH

#include <cstddef>
#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>

#include<dune/geometry/type.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/quadraturerules.hh>

#include"defaultimp.hh"
#include"pattern.hh"
#include"flags.hh"
#include "convectiondiffusionparameter.hh"

namespace Dune {
  namespace PDELab {

    // a local operator for solving the Poisson equation
    //     div sigma +a_0 u = f         in \Omega,
    //                sigma = -K grad u in \Omega,
    //                    u = g         on \partial\Omega_D
    //      sigma \cdot \nu = j         on \partial\Omega_N
    // with H(div) conforming (mixed) finite elements
    // param.A : diffusion tensor dependent on position
    // param.c : Helmholtz term
    // param.f : grid function type giving f
    // param.bctype : grid function type selecting boundary condition
    // param.g : grid function type giving g
    template<typename PARAM>
    class DiffusionMixed : public NumericalJacobianApplyVolume<DiffusionMixed<PARAM> >,
                           public NumericalJacobianVolume<DiffusionMixed<PARAM> >,
                           public FullVolumePattern,
                           public LocalOperatorDefaultFlags
    {

      using BCType = typename ConvectionDiffusionBoundaryConditions::Type;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doLambdaVolume = true };
      enum { doLambdaBoundary = true };

      DiffusionMixed ( const PARAM& param_,
                       int qorder_v_=2,
                       int qorder_p_=1 )
        : param(param_),
          qorder_v(qorder_v_),
          qorder_p(qorder_p_)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // Define types
        using VelocitySpace = typename LFSU::template Child<0>::Type;
        using PressureSpace = typename LFSU::template Child<1>::Type;
        using DF = typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType;
        using RF = typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using VelocityJacobianType = typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType;
        using VelocityRangeType = typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using PressureRangeType = typename PressureSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;

        // select the two components
        using namespace Indices;
        const auto& velocityspace = child(lfsu,_0);
        const auto& pressurespace = lfsu.template child<1>();

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // References to the cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // evaluate transformation which must be linear
        Dune::FieldVector<DF,dim> pos;
        pos=0.0;
        auto jac = geo.jacobianInverseTransposed(pos);
        jac.invert();
        auto det = geo.integrationElement(pos);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        auto ref_el = referenceElement(geo);
        auto localcenter = ref_el.position(0,0);
        auto tensor = param.A(cell,localcenter);
        tensor.invert(); // need iverse for mixed method

        // Initialize vectors outside for loop
        std::vector<VelocityRangeType> vbasis(velocityspace.size());
        std::vector<VelocityRangeType> vtransformedbasis(velocityspace.size());
        VelocityRangeType sigma;
        VelocityRangeType Kinvsigma;
        std::vector<VelocityJacobianType> vjacbasis(velocityspace.size());
        std::vector<PressureRangeType> pbasis(pressurespace.size());
        std::vector<RF> divergence(velocityspace.size(),0.0);


        // \sigma\cdot v term
        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder_v))
          {
            // evaluate shape functions at ip (this is a Galerkin method)
            velocityspace.finiteElement().localBasis().evaluateFunction(ip.position(),vbasis);

            // transform basis vectors
            for (std::size_t i=0; i<velocityspace.size(); i++)
              {
                vtransformedbasis[i] = 0.0;
                jac.umtv(vbasis[i],vtransformedbasis[i]);
              }

            // compute sigma
            sigma=0.0;
            for (std::size_t i=0; i<velocityspace.size(); i++)
              sigma.axpy(x(velocityspace,i),vtransformedbasis[i]);

            // K^{-1} * sigma
            tensor.mv(sigma,Kinvsigma);

            // integrate  (K^{-1}*sigma) * phi_i
            auto factor = ip.weight() / det;
            for (std::size_t i=0; i<velocityspace.size(); i++)
              r.accumulate(velocityspace,i,(Kinvsigma*vtransformedbasis[i])*factor);
          }

        // u div v term, div sigma q term, a0*u term
        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder_p))
          {
            // evaluate shape functions at ip (this is a Galerkin method)
            velocityspace.finiteElement().localBasis().evaluateJacobian(ip.position(),vjacbasis);
            pressurespace.finiteElement().localBasis().evaluateFunction(ip.position(),pbasis);

            // compute u
            PressureRangeType u;
            u=0.0;
            for (std::size_t i=0; i<pressurespace.size(); i++)
              u.axpy(x(pressurespace,i),pbasis[i]);

            // evaluate Helmholtz term (reaction term)
            auto a0value = param.c(cell,ip.position());

            // integrate a0 * u * q
            RF factor = ip.weight();
            for (std::size_t i=0; i<pressurespace.size(); i++)
              r.accumulate(pressurespace,i,-a0value*u*pbasis[i]*factor);

            // compute divergence of velocity basis functions on reference element
            for (std::size_t i=0; i<velocityspace.size(); i++){
              divergence[i] = 0;
              for (int j=0; j<dim; j++)
                divergence[i] += vjacbasis[i][j][j];
            }

            // integrate sigma * phi_i
            for (std::size_t i=0; i<velocityspace.size(); i++)
              r.accumulate(velocityspace,i,-u*divergence[i]*factor);

            // compute divergence of sigma
            RF divergencesigma = 0.0;
            for (std::size_t i=0; i<velocityspace.size(); i++)
              divergencesigma += x(velocityspace,i)*divergence[i];

            // integrate div sigma * q
            for (std::size_t i=0; i<pressurespace.size(); i++)
              r.accumulate(pressurespace,i,-divergencesigma*pbasis[i]*factor);
          }
      }

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // Define types
        using PressureSpace = typename LFSV::template Child<1>::Type;
        using PressureRangeType = typename PressureSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;

        // select the pressure component
        using namespace Indices;
        const auto& pressurespace = child(lfsv,_1);

        // References to the cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // Initialize vectors outside for loop
        std::vector<PressureRangeType> pbasis(pressurespace.size());

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder_p))
          {
            // evaluate shape functions
            pressurespace.finiteElement().localBasis().evaluateFunction(ip.position(),pbasis);

            // evaluate right hand side parameter function
            auto y = param.f(cell,ip.position());

            // integrate f
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            for (std::size_t i=0; i<pressurespace.size(); i++)
              r.accumulate(pressurespace,i,y*pbasis[i]*factor);
          }
      }

      // boundary integral independen of ansatz functions
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
        // Define types
        using VelocitySpace = typename LFSV::template Child<0>::Type;
        using DF = typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType;
        using VelocityRangeType = typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;

        // select the two velocity component
        using namespace Indices;
        const auto& velocityspace = child(lfsv,_0);

        // dimensions
        const int dim = IG::dimension;

        // References to the inside cell
        const auto& cell_inside = ig.inside();

        // Get geometry
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // Get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // evaluate transformation which must be linear
        Dune::FieldVector<DF,dim> pos;
        pos = 0.0;
        typename IG::Entity::Geometry::JacobianInverseTransposed jac;
        jac = geo_inside.jacobianInverseTransposed(pos);
        jac.invert();
        auto det = geo_inside.integrationElement(pos);

        // Declare vectors outside for loop
        std::vector<VelocityRangeType> vbasis(velocityspace.size());
        std::vector<VelocityRangeType> vtransformedbasis(velocityspace.size());

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : quadratureRule(geo,qorder_v))
          {
            // evaluate boundary condition type
            auto bctype = param.bctype(ig.intersection(),ip.position());

            // skip rest if we are on Neumann boundary
            if (bctype == ConvectionDiffusionBoundaryConditions::Neumann)
              continue;

            // position of quadrature point in local coordinates of element
            auto local = geo_in_inside.global(ip.position());

            // evaluate test shape functions
            velocityspace.finiteElement().localBasis().evaluateFunction(local,vbasis);

            // transform basis vectors
            for (std::size_t i=0; i<velocityspace.size(); i++)
              {
                vtransformedbasis[i] = 0.0;
                jac.umtv(vbasis[i],vtransformedbasis[i]);
              }

            // evaluate Dirichlet boundary condition
            auto y = param.g(cell_inside,local);

            // integrate g v*normal
            auto factor = ip.weight()*geo.integrationElement(ip.position())/det;
            for (std::size_t i=0; i<velocityspace.size(); i++)
              r.accumulate(velocityspace,i,y*(vtransformedbasis[i]*ig.unitOuterNormal(ip.position()))*factor);
          }
      }

    private:
      const PARAM& param;
      int qorder_v;
      int qorder_p;
    };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_DIFFUSIONMIXED_HH
