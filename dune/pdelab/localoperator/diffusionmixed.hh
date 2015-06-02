// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_DIFFUSIONMIXED_HH
#define DUNE_PDELAB_DIFFUSIONMIXED_HH

#include <cstddef>
#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>

#include<dune/geometry/type.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>

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

      typedef typename ConvectionDiffusionBoundaryConditions::Type BCType;

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
        // select the two components
        typedef typename LFSU::template Child<0>::Type VelocitySpace;
        const VelocitySpace& velocityspace = lfsu.template child<0>();
        typedef typename LFSU::template Child<1>::Type PressureSpace;
        const PressureSpace& pressurespace = lfsu.template child<1>();

        // domain and range field type
        typedef typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType VelocityJacobianType;
        typedef typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType VelocityRangeType;
        typedef typename PressureSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType PressureRangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // evaluate transformation which must be linear
        Dune::FieldVector<DF,dim> pos; pos = 0.0;
        typename EG::Geometry::JacobianInverseTransposed jac;
        jac = eg.geometry().jacobianInverseTransposed(pos);
        jac.invert();
        RF det = eg.geometry().integrationElement(pos);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        //typename K::Traits::RangeType tensor;
        Dune::GeometryType gt = eg.geometry().type();
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);

        typename PARAM::Traits::PermTensorType tensor;
        tensor = param.A(eg.entity(),localcenter);
        tensor.invert(); // need iverse for mixed method

        // \sigma\cdot v term
        const Dune::QuadratureRule<DF,dim>& vrule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder_v);
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=vrule.begin(); it!=vrule.end(); ++it)
          {
            // evaluate shape functions at ip (this is a Galerkin method)
            std::vector<VelocityRangeType> vbasis(velocityspace.size());
            velocityspace.finiteElement().localBasis().evaluateFunction(it->position(),vbasis);

            // transform basis vectors
            std::vector<VelocityRangeType> vtransformedbasis(velocityspace.size());
            for (std::size_t i=0; i<velocityspace.size(); i++)
              {
                vtransformedbasis[i] = 0.0;
                jac.umtv(vbasis[i],vtransformedbasis[i]);
              }

            // compute sigma
            VelocityRangeType sigma; sigma = 0.0;
            for (std::size_t i=0; i<velocityspace.size(); i++)
              sigma.axpy(x(velocityspace,i),vtransformedbasis[i]);

            // K^{-1} * sigma
            VelocityRangeType Kinvsigma;
            tensor.mv(sigma,Kinvsigma);

            // integrate  (K^{-1}*sigma) * phi_i
            RF factor = it->weight() / det;
            for (std::size_t i=0; i<velocityspace.size(); i++)
              r.accumulate(velocityspace,i,(Kinvsigma*vtransformedbasis[i])*factor);
          }

        // u div v term, div sigma q term, a0*u term
        const Dune::QuadratureRule<DF,dim>& prule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder_p);
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=prule.begin(); it!=prule.end(); ++it)
          {
            // evaluate shape functions at ip (this is a Galerkin method)
            std::vector<VelocityJacobianType> vbasis(velocityspace.size());
            velocityspace.finiteElement().localBasis().evaluateJacobian(it->position(),vbasis);
            std::vector<PressureRangeType> pbasis(pressurespace.size());
            pressurespace.finiteElement().localBasis().evaluateFunction(it->position(),pbasis);

            // compute u
            PressureRangeType u; u = 0.0;
            for (std::size_t i=0; i<pressurespace.size(); i++)
              u.axpy(x(pressurespace,i),pbasis[i]);

            // evaluate Helmholtz term (reaction term)
            typename PARAM::Traits::RangeFieldType a0value = param.c(eg.entity(),it->position());

            // integrate a0 * u * q
            RF factor = it->weight();
            for (std::size_t i=0; i<pressurespace.size(); i++)
              r.accumulate(pressurespace,i,-a0value*u*pbasis[i]*factor);

            // compute divergence of velocity basis functions on reference element
            std::vector<RF> divergence(velocityspace.size(),0.0);
            for (std::size_t i=0; i<velocityspace.size(); i++)
              for (int j=0; j<dim; j++)
                divergence[i] += vbasis[i][j][j];

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
        // select the two components
        typedef typename LFSV::template Child<1>::Type PressureSpace;
        const PressureSpace& pressurespace = lfsv.template child<1>();

        // domain and range field type
        typedef typename PressureSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename PressureSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename PressureSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType PressureRangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder_p);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate shape functions
            std::vector<PressureRangeType> pbasis(pressurespace.size());
            pressurespace.finiteElement().localBasis().evaluateFunction(it->position(),pbasis);

            // evaluate right hand side parameter function
            RF y = param.f(eg.entity(),it->position());

            // integrate f
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (std::size_t i=0; i<pressurespace.size(); i++)
              r.accumulate(pressurespace,i,y*pbasis[i]*factor);
          }
      }

      // boundary integral independen of ansatz functions
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<0>::Type VelocitySpace;
        const VelocitySpace& velocityspace = lfsv.template child<0>();

        // domain and range field type
        typedef typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename VelocitySpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType VelocityRangeType;

        // dimensions
        const int dim = IG::dimension;

        auto inside_cell = ig.inside();

        // evaluate transformation which must be linear
        Dune::FieldVector<DF,dim> pos; pos = 0.0;
        typename IG::Entity::Geometry::JacobianInverseTransposed jac;
        jac = inside_cell.geometry().jacobianInverseTransposed(pos);
        jac.invert();
        RF det = inside_cell.geometry().integrationElement(pos);

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder_v);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate boundary condition type
            BCType bctype = param.bctype(ig.intersection(),it->position());
            // skip rest if we are on Neumann boundary

            if (bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann)
            //if (bctype == ConvectionDiffusionBoundaryConditions::Neumann)
              continue;

            // position of quadrature point in local coordinates of element
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

            // evaluate test shape functions
            std::vector<VelocityRangeType> vbasis(velocityspace.size());
            velocityspace.finiteElement().localBasis().evaluateFunction(local,vbasis);

            // transform basis vectors
            std::vector<VelocityRangeType> vtransformedbasis(velocityspace.size());
            for (std::size_t i=0; i<velocityspace.size(); i++)
              {
                vtransformedbasis[i] = 0.0;
                jac.umtv(vbasis[i],vtransformedbasis[i]);
              }

            // evaluate Dirichlet boundary condition
            RF y = param.g(inside_cell,local);

            // integrate g v*normal
            RF factor = it->weight()*ig.geometry().integrationElement(it->position())/det;
            for (std::size_t i=0; i<velocityspace.size(); i++)
              r.accumulate(velocityspace,i,y*(vtransformedbasis[i]*ig.unitOuterNormal(it->position()))*factor);
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

#endif
