// -*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_LINEARELASTICITY_HH
#define DUNE_PDELAB_LOCALOPERATOR_LINEARELASTICITY_HH

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/quadraturerules.hh>

#include "defaultimp.hh"
#include "pattern.hh"
#include "flags.hh"
#include "idefault.hh"

#include "linearelasticityparameter.hh"

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** a local operator for solving the the linear elasticity problem using conforming FEM
     *
     * \note we only support Dirichlet and homogeneous Neumann boundary conditions
     *
     * \tparam T model of LinearElasticityParameterInterface
     *
     * \todo check LFSU size
     * \todo check LFSU == LFSV
     */
    template<typename T>
    class LinearElasticity : public FullVolumePattern,
                             public LocalOperatorDefaultFlags,
                             public InstationaryLocalOperatorDefaultMethods<typename T::Traits::DomainType>
    {
    public:

      using ParameterType = T;

      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doLambdaVolume = true };
      enum { doLambdaBoundary = true };

      LinearElasticity (const ParameterType & p, int intorder=4)
        : intorder_(intorder), param_(p)
      {}

      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat) const
      {
        // Define types
        using LFSU_SUB = TypeTree::Child<LFSU,0>;
        using RF = typename M::value_type;
        using JacobianType = typename LFSU_SUB::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType;
        using size_type = typename LFSU_SUB::Traits::SizeType;

        // dimensions
        const int dim = EG::Entity::dimension;
        const int dimw = EG::Geometry::coorddimension;
        static_assert(dim == dimw, "doesn't work on manifolds");

        // Reference to cell
        const auto& cell = eg.entity();

        // get geometry
        auto geo = eg.geometry();

        // Transformation
        typename EG::Geometry::JacobianInverseTransposed jac;

        // Initialize vectors outside for loop
        std::vector<JacobianType> js(lfsu.child(0).size());
        std::vector<FieldVector<RF,dim> > gradphi(lfsu.child(0).size());

        // loop over quadrature points
        for (const auto& qp : quadratureRule(geo,intorder_))
        {
          // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
          lfsu.child(0).finiteElement().localBasis().evaluateJacobian(qp.position(),js);

          // transform gradient to real element
          jac = geo.jacobianInverseTransposed(qp.position());
          for (size_type i=0; i<lfsu.child(0).size(); i++)
          {
            gradphi[i] = 0.0;
            jac.umv(js[i][0],gradphi[i]);
          }

          // material parameters
          auto mu = param_.mu(cell,qp.position());
          auto lambda = param_.lambda(cell,qp.position());

          // geometric weight
          auto factor = qp.weight() * geo.integrationElement(qp.position());

          for(int d=0; d<dim; ++d)
          {
            for (size_type i=0; i<lfsu.child(0).size(); i++)
            {
              for (int k=0; k<dim; k++)
              {
                for (size_type j=0; j<lfsv.child(k).size(); j++)
                {
                  // integrate \mu (grad u + (grad u)^T) * (grad phi_i + (grad phi_i)^T)
                  mat.accumulate(lfsv.child(k),j,lfsu.child(k),i,
                    // mu (d u_k / d x_d) (d v_k / d x_d)
                    mu * gradphi[i][d] * gradphi[j][d]
                    *factor);
                  mat.accumulate(lfsv.child(k),j,lfsu.child(d),i,
                    // mu (d u_d / d x_k) (d v_k / d x_d)
                    mu * gradphi[i][k] * gradphi[j][d]
                    *factor);
                  // integrate \lambda sum_(k=0..dim) (d u_d / d x_d) * (d v_k / d x_k)
                  mat.accumulate(lfsv.child(k),j,lfsu.child(d),i,
                    lambda * gradphi[i][d]*gradphi[j][k]
                    *factor);
                }
              }
            }
          }
        }
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU_HAT, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU_HAT& lfsu_hat, const X& x, const LFSV& lfsv, R& r) const
      {
        // Define types
        using LFSU = TypeTree::Child<LFSU_HAT,0>;
        using RF = typename R::value_type;
        using JacobianType = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType;
        using size_type = typename LFSU::Traits::SizeType;

        // dimensions
        const int dim = EG::Entity::dimension;
        const int dimw = EG::Geometry::coorddimension;
        static_assert(dim == dimw, "doesn't work on manifolds");

        // Reference to cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // Transformation
        typename EG::Geometry::JacobianInverseTransposed jac;

        // Initialize vectors outside for loop
        std::vector<JacobianType> js(lfsu_hat.child(0).size());
        std::vector<FieldVector<RF,dim> > gradphi(lfsu_hat.child(0).size());
        Dune::FieldVector<RF,dim> gradu(0.0);

        // loop over quadrature points
        for (const auto& qp : quadratureRule(geo,intorder_))
        {
          // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
          lfsu_hat.child(0).finiteElement().localBasis().evaluateJacobian(qp.position(),js);

          // transform gradient to real element
          jac = geo.jacobianInverseTransposed(qp.position());
          for (size_type i=0; i<lfsu_hat.child(0).size(); i++)
          {
            gradphi[i] = 0.0;
            jac.umv(js[i][0],gradphi[i]);
          }

          // material parameters
          auto mu = param_.mu(cell,qp.position());
          auto lambda = param_.lambda(cell,qp.position());

          // geometric weight
          auto factor = qp.weight() * geo.integrationElement(qp.position());

          for(int d=0; d<dim; ++d)
          {
            const LFSU & lfsu = lfsu_hat.child(d);

            // compute gradient of u
            gradu = 0.0;
            for (size_t i=0; i<lfsu.size(); i++)
            {
              gradu.axpy(x(lfsu,i),gradphi[i]);
            }

            for (size_type i=0; i<lfsv.child(d).size(); i++)
            {
              for (int k=0; k<dim; k++)
              {
                // integrate \mu (grad u + (grad u)^T) * (grad phi_i + (grad phi_i)^T)
                r.accumulate(lfsv.child(d),i,
                  // mu (d u_d / d x_k) (d phi_i_d / d x_k)
                  mu * gradu[k] * gradphi[i][k]
                  *factor);
                r.accumulate(lfsv.child(k),i,
                  // mu (d u_d / d x_k) (d phi_i_k / d x_d)
                  mu * gradu[k] * gradphi[i][d]
                  *factor);
                // integrate \lambda sum_(k=0..dim) (d u / d x_d) * (d phi_i / d x_k)
                r.accumulate(lfsv.child(k),i,
                  lambda * gradu[d]*gradphi[i][k]
                  *factor);
              }
            }
          }
        }
      }

      // volume integral depending only on test functions
      template<typename EG, typename LFSV_HAT, typename R>
      void lambda_volume (const EG& eg, const LFSV_HAT& lfsv_hat, R& r) const
      {
        // Define types
        using LFSV = TypeTree::Child<LFSV_HAT,0>;
        using RF = typename R::value_type;
        using RangeType = typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using size_type = typename LFSV::Traits::SizeType;

        // dimensions
        const int dim = EG::Entity::dimension;

        // Reference to cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // Initialize vectors outside for loop
        std::vector<RangeType> phi(lfsv_hat.child(0).size());
        FieldVector<RF,dim> y(0.0);

        // loop over quadrature points
        for (const auto& qp : quadratureRule(geo,intorder_))
        {
          // evaluate shape functions
          lfsv_hat.child(0).finiteElement().localBasis().evaluateFunction(qp.position(),phi);

          // evaluate right hand side parameter function
          y = 0.0;
          param_.f(cell,qp.position(),y);

          // weight
          auto factor = qp.weight() * geo.integrationElement(qp.position());

          for(int d=0; d<dim; ++d)
          {
            const auto& lfsv = lfsv_hat.child(d);

            // integrate f
            for (size_type i=0; i<lfsv.size(); i++)
              r.accumulate(lfsv,i, -y[d]*phi[i] * factor);
          }
        }
      }

      // jacobian of boundary term
      template<typename IG, typename LFSV_HAT, typename R>
      void lambda_boundary (const IG& ig, const LFSV_HAT& lfsv_hat, R& r) const
      {
        // Define types
        using namespace Indices;
        using LFSV = TypeTree::Child<LFSV_HAT,0>;
        using RF = typename R::value_type;
        using RangeType = typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using size_type = typename LFSV::Traits::SizeType;

        // dimensions
        const int dim = IG::Entity::dimension;

        // get geometry
        auto geo = ig.geometry();

        // Get geometry of intersection in local coordinates of inside cell
        auto geo_in_inside = ig.geometryInInside();

        // Initialize vectors outside for loop
        std::vector<RangeType> phi(lfsv_hat.child(0).size());
        FieldVector<RF,dim> y(0.0);

        // loop over quadrature points
        for (const auto& qp : quadratureRule(geo,intorder_))
        {
          // position of quadrature point in local coordinates of element
          auto local = geo_in_inside.global(qp.position());

          // evaluate boundary condition type
          // skip rest if we are on Dirichlet boundary
          if( param_.isDirichlet( ig.intersection(), qp.position() ) )
            continue;

          // evaluate shape functions
          lfsv_hat.child(0).finiteElement().localBasis().evaluateFunction(local,phi);

          // evaluate surface force
          y = 0.0;
          // currently we only implement homogeneous Neumann (e.g. Stress) BC
          // param_.g(eg.entity(),qp.position(),y);

          // weight
          auto factor = qp.weight() * geo.integrationElement(qp.position());

          for(int d=0; d<dim; ++d)
          {
            const auto& lfsv = lfsv_hat.child(d);

            // integrate f
            for (size_type i=0; i<lfsv.size(); i++)
              r.accumulate(lfsv,i, y[d]*phi[i] * factor);
          }
        }
      }

    protected:
      int intorder_;
      const ParameterType & param_;
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_LINEARELASTICITY_HH
