// -*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LINEARELASTICITY_HH
#define DUNE_PDELAB_LINEARELASTICITY_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/grid/common/quadraturerules.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>

#include"../common/geometrywrapper.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"pattern.hh"
#include"flags.hh"
#include"idefault.hh"

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    class LinearElasticity : public FullVolumePattern,
                             public LocalOperatorDefaultFlags,
                             public InstationaryLocalOperatorDefaultMethods<double>,
                             public JacobianBasedAlphaVolume<LinearElasticity>,
                             public NumericalJacobianVolume<LinearElasticity>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doLambdaVolume = true };
      enum { doLambdaBoundary = false };

      LinearElasticity (double m, double l, int intorder_=4)
        : intorder(intorder_), mu(m), lambda(l), g(0.0)
      {}

      template<typename EG, typename LFSU_HAT, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU_HAT& lfsu_hat, const X& x, const LFSV& lfsv, M & mat) const
      {
        // extract local function spaces
        typedef typename LFSU_HAT::template Child<0>::Type LFSU;

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
        dune_static_assert(dim == dimw, "doesn't work on manifolds");

        // select quadrature rule
        GeometryType gt = eg.geometry().type();
        const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        {
          // evaluate basis functions
          std::vector<RangeType> phi(lfsu_hat.child(0).size());
          lfsu_hat.child(0).finiteElement().localBasis().evaluateFunction(it->position(),phi);
            
          // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
          std::vector<JacobianType> js(lfsu_hat.child(0).size());
          lfsu_hat.child(0).finiteElement().localBasis().evaluateJacobian(it->position(),js);
            
          // transform gradient to real element
          const FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
          std::vector<FieldVector<RF,dim> > gradphi(lfsu_hat.child(0).size());
          for (size_type i=0; i<lfsu_hat.child(0).size(); i++)
          {
            gradphi[i] = 0.0;
            jac.umv(js[i][0],gradphi[i]);
          }
            
          for(int d=0; d<dim; ++d)
          {
            const LFSU & lfsu = lfsu_hat.child(d);
            
            // geometric weight 
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());

            for (size_type i=0; i<lfsu.size(); i++)
            {
              for (int k=0; k<dim; k++)
              {
                for (size_type j=0; j<lfsv.child(k).size(); j++)
                {
                  // integrate \mu (grad u + (grad u)^T) * (grad phi_i + (grad phi_i)^T)
                  mat.accumulate(lfsv.child(d),j,lfsu,i,
                    // mu (d u_d / d x_k) (d phi_i_d / d x_k)
                    mu * gradphi[i][k] * gradphi[j][k]
                    *factor);
                  mat.accumulate(lfsv.child(k),j,lfsu,i,
                    // mu (d u_d / d x_k) (d phi_i_k / d x_d)
                    mu * gradphi[i][k] * gradphi[j][d]
                    *factor);
                  // integrate \lambda sum_(k=0..dim) (d u / d x_d) * (d phi_i / d x_k)
                  mat.accumulate(lfsv.child(k),j,lfsu,i,
                    lambda * gradphi[i][d]*gradphi[j][k]
                    *factor);
                }
              }
            }
          }
        }
      }

      // volume integral depending on test and ansatz functions
#warning TODO: check LFSU size
#warning TODO: check LFSU == LFSV
      template<typename EG, typename LFSU_HAT, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU_HAT& lfsu_hat, const X& x, const LFSV& , R& r) const
      {
        // extract local function spaces
        typedef typename LFSU_HAT::template Child<0>::Type LFSU;

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
        dune_static_assert(dim == dimw, "doesn't work on manifolds");

        // select quadrature rule
        GeometryType gt = eg.geometry().type();
        const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        {
          // evaluate basis functions
          std::vector<RangeType> phi(lfsu_hat.child(0).size());
          lfsu_hat.child(0).finiteElement().localBasis().evaluateFunction(it->position(),phi);
            
          // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
          std::vector<JacobianType> js(lfsu_hat.child(0).size());
          lfsu_hat.child(0).finiteElement().localBasis().evaluateJacobian(it->position(),js);
            
          // transform gradient to real element
          const FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
          std::vector<FieldVector<RF,dim> > gradphi(lfsu_hat.child(0).size());
          for (size_type i=0; i<lfsu_hat.child(0).size(); i++)
          {
            gradphi[i] = 0.0;
            jac.umv(js[i][0],gradphi[i]);
          }
            
          for(int d=0; d<dim; ++d)
          {
            const LFSU & lfsu = lfsu_hat.child(d);
            
            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (size_t i=0; i<lfsu.size(); i++)
            {
              gradu.axpy(x(lfsu,i),gradphi[i]);
            }

            // geometric weight 
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());

            for (size_type i=0; i<lfsu.size(); i++)
            {
              RF a = 0.0;
              RF b = 0.0;
              
              // integrate \mu (grad u + (grad u)^T) * (grad phi_i + (grad phi_i)^T)
              for (int k=0; k<dim; k++)
              {
                a += mu * gradu[k] * (gradphi[i][k] + gradphi[i][d]);
              }
              // integrate \lambda sum_(k=0..dim) (d u / d x_d) * (d phi_i / d x_k)
              for (int k=0; k<dim; k++)
                b += lambda * gradu[d]*gradphi[i][k];

              r.accumulate(lfsu,i, (a+b) * factor);
            }
          }
        }
      }

      // volume integral depending only on test functions
      template<typename EG, typename LFSV_HAT, typename R>
      void lambda_volume (const EG& eg, const LFSV_HAT& lfsv_hat, R& r) const
      {
        // extract local function spaces
        typedef typename LFSV_HAT::template Child<0>::Type LFSV;

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
        
        // select quadrature rule
        GeometryType gt = eg.geometry().type();
        const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        {
          // evaluate shape functions 
          std::vector<RangeType> phi(lfsv_hat.child(0).size());
          lfsv_hat.child(0).finiteElement().localBasis().evaluateFunction(it->position(),phi);
            
          // evaluate right hand side parameter function
          //typename F::Traits::RangeType y;
          //f.evaluate(eg.entity(),it->position(),y);
          FieldVector<RF,dim> y(0.0);
          y[dim-1] = -g;

          // weight
          RF factor = it->weight() * eg.geometry().integrationElement(it->position());
          
          for(int d=0; d<dim; ++d)
          {
            const LFSV & lfsv = lfsv_hat.child(d);
            
            // integrate f
            for (size_type i=0; i<lfsv.size(); i++)
              r.accumulate(lfsv,i, y[d]*phi[i] * factor);
          }
        }
      }

    private:
      int intorder;
      double mu;
      double lambda;
      double g;
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
