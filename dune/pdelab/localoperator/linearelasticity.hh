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

    class LinearElasticity : public NumericalJacobianApplyVolume<LinearElasticity>,
                             public FullVolumePattern,
                             public LocalOperatorDefaultFlags,
                             public InstationaryLocalOperatorDefaultMethods<double>,
                             public NumericalJacobianVolume<LinearElasticity>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doLambdaVolume = true };
      enum { doLambdaBoundary = false };

      LinearElasticity (int intorder_=2)
        : intorder(intorder_), mu(1.0), lambda(1.0)
      {}

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
              
              // integrate 1/2 \mu (grad u)^2 * (grad phi_i)^2
              for (int k=0; k<dim; k++)
                a += 0.5 * mu * (gradu[k] * gradphi[i][k] + gradu[k] * gradphi[i][d]);
              // integrate \lambda sum_(k=0..dim) (d u / d x_d) * (d phi_i / d x_k)
              for (int k=0; k<dim; k++)
                b += lambda * gradu[d]*gradphi[i][k];

              r.accumulate(lfsu,i, (a+b) * factor);
            }
          }
        }
      }

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // // domain and range field type
        // typedef typename LFSV::Traits::FiniteElementType::
        //   Traits::LocalBasisType::Traits::DomainFieldType DF;
        // typedef typename LFSV::Traits::FiniteElementType::
        //   Traits::LocalBasisType::Traits::RangeFieldType RF;
        // typedef typename LFSV::Traits::FiniteElementType::
        //   Traits::LocalBasisType::Traits::RangeType RangeType;

        // typedef typename LFSV::Traits::SizeType size_type;
        
        // // dimensions
        // const int dim = EG::Geometry::dimension;
        
        // // select quadrature rule
        // GeometryType gt = eg.geometry().type();
        // const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);

        // // loop over quadrature points
        // for (typename QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        // {
        //   // evaluate shape functions 
        //   std::vector<RangeType> phi(lfsv.size());
        //   lfsv.finiteElement().localBasis().evaluateFunction(it->position(),phi);

        //   // evaluate right hand side parameter function
        //   typename F::Traits::RangeType y;
        //   f.evaluate(eg.entity(),it->position(),y);

        //   // integrate f
        //   RF factor = it->weight() * eg.geometry().integrationElement(it->position());
        //   for (size_type i=0; i<lfsv.size(); i++)
        //     r[i] -= y*phi[i]*factor;
        // }
      }

    private:
      int intorder;
      double mu;
      double lambda;
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
