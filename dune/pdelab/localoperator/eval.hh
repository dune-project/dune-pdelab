// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_EVAL_HH
#define DUNE_PDELAB_LOCALOPERATOR_EVAL_HH

#include <vector>

#include <dune/common/fmatrix.hh>

namespace Dune {
  namespace PDELab {

    /** two simple functions used again and again to evaluate
     *  the function value or the gradient vector of the solution on an element
     */

    template< class DomainType,
              class LFS,
              class LV,
              class RangeFieldType>
    inline void evalFunction( const DomainType& location,  // expects element local coordinates!
                              const LFS& lfs,
                              const LV& xlocal,
                              RangeFieldType &valU
                              ){

      typedef typename LFS::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeType RangeType;  // normally = Dune::FieldVector<double,1>

      std::vector<RangeType> phi(lfs.size());
      lfs.finiteElement().localBasis().evaluateFunction( location, phi );

      // evaluate u
      valU = RangeFieldType(0.0);
      for( std::size_t i=0; i<lfs.size(); i++ )
        valU += xlocal( lfs, i ) * phi[i];
    }




    template< class DomainType,
              class EG,
              class LFS,
              class LV,
              class RangeType
              >
    inline void evalGradient( const DomainType& location, // expects element local coordinates!
                              const EG& eg,
                              const LFS& lfs,
                              const LV& xlocal,
                              RangeType &gradu
                              ){

      typedef typename LFS::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::JacobianType JacobianType;

      typedef typename LFS::Traits::SizeType size_type;

      std::vector<JacobianType> js(lfs.size());
      lfs.finiteElement().localBasis().evaluateJacobian( location, js );

      enum { dim = RangeType::dimension };
      // transformation
      typename EG::Geometry::JacobianInverseTransposed jac;

      // transform gradients of shape functions to real element
      jac = eg.geometry().jacobianInverseTransposed( location );
      std::vector<RangeType> gradphi(lfs.size());
      for (size_type i=0; i<lfs.size(); i++)
        jac.mv( js[i][0], gradphi[i] );

      // compute gradient of u
      gradu = RangeType(0);
      for (size_type i=0; i<lfs.size(); i++)
        gradu.axpy( xlocal(lfs,i), gradphi[i] );

    }


  }
}
#endif // DUNE_PDELAB_LOCALOPERATOR_EVAL_HH
