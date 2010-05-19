// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_L2_HH
#define DUNE_PDELAB_L2_HH

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

    /** a local operator for the mass operator (L_2 integral)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     */
	class L2 : public NumericalJacobianApplyVolume<L2>,
               public FullVolumePattern,
               public LocalOperatorDefaultFlags,
               public InstationaryLocalOperatorDefaultMethods<double>
	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };

	  // residual assembly flags
      enum { doAlphaVolume = true };

      L2 (int intorder_=2)
        : intorder(intorder_)
      {}

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSU::Traits::SizeType size_type;
        
        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.localFiniteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x[i]*phi[i];

            // u*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsu.size(); i++)
              r[i] += u*phi[i]*factor;
          }
	  }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
                            LocalMatrix<R>& mat) const
      {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {            
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.localFiniteElement().localBasis().evaluateFunction(it->position(),phi);

            // integrate phi_j*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat(i,j) += phi[j]*phi[i]*factor;
          }
      }

    private:
      int intorder;
	};

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
