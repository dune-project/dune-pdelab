// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACEUTILITIES_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACEUTILITIES_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/referenceelements.hh>

#include"../common/countingptr.hh"
#include"../common/multitypetree.hh"
#include"../common/cpstoragepolicy.hh"
#include"../common/function.hh"

#include"gridfunctionspace.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

	// convert a single component function space into a grid function
	template<typename T, typename X>
	class DiscreteGridFunction
	{
	};

	// 
	template<typename GV, typename LFEM, typename B, typename P, typename X>
	class DiscreteGridFunction<GridFunctionSpace<GV,LFEM,B,P>,X>
	  : public GridFunctionInterface<GridFunctionTraits<GV,
														typename LFEM::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
														LFEM::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
														typename LFEM::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType>,
									 DiscreteGridFunction<GridFunctionSpace<GV,LFEM,B,P>,X> >
	{
	  typedef GridFunctionSpace<GV,LFEM,B,P> GFS;

	  typedef GridFunctionInterface<GridFunctionTraits<GV,
													   typename LFEM::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
													   LFEM::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
													   typename LFEM::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType>,
									DiscreteGridFunction<GridFunctionSpace<GV,LFEM,B,P>,X> > BaseT;
	
	public:
	  typedef typename BaseT::Traits Traits;
	  
	  DiscreteGridFunction (const GFS& gfs, const X& x_)
		: pgfs(&gfs), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), yb(gfs.maxLocalSize())
	  {
	  }

	  inline void evaluate (const typename Traits::ElementType& e, 
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {  
		lfs.bind(e);
		lfs.vread(xg,xl);
		lfs.localFiniteElement().localBasis().evaluateFunction(x,yb);
		y = 0;
		for (unsigned int i=0; i<yb.size(); i++)
		  y.axpy(xl[i],yb[i]);
	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView ()
	  {
		return pgfs->gridview();
	  }

	private:
	  CP<GFS const> pgfs;
	  const X& xg;
	  mutable typename GFS::LocalFunctionSpace lfs;
	  mutable std::vector<typename Traits::RangeFieldType> xl;
	  mutable std::vector<typename Traits::RangeType> yb;
	};



    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
