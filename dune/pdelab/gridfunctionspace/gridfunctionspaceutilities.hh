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

    //===============================================================
    // output: convert grid function space to discrete grid function
    //===============================================================


	// convert a single component function space into a grid function
    // the functions can be vector-valued
    // this is just an intermediate solution to provide VTK output
	template<typename T, typename X>
	class DiscreteGridFunction
	  : public GridFunctionInterface<GridFunctionTraits<typename T::Traits::GridViewType,
														typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
														T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
														typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType>,
									 DiscreteGridFunction<T,X> >
	{
	  typedef T GFS;

	  typedef GridFunctionInterface<GridFunctionTraits<typename T::Traits::GridViewType,
														typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
														T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
														typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType>,
                                    DiscreteGridFunction<T,X> > BaseT;
	
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


	// convert a power function space of scalar function spaces into a vector-valued grid function
    // this is just an intermediate solution to provide VTK output
	template<typename T, typename X>
	class VectorDiscreteGridFunction
	  : public GridFunctionInterface<GridFunctionTraits<typename T::Traits::GridViewType,
                typename T::template Child<0>::Type::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
                T::CHILDREN,
                Dune::FieldVector<typename T::template Child<0>::Type::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,T::CHILDREN> >,
			   VectorDiscreteGridFunction<T,X> >
	{
	  typedef T GFS;

      typedef GridFunctionInterface<GridFunctionTraits<typename T::Traits::GridViewType,
                typename T::template Child<0>::Type::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
                T::CHILDREN,
                Dune::FieldVector<typename T::template Child<0>::Type::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,T::CHILDREN> >,
                                    VectorDiscreteGridFunction<T,X> > BaseT;

	public:
	  typedef typename BaseT::Traits Traits;
      typedef typename T::template Child<0>::Type ChildType;
      typedef typename ChildType::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
      typedef typename ChildType::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType RT;
	  
	  VectorDiscreteGridFunction (const GFS& gfs, const X& x_)
		: pgfs(&gfs), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), yb(gfs.maxLocalSize())
	  {
	  }

	  inline void evaluate (const typename Traits::ElementType& e, 
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {  
		lfs.bind(e);
		lfs.vread(xg,xl);
        for (int k=0; k<T::CHILDREN; k++)
          {
            lfs.getChild(k).localFiniteElement().localBasis().evaluateFunction(x,yb);
            y[k] = 0.0;
            for (unsigned int i=0; i<yb.size(); i++)
              y[k] += xl[lfs.getChild(k).localIndex(i)]*yb[i];
          }
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
	  mutable std::vector<RF> xl;
	  mutable std::vector<RT> yb;
	};

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
