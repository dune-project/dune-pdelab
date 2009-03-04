// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CONSTRAINTS_HH
#define DUNE_PDELAB_CONSTRAINTS_HH

#include<dune/common/exceptions.hh>
#include<dune/grid/common/referenceelements.hh>

#include"../common/function.hh"
#include"../common/geometrywrapper.hh"

#include"gridfunctionspace.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    // do method invocation only if class has the method

    template<typename C, bool doIt>
    struct ConstraintsCallBoundary
    {
      template<typename F, typename I, typename LFS, typename T>
      static void boundary (const C& c, const F& f, const IntersectionGeometry<I>& ig, const LFS& lfs, T& trafo)
      {
      }
    };
    template<typename C, bool doIt>
    struct ConstraintsCallSkeleton
    {
      template<typename I, typename LFS, typename T>
      static void skeleton (const C& c, const IntersectionGeometry<I>& ig, const LFS& lfs, T& trafo)
      {
      }
    };
    template<typename C, bool doIt>
    struct ConstraintsCallVolume
    {
      template<typename E, typename LFS, typename T>
      static void volume (const C& c, const ElementGeometry<E>& eg, const LFS& lfs, T& trafo)
      {
      }
    };


    template<typename C>
    struct ConstraintsCallBoundary<C,true>
    {
      template<typename F, typename I, typename LFS, typename T>
      static void boundary (const C& c, const F& f, const IntersectionGeometry<I>& ig, const LFS& lfs, T& trafo)
      {
        c.boundary(f,ig,lfs,trafo);
      }
    };
    template<typename C>
    struct ConstraintsCallSkeleton<C,true>
    {
      template<typename I, typename LFS, typename T>
      static void skeleton (const C& c, const IntersectionGeometry<I>& ig, const LFS& lfs, T& trafo)
      {
        c.skeleton(ig,lfs,trafo);
      }
    };
    template<typename C>
    struct ConstraintsCallVolume<C,true>
    {
      template<typename E, typename LFS, typename T>
      static void volume (const C& c, const ElementGeometry<E>& eg, const LFS& lfs, T& trafo)
      {
        c.volume(eg,lfs,trafo);
      }
    };

    // meta program to evaluate boundary constraints

	template<typename F, bool FisLeaf, typename LFS, bool LFSisLeaf> 
	struct ConstraintsVisitNodeMetaProgram;

	template<typename F, typename LFS, int n, int i> 
	struct ConstraintsVisitChildMetaProgram // visit i'th child of inner node
	{
	  template<typename CG, typename I>
	  static void boundary (const F& f, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
	  {
        // vist children of both nodes in pairs
		typedef typename F::template Child<i>::Type FC;
		typedef typename LFS::template Child<i>::Type LFSC;

        const FC& fc=f.template getChild<i>();
        const LFSC& lfsc=lfs.template getChild<i>();

        ConstraintsVisitNodeMetaProgram<FC,FC::isLeaf,LFSC,LFSC::isLeaf>::boundary(fc,lfsc,cg,ig);
		ConstraintsVisitChildMetaProgram<F,LFS,n,i+1>::boundary(f,lfs,cg,ig);
	  }
	};

	template<typename F, typename LFS, int n> 
	struct ConstraintsVisitChildMetaProgram<F,LFS,n,n> // end of child recursion
	{
      // end of child recursion
	  template<typename CG, typename I>
	  static void boundary (const F& f, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
	  {
        return;
	  }
	};

	template<typename F, bool FisLeaf, typename LFS, bool LFSisLeaf> 
	struct ConstraintsVisitNodeMetaProgram // visit inner node
	{
	  template<typename CG, typename I>
	  static void boundary (const F& f, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
	  {
        // both are inner nodes, visit all children
        // check that both have same number of children
		dune_static_assert((static_cast<int>(F::CHILDREN)==static_cast<int>(LFS::CHILDREN)),  
						   "both nodes must have same number of children");
 
        // start child recursion
		ConstraintsVisitChildMetaProgram<F,LFS,F::CHILDREN,0>::boundary(f,lfs,cg,ig);
	  }
	};

	template<typename F, typename LFS> 
	struct ConstraintsVisitNodeMetaProgram<F,true,LFS,false> // try to interpolate components from vector valued function
	{
	  template<typename CG, typename I>
	  static void boundary (const F& f, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
	  {
		dune_static_assert((static_cast<int>(LFS::isPower)==1),  
						   "specialization only for power");
		dune_static_assert((static_cast<int>(LFS::template Child<0>::Type::isLeaf)==1),  
						   "children must be leaves");
		dune_static_assert((static_cast<int>(F::Traits::dimRange)==static_cast<int>(LFS::CHILDREN)),  
						   "number of components must coincide with number of children");

        // extract constraints type 
        typedef typename LFS::template Child<0>::Type::Traits::ConstraintsType C;

        for (int k=0; k<LFS::CHILDREN; k++)
          {
            // allocate empty local constraints map
            CG cl;

            // call boundary condition evaluation of child k with component k
            typedef BoundaryGridFunctionSelectComponentAdapter<F> FCOMP;
            FCOMP fcomp(f,k);

            ConstraintsCallBoundary<C,C::doBoundary>::boundary(lfs.getChild(k).constraints(),
                                                               fcomp,ig,lfs.getChild(k),cl);

            // write coefficients into local vector 
            lfs.getChild(k).mwrite(cl,cg);
          }
	  }
	};

	template<typename F, typename LFS> 
	struct ConstraintsVisitNodeMetaProgram<F,true,LFS,true> // leaf node in both trees 
	{
	  template<typename CG, typename I>
	  static void boundary (const F& f, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
	  {
        // now we are at a single component local function space
        // which is part of a multi component local function space

		// allocate local constraints map
		CG cl;

        // extract constraints type 
        typedef typename LFS::Traits::ConstraintsType C;

		// iterate over boundary, need intersection iterator
        ConstraintsCallBoundary<C,C::doBoundary>::boundary(lfs.constraints(),f,ig,lfs,cl);

		// write coefficients into local vector 
		lfs.mwrite(cl,cg);
	  }
	};


    // second metaprogram that iterates over local function space only

	template<typename LFS, bool LFSisLeaf> 
	struct ConstraintsVisitNodeMetaProgram2;

	template<typename LFS, int n, int i> 
	struct ConstraintsVisitChildMetaProgram2 // visit i'th child of inner node
	{
	  template<typename CG, typename I>
	  static void skeleton (const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
	  {
		typedef typename LFS::template Child<i>::Type LFSC;
        const LFSC& lfsc=lfs.template getChild<i>();

        ConstraintsVisitNodeMetaProgram2<LFSC,LFSC::isLeaf>::skeleton(lfsc,cg,ig);
		ConstraintsVisitChildMetaProgram2<LFS,n,i+1>::skeleton(lfs,cg,ig);
	  }
	  template<typename CG, typename E>
	  static void volume (const LFS& lfs, CG& cg, const ElementGeometry<E>& eg)
	  {
		typedef typename LFS::template Child<i>::Type LFSC;
        const LFSC& lfsc=lfs.template getChild<i>();

        ConstraintsVisitNodeMetaProgram2<LFSC,LFSC::isLeaf>::volume(lfsc,cg,eg);
		ConstraintsVisitChildMetaProgram2<LFS,n,i+1>::volume(lfs,cg,eg);
	  }
	};

	template<typename LFS, int n> 
	struct ConstraintsVisitChildMetaProgram2<LFS,n,n> // end of child recursion
	{
 	  template<typename CG, typename I>
	  static void skeleton (const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
	  {
        return;
	  }
	  template<typename CG, typename E>
	  static void volume (const LFS& lfs, CG& cg, const ElementGeometry<E>& eg)
	  {
        return;
	  }
	};


	template<typename LFS, bool LFSisLeaf> 
	struct ConstraintsVisitNodeMetaProgram2 // visit inner node
	{
	  template<typename CG, typename I>
	  static void skeleton (const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
	  {
        // start child recursion
		ConstraintsVisitChildMetaProgram2<LFS,LFS::CHILDREN,0>::skeleton(lfs,cg,ig);
	  }
	  template<typename CG, typename E>
	  static void volume (const LFS& lfs, CG& cg, const ElementGeometry<E>& eg)
	  {
        // start child recursion
		ConstraintsVisitChildMetaProgram2<LFS,LFS::CHILDREN,0>::volume(lfs,cg,eg);
	  }
	};


	template<typename LFS> 
	struct ConstraintsVisitNodeMetaProgram2<LFS,true> // leaf node
	{
	  template<typename CG, typename I>
	  static void skeleton (const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
	  {
        // now we are at a single component local function space
        // which is part of a multi component local function space

		// allocate local constraints map
		CG cl;

        // extract constraints type 
        typedef typename LFS::Traits::ConstraintsType C;

		// iterate over boundary, need intersection iterator
        ConstraintsCallSkeleton<C,C::doSkeleton>::skeleton(lfs.constraints(),ig,lfs,cl);

		// write coefficients into local vector 
		lfs.mwrite(cl,cg);
	  }
	  template<typename CG, typename E>
	  static void volume (const LFS& lfs, CG& cg, const ElementGeometry<E>& eg)
	  {
        // now we are at a single component local function space
        // which is part of a multi component local function space

		// allocate local constraints map
		CG cl;

        // extract constraints type 
        typedef typename LFS::Traits::ConstraintsType C;

		// iterate over boundary, need intersection iterator
        ConstraintsCallVolume<C,C::doVolume>::volume(lfs.constraints(),eg,lfs,cl);

		// write coefficients into local vector 
		lfs.mwrite(cl,cg);
	  }
	};


    // construct constraints from given boundary condition function
    template<typename F, typename GFS, typename CG>
    void constraints (F& f, const GFS& gfs, CG& cg)
    {
      // clear global constraints
	  cg.clear();

      // get some types
      typedef typename GFS::Traits::GridViewType GV;
      typedef typename GV::Traits::template Codim<0>::Entity Element;
      typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
	  typedef typename GV::IntersectionIterator IntersectionIterator;
	  typedef typename IntersectionIterator::Intersection Intersection;

      // make local function space
      typedef typename GFS::LocalFunctionSpace LFS;
      LFS lfs(gfs);

      // loop once over the grid
      for (ElementIterator it = gfs.gridview().template begin<0>();
           it!=gfs.gridview().template end<0>(); ++it)
        {
          // bind local function space to element
          lfs.bind(*it);

          ConstraintsVisitNodeMetaProgram2<LFS,LFS::isLeaf>
            ::volume(lfs,cg,ElementGeometry<Element>(*it));
          
		  // iterate over intersections and call metaprogram
		  IntersectionIterator endit = gfs.gridview().iend(*it);
		  for (IntersectionIterator iit = gfs.gridview().ibegin(*it); iit!=endit; ++iit)
			{
			  if (iit->boundary())
				ConstraintsVisitNodeMetaProgram<F,F::isLeaf,LFS,LFS::isLeaf>
				  ::boundary(f,lfs,cg,IntersectionGeometry<Intersection>(*iit));
			  if (iit->neighbor())
				ConstraintsVisitNodeMetaProgram2<LFS,LFS::isLeaf>
				  ::skeleton(lfs,cg,IntersectionGeometry<Intersection>(*iit));
			}
		}

	  // print result
	  std::cout << "constraints:" << std::endl;
	  typedef typename CG::iterator global_col_iterator;
	  typedef typename CG::value_type::second_type global_row_type;
	  typedef typename global_row_type::iterator global_row_iterator;
	  
      std::cout << cg.size() << " constrained degrees of freedom" << std::endl;

// 	  for (global_col_iterator cit=cg.begin(); cit!=cg.end(); ++cit)
// 		{
// 		  std::cout << cit->first << ": ";
// 		  for (global_row_iterator rit=(cit->second).begin(); rit!=(cit->second).end(); ++rit)
// 			std::cout << "(" << rit->first << "," << rit->second << ") ";
// 		  std::cout << std::endl;
// 		}
	}

    // construct constraints from given boundary condition function
    template<typename CG, typename XG>
    void set_constrained_dofs (const CG& cg, typename XG::ElementType x, XG& xg)
    {
	  typedef typename CG::const_iterator global_col_iterator;	  
	  for (global_col_iterator cit=cg.begin(); cit!=cg.end(); ++cit)
		xg[cit->first] = x;
	}

    // construct constraints from given boundary condition function
    template<typename CG, typename XG>
    void copy_constrained_dofs (const CG& cg, const XG& xgin, XG& xgout)
    {
	  typedef typename CG::const_iterator global_col_iterator;	  
	  for (global_col_iterator cit=cg.begin(); cit!=cg.end(); ++cit)
		xgout[cit->first] = xgin[cit->first];
	}

    // construct constraints from given boundary condition function
    template<typename CG, typename XG>
    void set_nonconstrained_dofs (const CG& cg, typename XG::ElementType x, XG& xg)
    {
      for (typename XG::size_type i=0; i<xg.size(); ++i)
        if (cg.find(i)==cg.end())
          xg[i] = x;
	}


    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
