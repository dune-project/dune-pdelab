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

	template<typename F, bool FisLeaf, typename LFS, bool LFSisLeaf> 
	struct ConstraintsVisitNodeMetaProgram;

	template<typename F, typename LFS, int n, int i> 
	struct ConstraintsVisitChildMetaProgram // visit i'th child of inner node
	{
	  template<typename CG, typename E, typename IG>
	  static void boundary (const F& f, const LFS& lfs, CG& cg, const E& e, const IG& ig)
	  {
        // vist children of both nodes in pairs
		typedef typename F::template Child<i>::Type FC;
		typedef typename LFS::template Child<i>::Type LFSC;

        const FC& fc=f.template getChild<i>();
        const LFSC& lfsc=lfs.template getChild<i>();

        ConstraintsVisitNodeMetaProgram<FC,FC::isLeaf,LFSC,LFSC::isLeaf>::boundary(fc,lfsc,cg,e,ig);
		ConstraintsVisitChildMetaProgram<F,LFS,n,i+1>::boundary(f,lfs,cg,e,ig);
	  }
	};

	template<typename F, typename LFS, int n> 
	struct ConstraintsVisitChildMetaProgram<F,LFS,n,n> // end of child recursion
	{
      // end of child recursion
	  template<typename CG, typename E, typename IG>
	  static void boundary (const F& f, const LFS& lfs, CG& cg, const E& e, const IG& ig)
	  {
        return;
	  }
	};

	template<typename F, bool FisLeaf, typename LFS, bool LFSisLeaf> 
	struct ConstraintsVisitNodeMetaProgram // visit inner node
	{
	  template<typename CG, typename E, typename IG>
	  static void boundary (const F& f, const LFS& lfs, CG& cg, const E& e, const IG& ig)
	  {
        // both are inner nodes, visit all children
        // check that both have same number of children
		dune_static_assert((static_cast<int>(F::CHILDREN)==static_cast<int>(LFS::CHILDREN)),  
						   "both nodes must have same number of children");
 
        // start child recursion
		ConstraintsVisitChildMetaProgram<F,LFS,F::CHILDREN,0>::boundary(f,lfs,cg,e,ig);
	  }
	};

	template<typename F, typename LFS> 
	struct ConstraintsVisitNodeMetaProgram<F,true,LFS,false> // try to interpolate components from vector valued function
	{
	  template<typename CG, typename E, typename IG>
	  static void boundary (const F& f, const LFS& lfs, CG& cg, const E& e, const IG& ig)
	  {
		dune_static_assert((static_cast<int>(LFS::isPower)==1),  
						   "specialization only for power");
		dune_static_assert((static_cast<int>(LFS::template Child<0>::Type::isLeaf)==1),  
						   "children must be leaves");
		dune_static_assert((static_cast<int>(F::Traits::dimRange)==static_cast<int>(LFS::CHILDREN)),  
						   "number of components must coincide with number of children");
        for (int k=0; k<LFS::CHILDREN; k++)
          {
            // allocate empty local constraints map
            CG cl;

            // call boundary condition evaluation of child k with component k
            typedef BoundaryGridFunctionSelectComponentAdapter<F> FCOMP;
            FCOMP fcomp(f,k);
            lfs.getChild(k).constraints().boundary(fcomp,ig,lfs.getChild(k),cl);

            // write coefficients into local vector 
            lfs.getChild(k).mwrite(cl,cg);
          }
	  }
	};

	template<typename F, typename LFS> 
	struct ConstraintsVisitNodeMetaProgram<F,true,LFS,true> // leaf node in both trees 
	{
	  template<typename CG, typename E, typename IG>
	  static void boundary (const F& f, const LFS& lfs, CG& cg, const E& e, const IG& ig)
	  {
        // now we are at a single component local function space
        // which is part of a multi component local function space

		// allocate local constraints map
		CG cl;

		// iterate over boundary, need intersection iterator
		lfs.constraints().boundary(f,ig,lfs,cl);

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

		  // iterate over intersections and call metaprogram
		  IntersectionIterator endit = gfs.gridview().iend(*it);
		  for (IntersectionIterator iit = gfs.gridview().ibegin(*it); iit!=endit; ++iit)
			{
			  if (iit->boundary())
				ConstraintsVisitNodeMetaProgram<F,F::isLeaf,LFS,LFS::isLeaf>
				  ::boundary(f,lfs,cg,*it,IntersectionGeometry<Intersection>(*iit));
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
