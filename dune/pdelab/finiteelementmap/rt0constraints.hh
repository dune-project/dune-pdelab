// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_RT0CONSTRAINTS_HH
#define DUNE_PDELAB_RT0CONSTRAINTS_HH

#include<dune/common/exceptions.hh>
#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/grid/common/grid.hh>
#include<dune/common/geometrytype.hh>
#include<dune/pdelab/common/geometrywrapper.hh>

namespace Dune {
  namespace PDELab {

	class RT0Constraints {
	public:
	  enum{doBoundary=true};enum{doProcessor=false};
	  enum{doSkeleton=false};enum{doVolume=false};
	  
	  template<typename B, typename I, typename LFS, typename T>
	  void boundary (const B& b, const I& ig, const LFS& lfs, T& trafo) const
	  {
        typedef typename IntersectionGeometry<I>::ctype DT;
        const int dim = IntersectionGeometry<I>::Entity::Geometry::dimension;
        const int face = ig.indexInInside();
	    const Dune::GenericReferenceElement<DT,dim-1> & 
	      face_refelem = Dune::GenericReferenceElements<DT,dim-1>::general(ig.geometry().type()); 
        const typename B::Traits::DomainType ip = face_refelem.position(0,0);
		typename B::Traits::RangeType bctype;   // return value
		b.evaluate(ig,ip,bctype);               // eval condition type
		if (bctype>0) return;                   // done
		typename T::RowType empty;              // need not interpolate
		trafo[face]=empty;
	  }
	};

  }
}

#endif
