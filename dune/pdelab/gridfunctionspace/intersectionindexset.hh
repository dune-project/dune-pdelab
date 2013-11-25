// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_INTERSECTIONINDEXSET_HH
#define DUNE_PDELAB_INTERSECTIONINDEXSET_HH

//===============================================================
// Index set for intersections
// implementation only for 
//  - hanging node grids, but conforming macro grid
//  - no multiple element types !
//===============================================================

#include<vector>
#include <iostream>

#include<dune/common/exceptions.hh>

namespace Dune {
  namespace PDELab {

	template<typename GV>
	class IntersectionIndexSet
	{
	public:
	  typedef typename GV::IndexSet IndexSet;
	  typedef typename IndexSet::IndexType IndexType;
	  typedef typename GV::Intersection Intersection;
	  typedef typename GV::Traits::template Codim<0>::Entity Element; 
 
	  IntersectionIndexSet (const GV& gv_)
		: gv(gv_), is(gv.indexSet())
	  {
		typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator; 

		pre();
		for (ElementIterator it = gv.template begin<0>(); 
			 it!=gv.template end<0>(); ++it)
		  {
			visit(*it);
		  }
		post();
	  }

	  // number of intersections in index set 
	  // (intersections do not have a geometry type)
	  IndexType size () const
	  {
		return intersection_counter;
	  }

	  // number of intersections associated with given element
	  IndexType size (const Element& element) const
	  {
		return number_of_intersections[is.index(element)];
	  }

	  // get index assigned to intersection
	  IndexType index (const Intersection& intersection) const
	  {
		// on the boundary, take own codim 1 subentity
		if (intersection.boundary() && (!intersection.neighbor()))
		  return codim1index_to_intersectionindex[is.subIndex(*(intersection.inside()),intersection.indexInInside(),1)];
    
		// if we have a neighbor, take higher level
		if (intersection.neighbor())
		  {
			if (intersection.inside()->level()>=intersection.outside()->level())
			  return codim1index_to_intersectionindex[is.subIndex(*(intersection.inside()),intersection.indexInInside(),1)];
			else
			  return codim1index_to_intersectionindex[is.subIndex(*(intersection.outside()),intersection.indexInOutside(),1)];
		  }

		// we are at a processor boundary
		DUNE_THROW(Dune::Exception,"intersection index at processor boundary requested");
	  }

	  // get index of i'th intersection of element 
	  // (in order they are visited by intersection iterator)
	  IndexType subIndex (const Element& element, int i) const
	  {
		return element_intersection_subindex[entry[is.index(element)]+i];
	  }

	private:

	  // prepare loop over elements
	  void pre ()
	  {
		codim1index_to_intersectionindex.resize(gv.size(1));
		invalidIndex = gv.size(1)+1;
		for (size_t i=0; i<codim1index_to_intersectionindex.size(); ++i)
		  codim1index_to_intersectionindex[i] = invalidIndex;
		number_of_intersections.resize(gv.size(0));
		entry.resize(gv.size(0));
		element_intersection_subindex.resize(2*gv.size(1));
		intersection_counter = 0;
		oriented_intersection_counter = 0;
		std::cout << "number of codim 1 entities is " << gv.size(1) << std::endl;
	  }

	  // process given element
	  void visit (const Element& element)
	  {
		typedef typename GV::IntersectionIterator IntersectionIterator;
		size_t count = 0;
		entry[is.index(element)] = oriented_intersection_counter;
		IntersectionIterator endit = gv.iend(element);
		for (IntersectionIterator iit = gv.ibegin(element); iit!=endit; ++iit)
		  {
			if (iit->neighbor())
			  {
				IndexType c1index;
				if (iit->inside()->level()>=iit->outside()->level())
				  c1index = is.subIndex(*(iit->inside()),iit->indexInInside(),1);
				else
				  c1index = is.subIndex(*(iit->outside()),iit->indexInOutside(),1);
				if (codim1index_to_intersectionindex[c1index]==invalidIndex)
				  codim1index_to_intersectionindex[c1index]=intersection_counter++;
				element_intersection_subindex[oriented_intersection_counter] = codim1index_to_intersectionindex[c1index];
			  }
			else if (iit->boundary())
			  {
				IndexType c1index = is.subIndex(*(iit->inside()),iit->indexInInside(),1);
				if (codim1index_to_intersectionindex[c1index]==invalidIndex)
				  codim1index_to_intersectionindex[c1index]=intersection_counter++;
				element_intersection_subindex[oriented_intersection_counter] = codim1index_to_intersectionindex[c1index];
			  }
			count++;
			oriented_intersection_counter++;
		  }
		number_of_intersections[is.index(element)] = static_cast<unsigned char>(count);
	  }

	  // finalize computation of index set
	  void post ()
	  {
		std::cout << "number of oriented intersections " << oriented_intersection_counter << std::endl;
		std::cout << "number of intersections " << intersection_counter << std::endl;
	  }

      GV gv;       // our grid view
	  const IndexSet& is; // index set of the grid view

	  // we assume that the mesh is conforming in the
	  // sense that there is a codim 1 entity for each intersection
	  // the following vector assigns intersection to the codim 1 entity on the "higher" side
	  std::vector<IndexType> codim1index_to_intersectionindex;

	  // number of intersections of an element
	  std::vector<unsigned char> number_of_intersections;

	  // map (element index,interection number) to index 
	  std::vector<IndexType> element_intersection_subindex;

	  // entry point of element into element_intersection_subindex
	  std::vector<size_t> entry;

	  size_t intersection_counter;
	  size_t oriented_intersection_counter;
	  IndexType invalidIndex;
	};
  }
}

#endif
