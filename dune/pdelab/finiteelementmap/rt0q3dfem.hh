// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_RT0Q3DFEM_HH
#define DUNE_PDELAB_RT0Q3DFEM_HH

#include<vector>
#include<dune/localfunctions/rt0q3d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap

    template<typename GV, typename D, typename R>
    class RT0Q3DLocalFiniteElementMap : 
	  public LocalFiniteElementMapInterface<LocalFiniteElementMapTraits< Dune::RT0Q3DLocalFiniteElement<D,R> >, 
											RT0Q3DLocalFiniteElementMap<GV,D,R> >,
      public Countable
    {
      typedef Dune::RT0Q3DLocalFiniteElement<D,R> FE;
      typedef typename GV::IndexSet IndexSet;
     
    public:
	  //! \brief export type of the signature
	  typedef LocalFiniteElementMapTraits<FE> Traits;  

	  //! \brief Use when Imp has a standard constructor
	  RT0Q3DLocalFiniteElementMap (const GV& gv_) 
        : gv(gv_), is(gv_.indexSet()), orient(gv_.size(0))
	  {
        // create all variants 
        for (int i=0; i<64; i++)
          variant[i] = FE(i);

        // compute orientation for all elements
        typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
        typedef typename GV::IntersectionIterator IntersectionIterator;

        // loop once over the grid
        for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
          {
            unsigned int myid = is.template index<0>(*it);
            orient[myid] = 0;

            IntersectionIterator endit = gv.iend(*it);
            for (IntersectionIterator iit = gv.ibegin(*it); iit!=endit; ++iit)
              if (iit->neighbor())
                {
                  if (is.template index<0>(*(iit->outside()))>myid)
                    orient[myid] |= 1<<iit->indexInInside();
                }
          }
      }

	  //! \brief get local basis functions for entity
	  template<class EntityType>
	  const typename Traits::LocalFiniteElementType& find (const EntityType& e) const
	  {
        return variant[orient[is.template index<0>(e)]];
	  }

	private:
      const GV& gv;
      FE variant[64];
      const IndexSet& is;
      std::vector<unsigned char> orient;
    };
  }
}

#endif
