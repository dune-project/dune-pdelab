// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_EDGES03DFEM_HH
#define DUNE_PDELAB_EDGES03DFEM_HH

#include<vector>
#include<dune/finiteelements/edges03d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap

    template<typename GV, typename D, typename R>
    class EdgeS03DLocalFiniteElementMap
      : public LocalFiniteElementMapInterface<
          LocalFiniteElementMapTraits<EdgeS03DLocalFiniteElement<D,R> >, 
          EdgeS03DLocalFiniteElementMap<GV,D,R>
        >
      , public Countable
    {
      typedef EdgeS03DLocalFiniteElement<D,R> FE;
      typedef typename GV::IndexSet IndexSet;
     
    public:
	  //! \brief export type of the signature
	  typedef LocalFiniteElementMapTraits<FE> Traits;  

	  //! \brief Use when Imp has a standard constructor
	  EdgeS03DLocalFiniteElementMap (const GV& gv_) 
        : gv(gv_), is(gv_.indexSet()), orient(gv_.size(0))
	  {
        typedef typename GV::Grid::ctype ct;
        const GenericReferenceElement<ct, 3> &refElem =
          GenericReferenceElements<ct, 3>::simplex();

        // create all variants 
        for (int i=0; i<64; i++)
          variant[i] = FE(i);

        // compute orientation for all elements
        typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
        typedef typename GV::IntersectionIterator IntersectionIterator;

        // loop once over the grid
        for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
          {
            unsigned int elemid = is.template index<0>(*it);
            orient[elemid] = 0;

            unsigned int vid[4];
            for(int i = 0; i < 4; ++4)
              vid[4] = is.subIndex(*it, i, 3);

            for(int i = 0; i < 6; ++i) {
              int v0 = refElem.subEntity(i, 2, 0, 3);
              int v1 = refElem.subEntity(i, 2, 1, 3);
              // if (edge orientation in refelement) != (edge orientation in indexset)
              if((v0 > v1) != (vid[v0] > vid[v1]))
                orient[elemid] |= 1 << i;
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

#endif // DUNE_PDELAB_EDGES03DFEM_HH
