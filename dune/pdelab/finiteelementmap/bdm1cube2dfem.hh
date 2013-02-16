// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_BDM1CUBE2DFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_BDM1CUBE2DFEM_HH

#include <vector>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube2d.hh>
#include "finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<typename GV, typename D, typename R>
    class BDM1Cube2DLocalFiniteElementMap :
      public LocalFiniteElementMapInterface<
        LocalFiniteElementMapTraits<Dune::BDM1Q2DLocalFiniteElement<D,R> >,
        BDM1Cube2DLocalFiniteElementMap<GV,D,R> >
    {
      typedef Dune::BDM12DLocalFiniteElement<D,R> FE;
      typedef typename GV::IndexSet IndexSet;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! \brief Use when Imp has a standard constructor
      BDM1Cube2DLocalFiniteElementMap(const GV& gv_)
        : gv(gv_), is(gv_.indexSet()), orient(gv_.size(0))
      {
        // create all variants
        for (int i = 0; i < 16; i++)
        {
          variant[i] = FE(i);
        }

        // compute orientation for all elements
        typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
        typedef typename GV::IntersectionIterator IntersectionIterator;

        // loop once over the grid
        for (ElementIterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
        {
          unsigned int myId = is.template index<0>(*it);
          orient[myId] = 0;

          IntersectionIterator endit = gv.iend(*it);
          for (IntersectionIterator iit = gv.ibegin(*it); iit != endit; ++iit)
          {
            if (iit->neighbor()
                && is.template index<0>(*(iit->outside())) > myId)
            {
              orient[myId] |= 1 << iit->indexInInside();
            }
          }
        }
      }

      //! \brief get local basis functions for entity
      template<class EntityType>
      const typename Traits::FiniteElementType& find(const EntityType& e) const
      {
        return variant[orient[is.template index<0>(e)]];
      }

    private:
      const GV& gv;
      FE variant[16];
      const IndexSet& is;
      std::vector<unsigned char> orient;
    };
  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_BDM1CUBE2DFEM_HH
