// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_RT2Q2DFEM_HH
#define DUNE_PDELAB_RT2Q2DFEM_HH

#include <vector>
#include <dune/localfunctions/raviartthomas/raviartthomas2q2d.hh>
#include "finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<typename GV, typename D, typename R>
    class RT2Q2DLocalFiniteElementMap :
      public LocalFiniteElementMapInterface<
        LocalFiniteElementMapTraits< Dune::RT2Q2DLocalFiniteElement<D,R> >,
        RT2Q2DLocalFiniteElementMap<GV,D,R> >
    {
      typedef Dune::RT2Q2DLocalFiniteElement<D,R> FE;
      typedef typename GV::IndexSet IndexSet;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! \brief Use when Imp has a standard constructor
      RT2Q2DLocalFiniteElementMap (const GV& gv_)
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
      const typename Traits::FiniteElementType& find (const EntityType& e) const
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

#endif // DUNE_PDELAB_RT2Q2DFEM_HH
