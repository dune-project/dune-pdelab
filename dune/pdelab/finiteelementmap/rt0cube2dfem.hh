// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_RT0CUBE2DFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_RT0CUBE2DFEM_HH

#include<vector>
#include<dune/localfunctions/raviartthomas/raviartthomas0cube2d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap

    template<typename GV, typename D, typename R>
    class RT0Cube2DLocalFiniteElementMap :
      public LocalFiniteElementMapInterface<LocalFiniteElementMapTraits< Dune::RT0Cube2DLocalFiniteElement<D,R> >,
                                            RT0Cube2DLocalFiniteElementMap<GV,D,R> >
    {
      typedef Dune::RT0Cube2DLocalFiniteElement<D,R> FE;
      typedef typename GV::IndexSet IndexSet;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! \brief Use when Imp has a standard constructor
      RT0Cube2DLocalFiniteElementMap (const GV& gv_)
        : gv(gv_), is(gv_.indexSet()), orient(gv_.size(0))
      {
        // create all variants
        for (int i=0; i<16; i++)
          variant[i] = FE(i);

        // compute orientation for all elements
        for (const auto& element : elements(gv))
        {
          unsigned int myId = is.index(element);
          orient[myId] = 0;

          for (const auto& intersection : intersections(gv,element))
          {
            if (intersection.neighbor()
                && is.index(intersection.outside()) > myId)
            {
              orient[myId] |= 1 << intersection.indexInInside();
            }
          }
        }
      }

      //! \brief get local basis functions for entity
      template<class EntityType>
      const typename Traits::FiniteElementType& find (const EntityType& e) const
      {
        return variant[orient[is.index(e)]];
      }

      bool fixedSize() const
      {
        return true;
      }

      std::size_t size(GeometryType gt) const
      {
        return gt.isLine() ? 1 : 0;
      }

      std::size_t maxLocalSize() const
      {
        return 4;
      }

    private:
      GV gv;
      FE variant[16];
      const IndexSet& is;
      std::vector<unsigned char> orient;
    };
  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RT0CUBE2DFEM_HH
