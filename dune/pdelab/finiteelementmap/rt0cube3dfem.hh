// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_RT0CUBE3DFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_RT0CUBE3DFEM_HH

#include<vector>
#include<dune/localfunctions/raviartthomas/raviartthomas0cube3d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap

    template<typename GV, typename D, typename R>
    class RT0Cube3DLocalFiniteElementMap :
      public LocalFiniteElementMapInterface<LocalFiniteElementMapTraits< Dune::RT0Cube3DLocalFiniteElement<D,R> >,
                                            RT0Cube3DLocalFiniteElementMap<GV,D,R> >
    {
      typedef Dune::RT0Cube3DLocalFiniteElement<D,R> FE;
      typedef typename GV::IndexSet IndexSet;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! \brief Use when Imp has a standard constructor
      RT0Cube3DLocalFiniteElementMap (const GV& gv_)
        : gv(gv_), is(gv_.indexSet()), orient(gv_.size(0))
      {
        // create all variants
        for (int i=0; i<64; i++)
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

      bool hasDOFs(int codim) const
      {
        return codim == 1;
      }

      std::size_t size(GeometryType gt) const
      {
        return gt.dim() == 2 && gt.isCube() ? 1 : 0;
      }

      std::size_t maxLocalSize() const
      {
        return 6;
      }

    private:
      GV gv;
      FE variant[64];
      const IndexSet& is;
      std::vector<unsigned char> orient;
    };
  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RT0CUBE3DFEM_HH
