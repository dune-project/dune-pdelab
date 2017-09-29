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
        LocalFiniteElementMapTraits<Dune::BDM1Cube2DLocalFiniteElement<D,R> >,
        BDM1Cube2DLocalFiniteElementMap<GV,D,R> >
    {
      typedef Dune::BDM1Cube2DLocalFiniteElement<D,R> FE;
      typedef typename GV::IndexSet IndexSet;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! The dimension of the finite elements returned by this map.
      static constexpr int dimension = GV::dimension;

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
        //--------------------------------------------
        // loop once over the grid
        for (const auto& cell : elements(gv)) {
          unsigned int myId = is.template index<0>(cell);
          orient[myId] = 0;

          for (const auto& intersection : intersections(gv,cell)) {
            if (intersection.neighbor()
                && is.template index<0>(intersection.outside()) > myId)
            {
              orient[myId] |= 1 << intersection.indexInInside();
            }
          }
        }
      }

      //! \brief get local basis functions for entity
      template<class EntityType>
      const typename Traits::FiniteElementType& find(const EntityType& e) const
      {
        return variant[orient[is.index(e)]];
      }

      static constexpr bool fixedSize()
      {
        return true;
      }

      static constexpr bool hasDOFs(int codim)
      {
        return codim == 1;
      }

      static constexpr std::size_t size(GeometryType gt)
      {
        switch (gt.dim())
          {
          case 1:
            return 2;
          default:
            return 0;
          }
      }

      static constexpr std::size_t maxLocalSize()
      {
        return 8;
      }

    private:
      GV gv;
      FE variant[16];
      const IndexSet& is;
      std::vector<unsigned char> orient;
    };
  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_BDM1CUBE2DFEM_HH
