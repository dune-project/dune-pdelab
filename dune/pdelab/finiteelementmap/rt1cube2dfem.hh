// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_RT1CUBE2DFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_RT1CUBE2DFEM_HH

#include <vector>
#include <dune/localfunctions/raviartthomas/raviartthomas1cube2d.hh>
#include "finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<typename GV, typename D, typename R>
    class RT1Cube2DLocalFiniteElementMap :
      public LocalFiniteElementMapInterface<
        LocalFiniteElementMapTraits< Dune::RT1Cube2DLocalFiniteElement<D,R> >,
        RT1Cube2DLocalFiniteElementMap<GV,D,R> >
    {
      typedef Dune::RT1Cube2DLocalFiniteElement<D,R> FE;
      typedef typename GV::IndexSet IndexSet;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! \brief Use when Imp has a standard constructor
      RT1Cube2DLocalFiniteElementMap(const GV& gv_)
        : gv(gv_), is(gv_.indexSet()), orient(gv_.size(0))
      {
        // create all variants
        for (int i = 0; i < 16; i++)
        {
          variant[i] = FE(i);
        }

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
      const typename Traits::FiniteElementType& find(const EntityType& e) const
      {
        return variant[orient[is.index(e)]];
      }

      bool fixedSize() const
      {
        return true;
      }

      bool hasDOFs(int codim) const
      {
        return codim == 0 || codim == 1;
      }

      std::size_t size(GeometryType gt) const
      {
        switch (gt.dim())
          {
          case 2:
            return 4;
          case 1:
            return 2;
          default:
            return 0;
          }
      }

      std::size_t maxLocalSize() const
      {
        return 12;
      }

    private:
      GV gv;
      FE variant[16];
      const IndexSet& is;
      std::vector<unsigned char> orient;
    };
  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RT1CUBE2DFEM_HH
