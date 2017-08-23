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
      public RTLocalFiniteElementMap<
        GV,
        Dune::RT1Cube2DLocalFiniteElement<D,R>,
        RT1Cube2DLocalFiniteElementMap<GV,D,R>,
        16>
    {
      typedef Dune::RT1Cube2DLocalFiniteElement<D,R> FE;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! \brief Use when Imp has a standard constructor
      RT1Cube2DLocalFiniteElementMap(const GV& gv)
        : RTLocalFiniteElementMap<
          GV,
          Dune::RT1Cube2DLocalFiniteElement<D,R>,
          RT1Cube2DLocalFiniteElementMap<GV,D,R>,
          16>(gv)
      {}

      static constexpr bool fixedSize()
      {
        return true;
      }

      static constexpr bool hasDOFs(int codim)
      {
        return codim == 0 || codim == 1;
      }

      static constexpr std::size_t size(GeometryType gt)
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

      static constexpr std::size_t maxLocalSize()
      {
        return 12;
      }

    };
  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RT1CUBE2DFEM_HH
