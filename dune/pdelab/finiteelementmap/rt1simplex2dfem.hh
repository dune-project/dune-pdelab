// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_RT1SIMPLEX2DFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_RT1SIMPLEX2DFEM_HH

#include <vector>
#include <dune/localfunctions/raviartthomas/raviartthomas12d.hh>
#include "finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<typename GV, typename D, typename R>
    class RT1Simplex2DLocalFiniteElementMap :
      public RTLocalFiniteElementMap<
        GV,
        Dune::RT12DLocalFiniteElement<D,R>,
        RT1Simplex2DLocalFiniteElementMap<GV,D,R>,
        8>
    {
      typedef Dune::RT12DLocalFiniteElement<D,R> FE;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! \brief Use when Imp has a standard constructor
      RT1Simplex2DLocalFiniteElementMap(const GV& gv)
        : RTLocalFiniteElementMap<
          GV,
          Dune::RT12DLocalFiniteElement<D,R>,
          RT1Simplex2DLocalFiniteElementMap<GV,D,R>,
          8>(gv)
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
            return 2;
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

    };
  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RT1SIMPLEX2DFEM_HH
