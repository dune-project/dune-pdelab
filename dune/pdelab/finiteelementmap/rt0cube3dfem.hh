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
      public RTLocalFiniteElementMap<
        GV,
        Dune::RT0Cube3DLocalFiniteElement<D,R>,
        RT0Cube3DLocalFiniteElementMap<GV,D,R>,
        64>
    {
      typedef Dune::RT0Cube3DLocalFiniteElement<D,R> FE;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! \brief Use when Imp has a standard constructor
      RT0Cube3DLocalFiniteElementMap (const GV& gv)
        : RTLocalFiniteElementMap<
          GV,
          Dune::RT0Cube3DLocalFiniteElement<D,R>,
          RT0Cube3DLocalFiniteElementMap<GV,D,R>,
          64>(gv)
      {}

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
        return gt == GeometryTypes::quadrilateral ? 1 : 0;
      }

      static constexpr std::size_t maxLocalSize()
      {
        return 6;
      }

    };
  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RT0CUBE3DFEM_HH
