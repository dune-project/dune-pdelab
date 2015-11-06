// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_RT1CUBE3DFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_RT1CUBE3DFEM_HH

#include <vector>
#include <dune/localfunctions/raviartthomas/raviartthomas1cube3d.hh>
#include "finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<typename GV, typename D, typename R>
    class RT1Cube3DLocalFiniteElementMap :
      public RTLocalFiniteElementMap<
        GV,
        Dune::RT1Cube3DLocalFiniteElement<D,R>,
        RT1Cube3DLocalFiniteElementMap<GV,D,R>,
        64>
    {
      typedef Dune::RT1Cube3DLocalFiniteElement<D,R> FE;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! \brief Use when Imp has a standard constructor
      RT1Cube3DLocalFiniteElementMap (const GV& gv)
        : RTLocalFiniteElementMap<
          GV,
          Dune::RT1Cube3DLocalFiniteElement<D,R>,
          RT1Cube3DLocalFiniteElementMap<GV,D,R>,
          64>(gv)
      {}

      bool fixedSize() const
      {
        return true;
      }

      bool hasDOFs(int codim) const
      {
        return (codim == 0 || codim == 1);
      }

      std::size_t size(GeometryType gt) const
      {
        switch (gt.dim())
          {
          case 3:
            return 12;
          case 2:
            return 4;
          default:
            return 0;
          }
      }

      std::size_t maxLocalSize() const
      {
        return 36;
      }

    };
  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RT1CUBE3DFEM_HH
