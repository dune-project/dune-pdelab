// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_RT2CUBE2DFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_RT2CUBE2DFEM_HH

#include <vector>
#include <dune/localfunctions/raviartthomas/raviartthomas2cube2d.hh>
#include "finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<typename GV, typename D, typename R>
    class RT2Cube2DLocalFiniteElementMap :
      public RTLocalFiniteElementMap<
        GV,
        Dune::RT2Cube2DLocalFiniteElement<D,R>,
        RT2Cube2DLocalFiniteElementMap<GV,D,R>,
        16>
    {
      typedef Dune::RT2Cube2DLocalFiniteElement<D,R> FE;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! \brief Use when Imp has a standard constructor
      RT2Cube2DLocalFiniteElementMap (const GV& gv)
        : RTLocalFiniteElementMap<
          GV,
          Dune::RT2Cube2DLocalFiniteElement<D,R>,
          RT2Cube2DLocalFiniteElementMap<GV,D,R>,
          16>(gv)
      {}

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
            return 12;
          case 1:
            return 3;
          default:
            return 0;
          }
      }

      std::size_t maxLocalSize() const
      {
        return 24;
      }

    };
  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RT2CUBE2DFEM_HH
