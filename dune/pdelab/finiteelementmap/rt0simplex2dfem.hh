// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_RT0SIMPLEX2DFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_RT0SIMPLEX2DFEM_HH

#include<vector>
#include<dune/localfunctions/raviartthomas/raviartthomas02d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<typename GV, typename D, typename R>
    class RT0Simplex2DLocalFiniteElementMap :
      public RTLocalFiniteElementMap<
        GV,
        Dune::RT02DLocalFiniteElement<D,R>,
        RT0Simplex2DLocalFiniteElementMap<GV,D,R>,
        8>
    {
      typedef Dune::RT02DLocalFiniteElement<D,R> FE;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! \brief Use when Imp has a standard constructor
      RT0Simplex2DLocalFiniteElementMap (const GV& gv)
        : RTLocalFiniteElementMap<
          GV,
          Dune::RT02DLocalFiniteElement<D,R>,
          RT0Simplex2DLocalFiniteElementMap<GV,D,R>,
          8> (gv)
      {}

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
        return 3;
      }

    };
  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RT0SIMPLEX2DFEM_HH
