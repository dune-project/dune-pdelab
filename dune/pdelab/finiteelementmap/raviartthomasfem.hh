// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_RAVIARTTHOMASFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_RAVIARTTHOMASFEM_HH

#include <dune/common/static_assert.hh>
#include <dune/grid/common/capabilities.hh>

#include <dune/pdelab/finiteelementmap/rt0simplex2dfem.hh>
#include <dune/pdelab/finiteelementmap/rt1simplex2dfem.hh>
#include <dune/pdelab/finiteelementmap/rt0cube2dfem.hh>
#include <dune/pdelab/finiteelementmap/rt1cube2dfem.hh>
#include <dune/pdelab/finiteelementmap/rt2cube2dfem.hh>
#include <dune/pdelab/finiteelementmap/rt0cube3dfem.hh>
#include <dune/pdelab/finiteelementmap/rt1cube3dfem.hh>


namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    namespace detail {

      template<typename GV, int dim, unsigned int topologyId, typename D, typename R, std::size_t k>
      struct RaviartThomasLocalFiniteElementMapBaseSelector
      {
        dune_static_assert((AlwaysFalse<GV>::value),"The requested type of Raviart-Thomas element is not implemented, sorry!");
      };


      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,2,0,D,R,0>
      {
        typedef RT0Simplex2DLocalFiniteElementMap<GV,D,R> type;
      };

      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,2,0,D,R,1>
      {
        typedef RT1Simplex2DLocalFiniteElementMap<GV,D,R> type;
      };


      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,2,(1<<2)-1,D,R,0>
      {
        typedef RT0Cube2DLocalFiniteElementMap<GV,D,R> type;
      };

      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,2,(1<<2)-1,D,R,1>
      {
        typedef RT1Cube2DLocalFiniteElementMap<GV,D,R> type;
      };

      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,2,(1<<2)-1,D,R,2>
      {
        typedef RT2Cube2DLocalFiniteElementMap<GV,D,R> type;
      };


      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,3,(1<<3)-1,D,R,0>
      {
        typedef RT0Cube3DLocalFiniteElementMap<GV,D,R> type;
      };

      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,3,(1<<3)-1,D,R,1>
      {
        typedef RT1Cube3DLocalFiniteElementMap<GV,D,R> type;
      };

    } // end namespace detail

#endif // DOXYGEN


    template<typename GV, typename D, typename R, std::size_t k, unsigned int topologyId = Capabilities::hasSingleGeometryType<typename GV::Grid>::topologyId>
    class RaviartThomasLocalFiniteElementMap :
      public detail::RaviartThomasLocalFiniteElementMapBaseSelector<GV,GV::dimension,topologyId,D,R,k>::type
    {

    public:

      RaviartThomasLocalFiniteElementMap(const GV& gv)
        : detail::RaviartThomasLocalFiniteElementMapBaseSelector<GV,GV::dimension,topologyId,D,R,k>::type(gv)
      {}

    };


  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RAVIARTTHOMASFEM_HH
