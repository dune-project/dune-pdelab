// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_RAVIARTTHOMASFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_RAVIARTTHOMASFEM_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/pdelab/common/topologyutility.hh>
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

      //! Template for determining the correct base implementation of RaviartThomasLocalFiniteElementMap.
      /**
       * This template must be specialized for supported variations of the Raviart-Thomas elements.
       * Its member type 'type' will be used as base class of RaviartThomasLocalFiniteElementMap. That
       * type must have a constructor that takes the GridView as single parameter.
       */
      template<typename GV, int dim, GeometryType::BasicType basic_type, typename D, typename R, std::size_t k>
      struct RaviartThomasLocalFiniteElementMapBaseSelector
      {
        static_assert((AlwaysFalse<GV>::value),"The requested type of Raviart-Thomas element is not implemented, sorry!");
      };


      // ********************************************************************************
      // Specializations
      // ********************************************************************************

      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,2,GeometryType::simplex,D,R,0>
      {
        typedef RT0Simplex2DLocalFiniteElementMap<GV,D,R> type;
      };

      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,2,GeometryType::simplex,D,R,1>
      {
        typedef RT1Simplex2DLocalFiniteElementMap<GV,D,R> type;
      };


      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,2,GeometryType::cube,D,R,0>
      {
        typedef RT0Cube2DLocalFiniteElementMap<GV,D,R> type;
      };

      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,2,GeometryType::cube,D,R,1>
      {
        typedef RT1Cube2DLocalFiniteElementMap<GV,D,R> type;
      };

      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,2,GeometryType::cube,D,R,2>
      {
        typedef RT2Cube2DLocalFiniteElementMap<GV,D,R> type;
      };


      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,3,GeometryType::cube,D,R,0>
      {
        typedef RT0Cube3DLocalFiniteElementMap<GV,D,R> type;
      };

      template<typename GV, typename D, typename R>
      struct RaviartThomasLocalFiniteElementMapBaseSelector<GV,3,GeometryType::cube,D,R,1>
      {
        typedef RT1Cube3DLocalFiniteElementMap<GV,D,R> type;
      };

    } // end namespace detail

#endif // DOXYGEN


    //! Raviart-Thomas elements of order k.
    /**
     * This LocalFiniteElementMap provides Raviart-Thomas elements of order k for the given
     * GridView GV. It currently supports
     *
     * * RT0, RT1 and RT2 for cubes in 2D,
     * * RT0 and RT1 for simplices in 2D,
     * * RT0 and RT1 for cubes in 3D.
     *
     * We try to infer the type of the reference element (cube / simplex) from the GridView, but
     * that only works for grids with a single element type that is fixed at compile time. For
     * potentially mixed grids like UGGrid, you need to provide the GeometryType::BasicType of the
     * cell reference element as an additional template parameter.
     *
     * \tparam GV          The GridView on which to construct the finite element map.
     * \tparam D           The domain field type of the elements.
     * \tparam R           The range field type of the elements.
     * \tparam k           The order of the finite elements.
     * \tparam basic_type  The GeometryType::BasicType of the grid cells. You only need to provide this
     *                     template parameter for mixed grids (if you don't provide the parameter for a
     *                     mixed grid, you get a compiler error).
     */
    template<typename GV,
             typename D,
             typename R,
             std::size_t k,
             GeometryType::BasicType basic_type = BasicTypeFromDimensionAndTopologyId<
               GV::dimension,
               Capabilities::hasSingleGeometryType<typename GV::Grid>::topologyId
               >::value
             >
    class RaviartThomasLocalFiniteElementMap :
      public detail::RaviartThomasLocalFiniteElementMapBaseSelector<GV,GV::dimension,basic_type,D,R,k>::type
    {

    public:

      //! The dimension of the finite elements returned by this map.
      static constexpr int dimension = GV::dimension;

      //! Constructs a finite element map on the GridView gv.
      RaviartThomasLocalFiniteElementMap(const GV& gv)
        : detail::RaviartThomasLocalFiniteElementMapBaseSelector<GV,GV::dimension,basic_type,D,R,k>::type(gv)
      {}

    };

#ifndef DOXYGEN

    // Specialization for grids that don't provide a valid topology id for their cells.
    template<typename GV, typename D, typename R, std::size_t k>
    class RaviartThomasLocalFiniteElementMap<GV,D,R,k,GeometryType::none>
    {
      static_assert((AlwaysFalse<GV>::value),
                    "Your chosen grid does not export a usable topology id for its cells."
                    "Please provide the correct GeometryType::BasicType as an additional template parameter.");
    };

#endif // DOXYGEN

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_RAVIARTTHOMASFEM_HH
