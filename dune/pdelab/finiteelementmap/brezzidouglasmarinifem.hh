// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_BREZZIDOUGLASMARINIFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_BREZZIDOUGLASMARINIFEM_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/pdelab/common/topologyutility.hh>
#include <dune/pdelab/finiteelementmap/bdm1simplex2dfem.hh>
#include <dune/pdelab/finiteelementmap/bdm1cube2dfem.hh>


namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    namespace detail {

      //! Template for determining the correct base implementation of BrezziDouglasMariniLocalFiniteElementMap.
      /**
       * This template must be specialized for supported variations of the Brezzi-Douglas-Marini elements.
       * Its member type 'type' will be used as base class of BrezziDouglasMariniLocalFiniteElementMap. That
       * type must have a constructor that takes the GridView as single parameter.
       */
      template<typename GV, int dim, GeometryType::BasicType basic_type, typename D, typename R, std::size_t k>
      struct BrezziDouglasMariniLocalFiniteElementMapBaseSelector
      {
        static_assert((AlwaysFalse<GV>::value),"The requested type of Brezzi-Douglas-Marini element is not implemented, sorry!");
      };


      // ********************************************************************************
      // Specializations
      // ********************************************************************************

      template<typename GV, typename D, typename R>
      struct BrezziDouglasMariniLocalFiniteElementMapBaseSelector<GV,2,GeometryType::simplex,D,R,1>
      {
        typedef BDM1Simplex2DLocalFiniteElementMap<GV,D,R> type;
      };

      template<typename GV, typename D, typename R>
      struct BrezziDouglasMariniLocalFiniteElementMapBaseSelector<GV,2,GeometryType::cube,D,R,1>
      {
        typedef BDM1Cube2DLocalFiniteElementMap<GV,D,R> type;
      };

    } // end namespace detail

#endif // DOXYGEN


    //! \ingroup FiniteElementMap
    //! Brezzi-Douglas-Marini elements of order k.
    /**
     * This LocalFiniteElementMap provides BDM elements of order k for the given
     * GridView GV. It currently supports BDM1 for both simplices and cubes in 2D.
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
    class BrezziDouglasMariniLocalFiniteElementMap :
      public detail::BrezziDouglasMariniLocalFiniteElementMapBaseSelector<GV,GV::dimension,basic_type,D,R,k>::type
    {

    public:

      //! Constructs a finite element map on the GridView gv.
      BrezziDouglasMariniLocalFiniteElementMap(const GV& gv)
        : detail::BrezziDouglasMariniLocalFiniteElementMapBaseSelector<GV,GV::dimension,basic_type,D,R,k>::type(gv)
      {}

    };

#ifndef DOXYGEN

    // Specialization for grids that don't provide a valid topology id for their cells.
    template<typename GV, typename D, typename R, std::size_t k>
    class BrezziDouglasMariniLocalFiniteElementMap<GV,D,R,k,GeometryType::none>
    {
      static_assert((AlwaysFalse<GV>::value),
                    "Your chosen grid does not export a usable topology id for its cells."
                    "Please provide the correct GeometryType::BasicType as an additional template parameter.");
    };

#endif // DOXYGEN

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_BREZZIDOUGLASMARINIFEM_HH
