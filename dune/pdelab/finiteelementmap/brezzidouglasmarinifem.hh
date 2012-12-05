// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITELEMENTMAP_BREZZIDOUGLASMARINIFEM_HH
#define DUNE_PDELAB_FINITELEMENTMAP_BREZZIDOUGLASMARINIFEM_HH

#include <dune/common/static_assert.hh>
#include <dune/grid/common/capabilities.hh>

#include <dune/pdelab/common/topologyutility.hh>
#include <dune/pdelab/finiteelementmap/bdm1simplex2dfem.hh>
#include <dune/pdelab/finiteelementmap/bdm1cube2dfem.hh>


namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    namespace detail {

      template<typename GV, int dim, GeometryType::BasicType basic_type, typename D, typename R, std::size_t k>
      struct BrezziDouglasMariniLocalFiniteElementMapBaseSelector
      {
        dune_static_assert((AlwaysFalse<GV>::value),"The requested type of Brezzi-Douglas-Marini element is not implemented, sorry!");
      };


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
      dune_static_assert((AlwaysFalse<GV>::value),
                         "Your chosen grid does not export a usable topology id for its cells."
                         "Please provide the correct GeometryType::BasicType as an additional template parameter.");
    };

#endif // DOXYGEN

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FINITELEMENTMAP_BREZZIDOUGLASMARINIFEM_HH
