// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_TOPOLOGYUTILITY_HH
#define DUNE_PDELAB_COMMON_TOPOLOGYUTILITY_HH

#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/conversion.hh>

namespace Dune {

  namespace PDELab {


    //! Utility TMP for determining the BasicType of a geometry from its dimension and topology id.
    /**
     * BasicTypeFromDimensionAndTopologyId invokes GenericGeometry-internal machinery to identify
     * the GeometryType::BasicType that belongs to the given combination of dimension and topologyId.
     *
     * This information is often useful when writing generic frontends for FiniteElementMaps that are
     * specialized on dimension and BasicType.
     *
     * \tparam dimension   The dimension of the geometry.
     * \tparam topologyId  The topologyId of the geometry.
     */
    template<int dimension, unsigned int topologyId>
    struct BasicTypeFromDimensionAndTopologyId
    {

      //! The topology described by dimension and topologyId.
      typedef typename GenericGeometry::Topology<topologyId,dimension>::type Topology;

      //! The BasicType of Topology.
      static const GeometryType::BasicType value =
        GenericGeometry::DuneGeometryType<Topology,GeometryType::simplex>::basicType;
    };


  } // namespace PDELab
} // namespace Dune


#endif // DUNE_PDELAB_COMMON_TOPOLOGYUTILITY_HH
