// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_TOPOLOGYUTILITY_HH
#define DUNE_PDELAB_COMMON_TOPOLOGYUTILITY_HH

#include <dune/geometry/type.hh>

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
      static const bool isCube =
        ((topologyId ^ ((1 << dimension)-1)) >> 1 == 0);

      static const bool isSimplex =
        (topologyId | 1) == 1;

      //! The BasicType of Topology.
      static const GeometryType::BasicType value =
        isSimplex ? GeometryType::simplex
        : (
          isCube ? GeometryType::cube
          : GeometryType::none
          );
    };


    /** Extract Dune::GeometryType from Dune::GeometryType::BasicType
     *
     * Note: Dune::GeometryType used to have constructors taking a BasicType
     * and the dimension. This utility function can be used as a replacement.
     */
    inline Dune::GeometryType geometryTypeFromBasicType(Dune::GeometryType::BasicType basicType, int dim){

      unsigned int topologyId(0);

      if (dim > 1){
        switch( basicType )
        {
        case GeometryType::simplex :
          topologyId = 0;
          break;
        case GeometryType::cube :
          topologyId = ((1 << dim) - 1);
          break;
        case GeometryType::pyramid :
          if (dim == 3)
            topologyId = 0b0011;
          else
            DUNE_THROW( RangeError,
                        "Invalid basic geometry type: no pyramids for dimension " << dim << "." );
          break;
        case GeometryType::prism :
          if (dim == 3)
            topologyId = 0b0101;
          else
            DUNE_THROW( RangeError,
                        "Invalid basic geometry type: no prisms for dimension " << dim << "." );
          break;
        case GeometryType::none :
          break;
        default :
          DUNE_THROW( RangeError,
                      "Invalid basic geometry type: " << basicType << " for dimension " << dim << "." );
        }
      }

      return Dune::GeometryType(topologyId, dim);
    }

  } // namespace PDELab
} // namespace Dune


#endif // DUNE_PDELAB_COMMON_TOPOLOGYUTILITY_HH
