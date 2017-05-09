// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_EDGES0_5FEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_EDGES0_5FEM_HH

#include <dune/localfunctions/whitney/edges0.5.hh>

#include <dune/pdelab/finiteelementmap/global.hh>

namespace Dune {
  namespace PDELab {

    //! Global-valued finite element map for EdgeS0_5 elements
    /**
     * \ingroup FiniteElementMap
     *
     * \tparam Geometry           Type of the geometry od the elements.
     * \tparam VertexOrderFactory Type of factory for extracting vertex
     *                            ordering information.
     * \tparam RF                 Range field type.
     */
    template<class Geometry, class VertexOrderFactory, class RF>
    class EdgeS0_5FiniteElementMap :
      public GeometryVertexOrderFiniteElementMap<
        EdgeS0_5FiniteElementFactory<Geometry, RF>, VertexOrderFactory
        >
    {
      typedef EdgeS0_5FiniteElementFactory<Geometry, RF> FEFactory;
      typedef GeometryVertexOrderFiniteElementMap<
        FEFactory, VertexOrderFactory
        > Base;

      static FEFactory &feFactory() {
        static FEFactory feFactory_;
        return feFactory_;
      }

    public:
      EdgeS0_5FiniteElementMap(const VertexOrderFactory &voFactory) :
        Base(feFactory(), voFactory)
      { }

      bool fixedSize() const
      {
        return true;
      }

      bool hasDOFs(int codim) const
      {
        return Geometry::mydimension - codim == 1;
      }

      std::size_t size(GeometryType gt) const
      {
        return gt.isLine() ? 1 : 0;
      }

      std::size_t maxLocalSize() const
      {
        return Dune::EdgeS0_5Common<Geometry::mydimension>::s;
      }

    };
  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_EDGES0_5FEM_HH
