// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_PK2DFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_PK2DFEM_HH

#warning dune/pdelab/finiteelementmap/pk2dfem.hh and Pk2DLocalFiniteElementMap are deprecated, please use dune/pdelab/finiteelementmap/pkfem.hh and PkLocalFiniteElementMap instead

#include <dune/common/deprecated.hh>

#include <dune/localfunctions/lagrange/pk2d.hh>

#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/global.hh>

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap

    template<typename GV, typename D, typename R, unsigned int k>
    class DUNE_DEPRECATED_MSG("Please use PkLocalFiniteElementMap instead") Pk2DLocalFiniteElementMap
      : public fem::PkLocalFiniteElementMapBase<GV,D,R,k,2>
    {

    public:

      Pk2DLocalFiniteElementMap(const GV& gv)
        : fem::PkLocalFiniteElementMapBase<GV,D,R,k,2>(gv)
        {}

    };

    //! Global-valued finite element map for Pk2D elements
    /**
     * \ingroup FiniteElementMap
     *
     * \tparam Geometry           Type of the geometry od the elements.
     * \tparam VertexOrderFactory Type of factory for extracting vertex
     *                            ordering information.
     * \tparam RF                 Range field type.
     * \tparam k                  Order of the elements.
     */
    template<class Geometry, class VertexOrderFactory, class RF, std::size_t k>
    class DUNE_DEPRECATED_MSG("Please use PkLocalFiniteElementMap instead") Pk2DFiniteElementMap :
      public GeometryVertexOrderFiniteElementMap<
      Pk2DFiniteElementFactory<Geometry, RF, k>, VertexOrderFactory
      >
    {
      typedef Pk2DFiniteElementFactory<Geometry, RF, k> FEFactory;
      typedef GeometryVertexOrderFiniteElementMap<
        FEFactory, VertexOrderFactory
        > Base;

      static FEFactory &feFactory() {
        static FEFactory feFactory_;
        return feFactory_;
      }

    public:
      Pk2DFiniteElementMap(const VertexOrderFactory &voFactory) :
        Base(feFactory(), voFactory)
        { }

      bool fixedSize() const
      {
        return true;
      }

      std::size_t size(GeometryType gt) const
      {
        if (gt.isVertex())
          return k > 0 ? 1 : 0;
        if (gt.isLine())
          return k > 1 ? k - 1 : 0;
        if (gt.isTriangle())
          return k > 2 ? (k-2)*(k-1)/2 : (k == 0);
        return 0;
      }

      std::size_t maxLocalSize() const
      {
        return (k+1)*(k+2)/2;
      }

    };
  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_PK2DFEM_HH
