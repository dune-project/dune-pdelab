// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_Q22DFEM_HH
#define DUNE_PDELAB_Q22DFEM_HH

#warning dune/pdelab/finiteelementmap/q22dfem.hh and Q22DLocalFiniteElementMap are deprecated, please use dune/pdelab/finiteelementmap/qkfem.hh and QkLocalFiniteElementMap instead

#include <dune/common/deprecated.hh>
#include <dune/localfunctions/lagrange/q2.hh>

#include "finiteelementmap.hh"
#include <dune/pdelab/finiteelementmap/global.hh>

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R>
    class DUNE_DEPRECATED_MSG("Please use QkLocalFiniteElementMap instead") Q22DLocalFiniteElementMap
      : public SimpleLocalFiniteElementMap< Dune::Q2LocalFiniteElement<D,R,2> >
    {

    public:

      bool fixedSize() const
      {
        return true;
      }

      std::size_t size(GeometryType gt) const
      {
        if (gt.isVertex() || gt.isLine() || gt.isQuadrilateral())
          return 1;
        else
          return 0;
      }

      std::size_t maxLocalSize() const
      {
        return 9;
      }

    };

    //! Global-valued finite element map for Q22D elements
    /**
     * \ingroup FiniteElementMap
     *
     * \tparam Geometry Type of the geometry od the elements.
     * \tparam RF       Range field type.
     */
    template<class Geometry, class RF>
    class DUNE_DEPRECATED_MSG("Please use QkLocalFiniteElementMap instead") Q22DFiniteElementMap
      : public GeometryFiniteElementMap<
          Q2FiniteElementFactory<Geometry, RF>
          >
    {
      typedef Q2FiniteElementFactory<Geometry, RF> FEFactory;
      typedef GeometryFiniteElementMap<FEFactory> Base;

      static FEFactory feFactory;

    public:
      Q22DFiniteElementMap() : Base(feFactory) { }

      bool fixedSize() const
      {
        return true;
      }

      std::size_t size(GeometryType gt) const
      {
        if (gt.isVertex() || gt.isLine() || gt.isQuadrilateral())
          return 1;
        else
          return 0;
      }

      std::size_t maxLocalSize() const
      {
        return 9;
      }

    };

    template<class GV, class RF>
    typename Q22DFiniteElementMap<GV, RF>::FEFactory
    Q22DFiniteElementMap<GV, RF>::feFactory;
  }
}

#endif
