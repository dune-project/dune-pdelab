// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_Q1FEM_HH
#define DUNE_PDELAB_Q1FEM_HH

#include <cstddef>

#include<dune/localfunctions/lagrange/q1.hh>

#include <dune/pdelab/finiteelementmap/global.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R, int d>
    class Q1LocalFiniteElementMap
      : public SimpleLocalFiniteElementMap< Dune::Q1LocalFiniteElement<D,R,d> >
    {

    public:

      bool fixedSize() const
      {
        return true;
      }

      std::size_t size(GeometryType gt) const
      {
        return gt.isVertex() ? 1 : 0;
      }

      std::size_t maxLocalSize() const
      {
        return 1 << d;
      }

    };

    //! Global-valued finite element map for Q1 elements
    /**
     * \ingroup FiniteElementMap
     *
     * \tparam Geometry Type of geometries of the elements.  Only used to
     *                  extract the dimension and the type for the domain
     *                  field.
     * \tparam RF       Range field type.
     */
    template<class Geometry, class RF>
    class Q1FiniteElementMap
      : public GeometryFiniteElementMap<Q1FiniteElementFactory<Geometry, RF> >
    {
      typedef Q1FiniteElementFactory<Geometry, RF> FEFactory;
      typedef GeometryFiniteElementMap<FEFactory> Base;

      static FEFactory feFactory;

    public:
      Q1FiniteElementMap() : Base(feFactory) { }

      bool fixedSize() const
      {
        return true;
      }

      std::size_t size(GeometryType gt) const
      {
        return gt.isVertex() ? 1 : 0;
      }

      std::size_t maxLocalSize() const
      {
        return 1 << Geometry::dimension;
      }
    };

    template<class GV, class RF>
    typename Q1FiniteElementMap<GV, RF>::FEFactory
    Q1FiniteElementMap<GV, RF>::feFactory;
  }
}

#endif
