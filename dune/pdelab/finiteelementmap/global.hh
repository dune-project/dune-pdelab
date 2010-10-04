// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_GLOBAL_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_GLOBAL_HH

#include <dune/pdelab/common/countingptr.hh>
#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

namespace Dune {
  namespace PDELab {

    template<class Factory>
    class GeometryFiniteElementMap :
      public Countable
    {
      Factory& factory;

    public:
      typedef FiniteElementMapTraits<typename Factory::FiniteElement> Traits;

      GeometryFiniteElementMap(Factory& factory_) : factory(factory_) {}

      //! Return finite element for the given entity.
      template<class Element>
      typename Traits::FiniteElementType find(const Element& e) const {
        return factory.make(e.geometry());
      }
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_GLOBAL_HH
