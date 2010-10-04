// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_GLOBALP1_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_GLOBALP1_HH

#include <dune/localfunctions/lagrange/p1.hh>

#include <dune/pdelab/finiteelement/localtoglobaladaptors.hh>
#include <dune/pdelab/finiteelementmap/global.hh>

namespace Dune {
  namespace PDELab {

    template<class Geometry, class RangeField>
    class P1FiniteElementMap :
      public GeometryFiniteElementMap
        <ScalarLocalToGlobalFiniteElementAdaptorFactory<
           P1LocalFiniteElement<typename Geometry::ctype,
                                RangeField,
                                Geometry::mydimension>,
           Geometry> >
    {
      typedef P1LocalFiniteElement<typename Geometry::ctype,
                                   RangeField,
                                   Geometry::mydimension> LocalFE;
      typedef ScalarLocalToGlobalFiniteElementAdaptorFactory<LocalFE, Geometry>
        Factory;
      typedef GeometryFiniteElementMap<Factory> Base;

      static Factory& getFactory() {
        static const LocalFE localFE;
        static Factory factory(localFE);
        return factory;
      }

    public:
      P1FiniteElementMap() : Base(getFactory()) { }

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_GLOBALP1_HH
