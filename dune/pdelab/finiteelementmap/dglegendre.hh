// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_DGLEGENDRE_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_DGLEGENDRE_HH

#include "finiteelementmap.hh"
#include <dune/pdelab/finiteelement/dglegendre.hh>

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R, int k, int d>
    class DGLegendreLocalFiniteElementMap
      : public Dune::PDELab::SimpleLocalFiniteElementMap< Dune::DGLegendreLocalFiniteElement<D,R,k,d> >
    {

    public:

      bool fixedSize() const
      {
        return true;
      }

      bool hasDOFs (int codim) const
      {
        return codim == 0;
      }

      std::size_t size (GeometryType gt) const
      {
        if (gt == GeometryType(GeometryType::cube,d))
          return Dune::QkStuff::QkSize<k,d>::value;
        else
          return 0;
      }

      std::size_t maxLocalSize() const
      {
        return Dune::LegendreStuff::LegendreSize<k,d>::value;
      }

    };

  }
}

#endif