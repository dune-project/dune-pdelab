// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_P1FEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_P1FEM_HH

#warning dune/pdelab/finiteelementmap/p1fem.hh and P1LocalFiniteElementMap are deprecated, please use dune/pdelab/finiteelementmap/pkfem.hh and PkLocalFiniteElementMap instead

#include <dune/common/deprecated.hh>
#include <dune/localfunctions/lagrange/p1.hh>
#include "finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R, int d>
    class DUNE_DEPRECATED_MSG("Please use PkLocalFiniteElementMap instead") P1LocalFiniteElementMap
      : public SimpleLocalFiniteElementMap< Dune::P1LocalFiniteElement<D,R,d> >
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
        return d+1;
      }
    };

  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_P1FEM_HH
