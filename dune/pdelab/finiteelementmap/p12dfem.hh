// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_P12DFEM_HH
#define DUNE_PDELAB_P12DFEM_HH

#include<dune/localfunctions/lagrange/p1.hh>
#include"../common/geometrywrapper.hh"
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R>
    class P12DLocalFiniteElementMap
      : public SimpleLocalFiniteElementMap< Dune::P1LocalFiniteElement<D,R,2> >
    {

    public:

      bool fixedSize() const
      {
        return true;
      }

      std::size_t size(GeometryType gt) const
      {
        if (gt.isVertex())
          return 1;
        return 0;
      }

      std::size_t maxLocalSize() const
      {
        return 3;
      }

    };

  }
}

#endif
