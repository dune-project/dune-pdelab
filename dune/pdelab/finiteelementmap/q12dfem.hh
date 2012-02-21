// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_Q12DFEM_HH
#define DUNE_PDELAB_Q12DFEM_HH

#include<dune/localfunctions/lagrange/q1.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R>
    class Q12DLocalFiniteElementMap
      : public SimpleLocalFiniteElementMap< Dune::Q1LocalFiniteElement<D,R,2> >
    {

    public:

      bool fixedSize() const
      {
        return true;
      }

      std::size_t size(const GeometryType& gt) const
      {
        return gt.isVertex() ? 1 : 0;
      }

      std::size_t maxLocalSize() const
      {
        return 4;
      }

    };

  }
}

#endif
