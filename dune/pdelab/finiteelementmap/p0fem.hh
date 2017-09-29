// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_P0FEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_P0FEM_HH

#include <dune/geometry/type.hh>

#include<dune/localfunctions/lagrange/p0.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<class D, class R, int d>
    class P0LocalFiniteElementMap
      : public SimpleLocalFiniteElementMap<Dune::P0LocalFiniteElement<D,R,d>,d>
    {
    public:

      P0LocalFiniteElementMap (const Dune::GeometryType& type)
        : SimpleLocalFiniteElementMap<Dune::P0LocalFiniteElement<D,R,d>,d>(Dune::P0LocalFiniteElement<D,R,d>(type))
        , _gt(type)
      {
      }

      static constexpr bool fixedSize()
      {
        return true;
      }

      static constexpr bool hasDOFs(int codim)
      {
        return codim == 0;
      }

      std::size_t size(GeometryType gt) const
      {
        return gt == _gt ? 1 : 0;
      }

      static constexpr std::size_t maxLocalSize()
      {
        return 1;
      }

    private:
      const GeometryType _gt;

    };

  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_P0FEM_HH
