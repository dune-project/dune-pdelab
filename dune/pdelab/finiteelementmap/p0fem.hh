// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_P0FEM_HH
#define DUNE_PDELAB_P0FEM_HH

#include <dune/geometry/type.hh>

#include<dune/localfunctions/lagrange/p0.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap
	template<class D, class R, int d>
	class P0LocalFiniteElementMap
	  : public SimpleLocalFiniteElementMap< Dune::P0LocalFiniteElement<D,R,d> >
	{
    public:

      P0LocalFiniteElementMap (const Dune::GeometryType& type)
        : SimpleLocalFiniteElementMap< Dune::P0LocalFiniteElement<D,R,d> >(Dune::P0LocalFiniteElement<D,R,d>(type))
        , _gt(type)
      {
      }

      bool fixedSize() const
      {
        return true;
      }

      std::size_t size(GeometryType gt) const
      {
        return gt == _gt ? 1 : 0;
      }

      std::size_t maxLocalSize() const
      {
        return 1;
      }

    private:
      const GeometryType _gt;

    };

  }
}

#endif
