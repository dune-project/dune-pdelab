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
      /**
       * \deprecated Use the constructor that takes a Dune::GeometryType
       *             argument.
       */
      P0LocalFiniteElementMap (Dune::GeometryType::BasicType basicType) DUNE_DEPRECATED
        : SimpleLocalFiniteElementMap< Dune::P0LocalFiniteElement<D,R,d> >(Dune::P0LocalFiniteElement<D,R,d>(basicType))
      {
      }
      P0LocalFiniteElementMap (const Dune::GeometryType& type)
        : SimpleLocalFiniteElementMap< Dune::P0LocalFiniteElement<D,R,d> >(Dune::P0LocalFiniteElement<D,R,d>(type))
      {
      }

    };

  }
}

#endif
