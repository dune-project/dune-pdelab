// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_Q1FEM_HH
#define DUNE_PDELAB_Q1FEM_HH

#include<dune/finiteelements/q1.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap
	template<class D, class R, int d>
	class Q1LocalFiniteElementMap
	  : public SimpleLocalFiniteElementMap< Dune::Q1LocalFiniteElement<D,R,d> >
	{};

  }
}

#endif
