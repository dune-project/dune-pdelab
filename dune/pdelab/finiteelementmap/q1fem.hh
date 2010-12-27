// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_Q1FEM_HH
#define DUNE_PDELAB_Q1FEM_HH

#include<dune/localfunctions/lagrange/q1.hh>
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
