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
	{};

  }
}

#endif
