// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_Q22DFEM_HH
#define DUNE_PDELAB_Q22DFEM_HH

#include<dune/localfunctions/q22d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap
	template<class D, class R>
	class Q22DLocalFiniteElementMap
	  : public SimpleLocalFiniteElementMap< Dune::Q22DLocalFiniteElement<D,R> >
	{};

  }
}

#endif
