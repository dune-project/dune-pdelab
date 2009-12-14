// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_Q12DFEM_HH
#define DUNE_PDELAB_Q12DFEM_HH

#include<dune/localfunctions/q1.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap
	template<class D, class R>
	class Q12DLocalFiniteElementMap
	  : public SimpleLocalFiniteElementMap< Dune::Q1LocalFiniteElement<D,R,2> >
	{};

  }
}

#endif
