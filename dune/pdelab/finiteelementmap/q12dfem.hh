#ifndef DUNE_PDELAB_Q12DFEM_HH
#define DUNE_PDELAB_Q12DFEM_HH

#include<dune/finiteelements/q12d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	// wrap up element from local functions
	template<class D, class R>
	class Q12DLocalFiniteElementMap
	  : public SimpleLocalFiniteElementMap< Dune::Q12DLocalFiniteElement<D,R> >
	{};

  }
}

#endif
