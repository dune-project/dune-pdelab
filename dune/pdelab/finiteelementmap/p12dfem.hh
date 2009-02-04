#ifndef DUNE_PDELAB_P12DFEM_HH
#define DUNE_PDELAB_P12DFEM_HH

#include<dune/finiteelements/p12d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	// wrap up element from local functions
	template<class D, class R>
	class P12DLocalFiniteElementMap
	  : public SimpleLocalFiniteElementMap< Dune::P12DLocalFiniteElement<D,R> >
	{};

  }
}

#endif
