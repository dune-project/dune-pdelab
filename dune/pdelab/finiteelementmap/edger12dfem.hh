// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_EDGER12DFEM_HH
#define DUNE_PDELAB_EDGER12DFEM_HH

#include<dune/finiteelements/edger12d.hh>

#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! Finite element map for EdgeR12D local finite elements
    //! \ingroup FiniteElementMap
	template<class D, class R>
	class EdgeR12DLocalFiniteElementMap
	  : public SimpleLocalFiniteElementMap< Dune::EdgeR12DLocalFiniteElement<D,R> >
	{};

  }
}

#endif // DUNE_PDELAB_EDGER12DFEM_HH
