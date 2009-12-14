// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_EDGER02DFEM_HH
#define DUNE_PDELAB_EDGER02DFEM_HH

#include<dune/localfunctions/edger02d.hh>

#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! Finite element map for EdgeR02D local finite elements
    //! \ingroup FiniteElementMap
	template<class D, class R>
	class EdgeR02DLocalFiniteElementMap
	  : public SimpleLocalFiniteElementMap< Dune::EdgeR02DLocalFiniteElement<D,R> >
	{};

  }
}

#endif // DUNE_PDELAB_EDGER02DFEM_HH
