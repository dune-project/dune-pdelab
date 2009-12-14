// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MONOMFEM_HH
#define DUNE_PDELAB_MONOMFEM_HH

#include <dune/common/geometrytype.hh>

#include<dune/localfunctions/monom.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap
	template<class D, class R, int d, int p>
	class MonomLocalFiniteElementMap
	  : public SimpleLocalFiniteElementMap< Dune::MonomLocalFiniteElement<D,R,d,p> >
	{
    public:
      MonomLocalFiniteElementMap (Dune::GeometryType::BasicType basicType)
        : SimpleLocalFiniteElementMap< Dune::MonomLocalFiniteElement<D,R,d,p> >(Dune::MonomLocalFiniteElement<D,R,d,p>(basicType))
      {
      }
    };

  }
}

#endif //DUNE_PDELAB_MONOMFEM_HH
