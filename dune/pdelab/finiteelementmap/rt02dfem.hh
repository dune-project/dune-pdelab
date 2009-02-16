// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_RT02DFEM_HH
#define DUNE_PDELAB_RT02DFEM_HH

#include<dune/finiteelements/rt02d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap

    template<typename GV, typename D, typename R>
    class RT02DLocalFiniteElementMap : 
	  public LocalFiniteElementMapInterface<LocalFiniteElementMapTraits< Dune::RT02DLocalFiniteElement<D,R> >, 
											RT02DLocalFiniteElementMap<GV,D,R> >,
      public Countable
    {
      typedef Dune::RT02DLocalFiniteElement<D,R> FE;
      typedef typename GV::IndexSet IndexSet;
     
    public:
	  //! \brief export type of the signature
	  typedef LocalFiniteElementMapTraits<FE> Traits;  

	  //! \brief Use when Imp has a standard constructor
	  RT02DLocalFiniteElementMap (const GV& gv_) : is(gv_.indexSet())
	  {
        // create all variants 
        for (int i=0; i<8; i++)
          variant[i] = FE(i);
      }

	  //! \brief get local basis functions for entity
	  template<class EntityType>
	  const typename Traits::LocalFiniteElementType& find (const EntityType& e) const
	  {
        unsigned int n0,n1,n2;
        n0 = is.template subIndex<2>(e,0);
        n1 = is.template subIndex<2>(e,1);
        n2 = is.template subIndex<2>(e,2);
        unsigned int j=0;
        if (n1>n2) j += 1;
        if (n0>n2) j += 2;
        if (n0>n1) j += 4;
		return variant[j];
	  }

	private:
      FE variant[8];
      const IndexSet& is;
    };

  }
}

#endif
