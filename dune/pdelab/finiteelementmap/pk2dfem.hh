// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_Pk2DFEM_HH
#define DUNE_PDELAB_Pk2DFEM_HH

#include<dune/common/exceptions.hh>

#include<dune/localfunctions/lagrange/pk2d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap

    template<typename GV, typename D, typename R, unsigned int k>
    class Pk2DLocalFiniteElementMap : 
	  public LocalFiniteElementMapInterface<LocalFiniteElementMapTraits< Dune::Pk2DLocalFiniteElement<D,R,k> >, 
											Pk2DLocalFiniteElementMap<GV,D,R,k> >,
      public Countable
    {
      typedef Dune::Pk2DLocalFiniteElement<D,R,k> FE;
      typedef typename GV::IndexSet IndexSet;
     
    public:
	  //! \brief export type of the signature
	  typedef LocalFiniteElementMapTraits<FE> Traits;  

	  //! \brief Use when Imp has a standard constructor
	  Pk2DLocalFiniteElementMap (const GV& gv_) : is(gv_.indexSet())
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
        n0 = is.subIndex(e,0,2);
        n1 = is.subIndex(e,1,2);
        n2 = is.subIndex(e,2,2);
        unsigned int j=0;
        if (n0>n1) j += 1;
        if (n0>n2) j += 2;
        if (n1>n2) j += 4;
		return variant[j];
	  }

	private:
      FE variant[8];
      const IndexSet& is;
    };
  }
}

#endif
