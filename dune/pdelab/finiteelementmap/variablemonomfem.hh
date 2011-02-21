// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_VARIABLEMONOMFEM_HH
#define DUNE_PDELAB_VARIABLEMONOMFEM_HH

#include <dune/common/geometrytype.hh>

#include <dune/localfunctions/common/virtualwrappers.hh>
#include <dune/localfunctions/monom.hh>
#include <dune/common/array.hh>
#include <dune/common/shared_ptr.hh>
#include "finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    namespace {
      template<class D, class R, int d, int p>
      struct InitVariableMonomLocalFiniteElementMap
      {
        template<typename C>
        static void init(C & c)
        {
          typedef Dune::MonomLocalFiniteElement<D,R,d,p> LFE;
          typedef typename C::value_type ptr;
          static GeometryType cube(d);
          c[p] = ptr(new LocalFiniteElementVirtualImp<LFE>(LFE(cube)));

          InitVariableMonomLocalFiniteElementMap<D,R,d,p-1>::init(c);
        }
      };
      template<class D, class R, int d>
      struct InitVariableMonomLocalFiniteElementMap<D,R,d,-1>
      {
        template<typename C>
        static void init(C & c) {}
      };
    }

	//! FiniteElementMap which provides MonomLocalFiniteElement instances, depending on the local polynomial degree
    //! \ingroup FiniteElementMap
	template<class M, class D, class R, int d, int maxP>
	class VariableMonomLocalFiniteElementMap
	{
      typedef typename FixedOrderLocalBasisTraits<
        typename MonomLocalFiniteElement<D,R,d,0>::Traits::LocalBasisType::Traits,0>::Traits T;
      //! Type of finite element from local functions
      typedef LocalFiniteElementVirtualInterface<T> FiniteElementType;
    public:
      typedef FiniteElementMapTraits<FiniteElementType> Traits;

      VariableMonomLocalFiniteElementMap (const M & m, unsigned int defaultP) :
        mapper_(m), polOrder_(mapper_.size(), defaultP)
      {
        InitVariableMonomLocalFiniteElementMap<D,R,d,maxP>::init(finiteElements_);
      }

	  //! \brief get local basis functions for entity
	  template<class EntityType>
	  const typename Traits::FiniteElementType& find (const EntityType& e) const
	  {
        unsigned int p = getOrder(e);
        assert(p <= maxP);
		return *(finiteElements_[p]);
	  }

	  template<class EntityType>
	  void setOrder (const EntityType& e, unsigned int p)
	  {
        assert(p <= maxP);
        unsigned int i = mapper_.map(e);
        polOrder_[i] = p;
	  }

	  template<class EntityType>
	  unsigned int getOrder (const EntityType& e) const
	  {
        unsigned int i = mapper_.map(e);
        return polOrder_[i];
	  }

	private:
      const M & mapper_;
      std::vector<unsigned char> polOrder_;
      Dune::array< Dune::shared_ptr<FiniteElementType>, maxP+1 > finiteElements_;
    };



  }
}

#endif //DUNE_PDELAB_VARIABLEMONOMFEM_HH
