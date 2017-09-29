// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_VARIABLEQKDGFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_VARIABLEQKDGFEM_HH

#include <memory>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualwrappers.hh>
#include <dune/common/array.hh>
#include "finiteelementmap.hh"
#include "qkdg.hh"

namespace Dune {
  namespace PDELab {

    namespace {
      template<class D, class R, int d, int p>
      struct InitVariableQkDGLocalFiniteElementMap
      {
        template<typename C>
        static void init(C & c)
        {
          typedef Dune::QkDGLocalFiniteElement<D,R,p,d> LFE;
          typedef typename C::value_type ptr;
          c[p] = ptr(new LocalFiniteElementVirtualImp<LFE>(LFE()));

          InitVariableQkDGLocalFiniteElementMap<D,R,d,p-1>::init(c);
        }
      };
      template<class D, class R, int d>
      struct InitVariableQkDGLocalFiniteElementMap<D,R,d,-1>
      {
        template<typename C>
        static void init(C & c) {}
      };
    }

    //! FiniteElementMap which provides QkDGLocalFiniteElement instances, depending on the local polynomial degree
    //! \ingroup FiniteElementMap
    template<class M, class D, class R, int d, int maxP=6>
    class VariableQkDGLocalFiniteElementMap
    {
      typedef typename FixedOrderLocalBasisTraits<
      typename QkDGLocalFiniteElement<D,R,0,d>::Traits::LocalBasisType::Traits,0>::Traits T;
      //! Type of finite element from local functions
      typedef LocalFiniteElementVirtualInterface<T> FiniteElementType;
    public:
      typedef FiniteElementMapTraits<FiniteElementType> Traits;

      //! The dimension of the finite elements returned by this map.
      static constexpr int dimension = d;

      VariableQkDGLocalFiniteElementMap (const M & m, unsigned int defaultP) :
        mapper_(m), polOrder_(mapper_.size(), defaultP), defaultP_(defaultP)
      {
        InitVariableQkDGLocalFiniteElementMap<D,R,d,maxP>::init(finiteElements_);
      }

      //! \brief get local basis functions for entity
      template<class EntityType>
      const typename Traits::FiniteElementType& find (const EntityType& e) const
      {
        return getFEM(getOrder(e));
      }

      //! \brief get local basis functions for a given polynomial order
      const typename Traits::FiniteElementType& getFEM (unsigned int p) const
      {
        return *(finiteElements_[p]);
      }

      //! \brief get local basis functions for the default order
      const typename Traits::FiniteElementType& getFEM () const
      {
        return *(finiteElements_[defaultP_]);
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
        unsigned int p = polOrder_[i];
        assert(p <= maxP);
        return p;
      }

      static constexpr bool fixedSize()
      {
        return false;
      }

      static constexpr bool hasDOFs(int codim)
      {
        return codim == 0;
      }

      std::size_t size(GeometryType gt) const
      {
        DUNE_THROW(Dune::Exception,"This should not be called!");
      }

      std::size_t maxLocalSize() const
      {
        return getFEM(maxP).localCoefficients().size();
      }

    private:
      const M & mapper_;
      std::vector<unsigned char> polOrder_;
      unsigned int defaultP_;
      std::array< std::shared_ptr<FiniteElementType>, maxP+1 > finiteElements_;
    };


  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_VARIABLEQKDGFEM_HH
