// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_VARIABLEMONOMFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_VARIABLEMONOMFEM_HH

#include <memory>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualwrappers.hh>
#include <dune/localfunctions/monomial.hh>
#include <dune/common/array.hh>
#include "finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    namespace {
      template<class D, class R, int d, int p>
      struct InitVariableMonomLocalFiniteElementMap
      {
        template<typename C>
        static void init(C & c, GeometryType gt)
        {
          typedef Dune::MonomialLocalFiniteElement<D,R,d,p> LFE;
          typedef typename C::value_type ptr;
          c[p] = ptr(new LocalFiniteElementVirtualImp<LFE>(LFE(gt)));

          InitVariableMonomLocalFiniteElementMap<D,R,d,p-1>::init(c,gt);
        }
      };
      template<class D, class R, int d>
      struct InitVariableMonomLocalFiniteElementMap<D,R,d,-1>
      {
        template<typename C>
          static void init(C &, GeometryType) {}
      };
    }

    //! FiniteElementMap which provides MonomialLocalFiniteElement instances, depending on the local polynomial degree
    //! \ingroup FiniteElementMap
    template<class M, class D, class R, int d, int maxP=6>
    class VariableMonomLocalFiniteElementMap
    {
      typedef typename FixedOrderLocalBasisTraits<
        typename MonomialLocalFiniteElement<D,R,d,0>::Traits::LocalBasisType::Traits,0>::Traits T;
      //! Type of finite element from local functions
      typedef LocalFiniteElementVirtualInterface<T> FiniteElementType;
    public:
      typedef FiniteElementMapTraits<FiniteElementType> Traits;

      //! The dimension of the finite elements returned by this map.
      static constexpr int dimension = d;

      /** Construct a VariableMonomLocalFiniteElementMap for GeometryType Dune::cube */
      VariableMonomLocalFiniteElementMap (const M & m, unsigned int defaultP) :
        gt_(Dune::GeometryType::cube,d), mapper_(m), polOrder_(mapper_.size(), defaultP), defaultP_(defaultP)
      {
        InitVariableMonomLocalFiniteElementMap<D,R,d,maxP>::init(finiteElements_, gt_);
      }

      /** Construct a VariableMonomLocalFiniteElementMap for a given GeometryType gt */
      VariableMonomLocalFiniteElementMap (const M & m, Dune::GeometryType gt, unsigned int defaultP) :
        gt_(gt), mapper_(m), polOrder_(mapper_.size(), defaultP), defaultP_(defaultP)
      {
        InitVariableMonomLocalFiniteElementMap<D,R,d,maxP>::init(finiteElements_, gt_);
      }

      //! \brief get local basis functions for entity
      template<class EntityType>
      const typename Traits::FiniteElementType& find (const EntityType& e) const
      {
        if (e.type() != gt_)
          DUNE_THROW(InvalidGeometryType,"Unsupported geometry type: Support only " << gt_ << ", but got " << e.type());
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
        unsigned int i = mapper_.index(e);
        polOrder_[i] = p;
      }

      template<class EntityType>
      unsigned int getOrder (const EntityType& e) const
      {
        unsigned int i = mapper_.index(e);
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
        DUNE_THROW(VariableElementSize,"VariableMonomLocalFiniteElementMap can contain elements of variable order.");
      }

      std::size_t maxLocalSize() const
      {
        DUNE_THROW(VariableElementSize,"VariableMonomLocalFiniteElementMap can contain elements of variable order.");
      }

    private:
      const Dune::GeometryType gt_;
      const M & mapper_;
      std::vector<unsigned char> polOrder_;
      unsigned int defaultP_;
      std::array< std::shared_ptr<FiniteElementType>, maxP+1 > finiteElements_;
    };



  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_VARIABLEMONOMFEM_HH
