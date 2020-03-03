// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_VECTORHELPERS_HH
#define DUNE_PDELAB_BACKEND_ISTL_VECTORHELPERS_HH

#include <dune/common/typetraits.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dune/pdelab/backend/istl/tags.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/finiteelementmap/utility.hh>
#include <dune/pdelab/ordering/orderingbase.hh>

#include <dune/pdelab/backend/istl/blocking.hh>
#include <dune/pdelab/backend/istl/blockinggfs.hh>

namespace Dune {
  namespace PDELab {
    namespace ISTL {

      template<typename CI, typename Block>
      typename Block::field_type&
      access_vector_element(tags::field_vector_1, Block& b, const CI& ci, int i)
      {
        // usually we are at the end of the multi-index (-1),
        // but we might be in a PowerFunctionSpace of size 1,
        // then we are at the lowest multi-index component (0)
        assert(i == -1 || i == 0);
        return b[0];
      }

      template<typename CI, typename Block>
      typename Block::field_type&
      access_vector_element(tags::field_vector_n, Block& b, const CI& ci, int i)
      {
        assert(i == 0);
        return b[ci[0]];
      }

      template<typename CI, typename Block>
      typename Block::field_type&
      access_vector_element(tags::block_vector, Block& b, const CI& ci, int i)
      {
        return access_vector_element(container_tag(b[ci[i]]),b[ci[i]],ci,i-1);
      }


      template<typename CI, typename Block>
      const typename Block::field_type&
      access_vector_element(tags::field_vector_1, const Block& b, const CI& ci, int i)
      {
        // usually we are at the end of the multi-index (-1),
        // but we might be in a PowerFunctionSpace of size 1,
        // then we are at the lowest multi-index component (0)
        assert(i == -1 || i == 0);
        return b[0];
      }

      template<typename CI, typename Block>
      const typename Block::field_type&
      access_vector_element(tags::field_vector_n, const Block& b, const CI& ci, int i)
      {
        assert(i == 0);
        return b[ci[0]];
      }

      template<typename CI, typename Block>
      const typename Block::field_type&
      access_vector_element(tags::block_vector, const Block& b, const CI& ci, int i)
      {
        return access_vector_element(container_tag(b[ci[i]]),b[ci[i]],ci,i-1);
      }

      // ********************************************************************************
      // size management, not used anymore, now done via dune-functions...
      // ********************************************************************************
#if 0
      template<typename Vector>
      void resize_vector(tags::block_vector, Vector& v, std::size_t size, bool copy_values)
      {
        v.resize(size);
      }

      template<typename Vector>
      void resize_vector(tags::field_vector, Vector& v, std::size_t size, bool copy_values)
      {
      }

      template<typename DI, typename CI, typename Container>
      void allocate_vector(tags::field_vector, const OrderingBase<DI,CI>& ordering, Container& c)
      {
      }

      template<typename DI, typename CI, typename Container>
      void allocate_vector(tags::block_vector, const OrderingBase<DI,CI>& ordering, Container& c)
      {
        for (std::size_t i = 0; i < ordering.childOrderingCount(); ++i)
          {
            if (ordering.containerBlocked())
              {
                resize_vector(container_tag(c[i]),c[i],ordering.childOrdering(i).blockCount(),false);
                allocate_vector(container_tag(c[i]),ordering.childOrdering(i),c[i]);
              }
            else
              allocate_vector(container_tag(c),ordering.childOrdering(i),c);
          }
      }

      template<typename Ordering, typename Container>
      void dispatch_vector_allocation(const Ordering& ordering, Container& c, HierarchicContainerAllocationTag tag)
      {
        allocate_vector(container_tag(c),ordering,c);
      }

      template<typename Ordering, typename Container>
      void dispatch_vector_allocation(const Ordering& ordering, Container& c, FlatContainerAllocationTag tag)
      {
        resize_vector(container_tag(c),c,ordering.blockCount(),false);
      }
#endif

      // ********************************************************************************
      // map blocking info to ISTL vector type
      // ********************************************************************************
      template<typename EntryType, typename Blocking>
      struct VectorType;

      template<typename E, std::size_t SZ>
      struct VectorType<E,Dune::PDELab::Blocking::tag::leafBlocked<SZ>>
      {
        using Type = Dune::BlockVector<FieldVector<E,SZ>>;
      };

      template<typename E>
      struct VectorType<E,Dune::PDELab::Blocking::tag::flat>
      {
        using Type = Dune::BlockVector<FieldVector<E,1>>;
        // using Type = Dune::BlockVector<E>;
      };

      template<typename E, typename First>
      struct VectorType<E,Dune::PDELab::Blocking::tag::blocked<First>>
      {
        using Type = Dune::BlockVector<typename VectorType<E,First>::Type>;
      };

      template<typename E, typename First, typename... Siblings>
      struct VectorType<E,Dune::PDELab::Blocking::tag::blocked<First, Siblings...>>
      {
        using Type =
          std::conditional_t<
            std::conjunction_v<std::is_same<First, Siblings>...>,
            Dune::BlockVector<typename VectorType<E,First>::Type>,
            Dune::MultiTypeBlockVector<
              typename VectorType<E,First>::Type,
              typename VectorType<E,Siblings>::Type...>>;
      };

    } // namespace ISTL
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_VECTORHELPERS_HH
