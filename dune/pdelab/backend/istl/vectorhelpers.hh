// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_VECTORHELPERS_HH
#define DUNE_PDELAB_BACKEND_ISTL_VECTORHELPERS_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <dune/common/typetraits.hh>

#include <dune/istl/bvector.hh>

#include <dune/pdelab/backend/istl/tags.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/finiteelementmap/utility.hh>
#include <dune/pdelab/ordering/orderingbase.hh>

namespace Dune {
  namespace PDELab {

    // Recursive accessors for vector entries using tag dispatch

#ifndef DOXYGEN // All of the following functions are mere implementation details

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


      template<typename Vector>
      void resize_vector(tags::block_vector, Vector& v, std::size_t size, bool copy_values)
      {
        v.resize(size,copy_values);
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


      // ********************************************************************************
      // TMPs for deducing ISTL block structure from GFS backends
      // ********************************************************************************

      // tag dispatch switch on GFS tag for per-node functor - general version
      template<typename E,typename Node, typename Tag, bool isLeafTag = std::is_base_of<LeafGridFunctionSpaceTag,Tag>::value >
      struct vector_descriptor_helper
      {
        // export backend type, as the actual TMP is in the parent reduction functor
        typedef Node type;
      };

      // descriptor for backends of leaf spaces collecting various information about
      // possible blocking structures
      template<typename E, typename GFS>
      struct leaf_vector_descriptor
      {

        using Backend = typename GFS::Traits::Backend;
        using FEM = typename GFS::Traits::FiniteElementMap;

        static_assert(Backend::Traits::block_type != Blocking::bcrs,
                      "Dynamically blocked leaf spaces are not supported by this backend.");

        // flag for sibling reduction - always true in the leaf case
        static const bool support_no_blocking = true;

        // flag indicating whether the associated vector type supports cascading
        // the static blocking further up the tree (i.e. create larger static blocks
        // at the parent node level. Due to ISTL limitations, this only works once in
        // the hierarchy, so we only support cascading if we don't already do static
        // blocking at the current level.
        static const bool support_cascaded_blocking =
          Backend::Traits::block_type == Blocking::none; // FIXME

        // The cumulative block size is used by the algorithm to calculate total block
        // size over several children for cascaded blocking. We try to extract this size
        // from the finite element map (that works if the number of DOFs per entity is
        // an identical constant for all types of entities that have DOFs attached), and
        // if we fail we fall back to 1, which always works.
        static const std::size_t detected_cumulative_block_size = Dune::Std::detected_or_t<
          std::integral_constant<std::size_t,0>,
          FiniteElementMapBlockSize,
          FEM
          >::value;

        // Has the user selected an explicit block size? If yes, give priority to that.
        static const std::size_t cumulative_block_size = Backend::Traits::block_size > 0
          ? Backend::Traits::block_size
          : detected_cumulative_block_size;

        static constexpr bool have_valid_block_size = cumulative_block_size > 0;

        static_assert(
          Backend::Traits::block_size == 0 or (detected_cumulative_block_size % cumulative_block_size) == 0,
          "The vector block size you specified is not compatible with the finite element map"
          );

        static_assert(
          Backend::Traits::block_type != Blocking::fixed or have_valid_block_size,
          "You requested static blocking, but we cannot extract a valid block size from the finite element map. Please specify the block size with the second template parameter of the vector backend."
          );

        // The static block size of the associated vector
        static const std::size_t block_size =
          Backend::Traits::block_type == Blocking::fixed ? cumulative_block_size : 1;

        // The element type for the vector.
        typedef E element_type;

        // The ISTL vector type associated with the current subtree.
        typedef Dune::BlockVector<FieldVector<E,block_size> > vector_type;

      };

      // Tag dispatch for leaf spaces - extract leaf descriptor.
      template<typename E, typename Node, typename Tag>
        struct vector_descriptor_helper<E,Node,Tag, /* is LeafTag */ true>
      {
        typedef leaf_vector_descriptor<E,Node> type;
      };

      // the actual functor
      template<typename E>
      struct extract_vector_descriptor
      {

        template<typename Node, typename TreePath>
        struct doVisit
        {
          // visit all nodes
          static const bool value = true;
        };

        template<typename Node, typename TreePath>
        struct visit
        {
          // forward to actual implementation via tag dispatch
          typedef typename vector_descriptor_helper<E,Node,TypeTree::ImplementationTag<Node>>::type type;
        };

      };

      // Descriptor for combining sibling nodes in the tree
      template<typename Sibling, typename Child>
      struct cascading_vector_descriptor
      {

        // We only support cascaded blocking if all children support it
        static const bool support_cascaded_blocking =
          Sibling::support_cascaded_blocking &&
          Child::support_cascaded_blocking;

        // ISTL requires a single, globally constant blocking structure
        // for its containers, so we make sure the siblings don't disagree
        // on it.
        static const bool support_no_blocking =
          (Sibling::support_no_blocking &&
           std::is_same<
            typename Sibling::vector_type,
            typename Child::vector_type
           >::value);

        static constexpr bool have_valid_block_size =
          Sibling::have_valid_block_size and Child::have_valid_block_size;

        // block size
        static const std::size_t block_size =
          support_no_blocking ? Sibling::block_size : 1;

        // The element type for the vector.
        typedef typename Sibling::element_type element_type;

        // Accumulate total block size of all siblings
        static const std::size_t cumulative_block_size =
          Sibling::cumulative_block_size + Child::cumulative_block_size;

        // The ISTL vector type associated with the current subtree.
        typedef Dune::BlockVector<FieldVector<element_type,block_size> > vector_type;

      };


      // Switch that turns off standard reduction for the first child of a node.
      // Default case: do the standard reduction.
      template<typename D1, typename D2>
      struct initial_reduction_switch
      {
        typedef cascading_vector_descriptor<D1,D2> type;
      };

      // specialization for first child
      template<typename D2>
      struct initial_reduction_switch<void,D2>
      {
        typedef D2 type;
      };

      // sibling reduction functor
      struct combine_vector_descriptor_siblings
      {

        template<typename D1, typename D2>
        struct reduce
          : public initial_reduction_switch<D1,D2>
        {};

      };

      // Data part of child -> parent reduction descriptor
      template<typename Child, typename GFS>
      struct parent_child_vector_descriptor_data
      {

        using Backend = typename GFS::Traits::Backend;

        static constexpr bool have_valid_block_size = Child::have_valid_block_size;

        // If all our have a common blocking structure, we can just
        // concatenate them without doing any blocking
        static const bool support_no_blocking =
          Child::support_no_blocking;

        // We support cascaded blocking if neither we nor any of our
        // children are blocked yet.
        static const bool support_cascaded_blocking =
          Child::support_cascaded_blocking &&
          Backend::Traits::block_type == Blocking::none;

        // It is not allowed to specify a block size on an interior node
        static_assert(
          Backend::Traits::block_size == 0,
          "You cannot specify a block size on interior nodes of the function space tree."
          );

        // Throw an assertion if the user requests static blocking at this level,
        // but we cannot support it.
        static_assert((Backend::Traits::block_type != Blocking::fixed) ||
                      Child::support_cascaded_blocking,
                      "invalid blocking structure.");

        static_assert(
          Backend::Traits::block_type != Blocking::fixed or have_valid_block_size,
          "You requested static blocking, but at least one leaf space has a finite element that does not support automatic block size extraction. Please specify the block size with the second template parameter of that space's vector backend."
          );

        // If we block statically, we create bigger blocks, otherwise the
        // block size doesn't change.
        static const std::size_t block_size =
          Backend::Traits::block_type == Blocking::fixed
          ? Child::cumulative_block_size
          : Child::block_size;

        // Just forward this...
        static const std::size_t cumulative_block_size =
          Child::cumulative_block_size;

        // The element type for the vector.
        typedef typename Child::element_type element_type;

        // The ISTL vector type associated with our subtrees.
        typedef typename Child::vector_type child_vector_type;

      };

      // dispatch switch on blocking type - prototype
      template<typename Data, Blocking>
      struct parent_child_vector_descriptor;

      // dispatch switch on blocking type - no blocking case
      template<typename Data>
      struct parent_child_vector_descriptor<
        Data,
        Blocking::none
        >
        : public Data
      {
        static_assert(Data::support_no_blocking,
                      "Cannot combine incompatible child block structures without static blocking. "
                      "Did you want to apply static blocking at this level?");

        // Just forward the child vector type
        typedef typename Data::child_vector_type vector_type;
      };

      // dispatch switch on blocking type - dynamic blocking case
      template<typename Data>
      struct parent_child_vector_descriptor<
        Data,
        Blocking::bcrs
        >
        : public Data
      {
        static_assert(Data::support_no_blocking,
                      "Incompatible child block structures detected, cannot perform dynamic blocking. "
                      "Did you want to apply static blocking at this level?");

        // Wrap the child vector type in another BlockVector
        typedef Dune::BlockVector<typename Data::child_vector_type> vector_type;
      };

      // dispatch switch on blocking type - static blocking case
      template<typename Data>
      struct parent_child_vector_descriptor<
        Data,
        Blocking::fixed
        >
        : public Data
      {
        // build new block vector with large field block size
        typedef Dune::BlockVector<
          FieldVector<
            typename Data::element_type,
            Data::block_size
            >
          > vector_type;
      };

      // Child - parent reduction functor
      struct combine_vector_descriptor_parent
      {

        template<typename Child, typename GFS>
        struct reduce
        {

          struct type
            : public parent_child_vector_descriptor<parent_child_vector_descriptor_data<
                                                      Child,
                                                      GFS>,
                                                    GFS::Traits::Backend::Traits::block_type
                                                    >
          {};
        };

      };

      // policy describing the GFS tree -> ISTL vector reduction
      template<typename E>
      struct vector_creation_policy
        : public TypeTree::TypeAccumulationPolicy<extract_vector_descriptor<E>,
                                                  combine_vector_descriptor_siblings,
                                                  void,
                                                  combine_vector_descriptor_parent,
                                                  TypeTree::bottom_up_reduction>
      {};

    } // namespace ISTL


#endif // DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_VECTORHELPERS_HH
