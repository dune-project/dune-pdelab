// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_ORDERING_INTERLEAVEDORDERING_HH
#define DUNE_PDELAB_ORDERING_INTERLEAVEDORDERING_HH

#include <string>

#include <dune/pdelab/common/exceptions.hh>
#include <dune/pdelab/common/typetree/compositenodemacros.hh>
#include <dune/pdelab/common/typetree/powernode.hh>

#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/orderingbase.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    namespace interleaved_ordering {

      //! Interface for merging index spaces
      template<typename DI, typename GDI, typename CI, typename Node>
      class Base
        : public OrderingBase<DI,GDI,CI>
      {

      public:

        typedef typename OrderingBase<DI,GDI,CI>::Traits Traits;

        typedef InterleavedOrderingTag OrderingTag;

        static const bool consume_tree_index = true;

        //! Construct ordering object
        /**
         * In general, an ordering object is not properly setup after
         * construction.  This must be done by a seperate call to update()
         * after all the children have been properly set up.
         */
        Base(Node& node, bool container_blocked, const OrderingTag& ordering_tag)
          : OrderingBase<DI,GDI,CI>(node,container_blocked,ordering_tag.offsets(),nullptr)
        {
          // This check looks a little weird, but there is always one offset more than
          // there are blocks (the first offsets is 0, and the last one is the "offset
          // beyond the end" to encode the size of the final child).
          if (node.CHILDREN != ordering_tag.offsets().size() - 1)
            DUNE_THROW(OrderingStructureError,
                       "Invalid block structure for InterleavedOrdering: "
                       << node.CHILDREN << " children, but "
                       << (ordering_tag.offsets().size() - 1) << " block sizes.");
        }

        template<typename ItIn, typename ItOut>
        void map_lfs_indices(const ItIn begin, const ItIn end, ItOut out) const
        {
          typedef typename Traits::SizeType size_type;
          if (this->_container_blocked)
            {
              for (ItIn in = begin; in != end; ++in, ++out)
                {
                  size_type child_index = in->treeIndex().back();
                  size_type child_block_offset = this->_child_block_merge_offsets[child_index];
                  size_type child_block_size = this->_child_block_merge_offsets[child_index + 1] - child_block_offset;
                  size_type index = out->back();
                  size_type block_index = index / child_block_size;
                  size_type offset = index % child_block_size;
                  out->back() = child_block_offset + offset;
                  out->push_back(block_index);
                }
            }
          else
            {
              for (ItIn in = begin; in != end; ++in, ++out)
                {
                  size_type child_index = in->treeIndex().back();
                  size_type child_block_offset = this->_child_block_merge_offsets[child_index];
                  size_type child_block_size = this->_child_block_merge_offsets[child_index + 1] - child_block_offset;
                  size_type block_size = this->_child_block_merge_offsets.back();
                  size_type index = out->back();
                  size_type block_index = index / child_block_size;
                  size_type offset = index % child_block_size;
                  out->back() = block_index * block_size + child_block_offset + offset;
                }
            }
        }

        template<typename CIOutIterator, typename DIOutIterator = DummyDOFIndexIterator>
        typename Traits::SizeType
        extract_entity_indices(const typename Traits::DOFIndex::EntityIndex& ei,
                               typename Traits::SizeType child_index,
                               CIOutIterator ci_out, const CIOutIterator ci_end) const
        {
          typedef typename Traits::SizeType size_type;
          if (this->_container_blocked)
            {
              for (; ci_out != ci_end; ++ci_out)
                {
                  size_type child_block_offset = this->_child_block_merge_offsets[child_index];
                  size_type child_block_size = this->_child_block_merge_offsets[child_index + 1] - child_block_offset;
                  size_type index = ci_out->back();
                  size_type block_index = index / child_block_size;
                  size_type offset =index % child_block_size;
                  ci_out->back() = child_block_offset + offset;
                  ci_out->push_back(block_index);
                }
            }
          else
            {
              for (; ci_out != ci_end; ++ci_out)
                {
                  size_type child_block_offset = this->_child_block_merge_offsets[child_index];
                  size_type child_block_size = this->_child_block_merge_offsets[child_index + 1] - child_block_offset;
                  size_type block_size = this->_child_block_merge_offsets.back();
                  size_type index = ci_out->back();
                  size_type block_index = index / child_block_size;
                  size_type offset =index % child_block_size;
                  ci_out->back() = block_index * block_size + child_block_offset + offset;
                }
            }

          // The return value is not used for non-leaf orderings.
          return 0;
        }

      };

    } // namespace interleaved_ordering

    //! Interface for merging index spaces
    template<typename DI, typename GDI, typename CI, typename Child, std::size_t k>
    class PowerInterleavedOrdering
      : public TypeTree::PowerNode<Child, k>
      , public interleaved_ordering::Base<DI,
                                          GDI,
                                          CI,
                                          PowerInterleavedOrdering<DI,GDI,CI,Child,k>
                                          >
    {
      typedef TypeTree::PowerNode<Child, k> Node;

      typedef interleaved_ordering::Base<DI,
                                         GDI,
                                         CI,
                                         PowerInterleavedOrdering<DI,GDI,CI,Child,k>
                                         > Base;

    public:

      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update() after
       * all the children have been properly set up.
       *
       * \note This constructor must be present for ordering objects not at
       *       the leaf of the tree.
       */
      PowerInterleavedOrdering(bool container_blocked, const InterleavedOrderingTag& ordering_tag, const typename Node::NodeStorage& children)
        : Node(children)
        , Base(*this,container_blocked,ordering_tag)
      {}

      void update()
      {
        for (std::size_t i = 0; i < Node::CHILDREN; ++i)
          {
            this->child(i).update();
          }
        Base::update();
      }

      std::string name() const { return "PowerInterleavedOrdering"; }
    };


    template<typename GFS, typename Transformation>
    struct power_gfs_to_ordering_descriptor<GFS,Transformation,InterleavedOrderingTag>
    {

      static const bool recursive = true;

      template<typename TC>
      struct result
      {

        typedef PowerInterleavedOrdering<
          typename Transformation::DOFIndex,
          typename TC::Traits::GlobalDOFIndex,
          typename Transformation::ContainerIndex,
          TC,
          GFS::CHILDREN
          > type;

        typedef shared_ptr<type> storage_type;

      };

      template<typename TC>
      static typename result<TC>::type transform(const GFS& gfs, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
      {
        return typename result<TC>::type(gfs.backend().blocked(),gfs.orderingTag(),children);
      }

      template<typename TC>
      static typename result<TC>::storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
      {
        return make_shared<typename result<TC>::type>(gfs->backend().blocked(),gfs->orderingTag(),children);
      }

    };

    //! Interface for merging index spaces
    template<typename DI, typename GDI, typename CI, DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeInterleavedOrdering :
      public DUNE_TYPETREE_COMPOSITENODE_BASETYPE,
      public interleaved_ordering::Base<DI,
                                        GDI,
                                        CI,
                                        CompositeInterleavedOrdering<
                                          DI,
                                          GDI,
                                          CI,
                                          DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES
                                          >
                                        >
    {
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE Node;

      typedef interleaved_ordering::Base<
        DI,
        GDI,
        CI,
        CompositeInterleavedOrdering<
          DI,
          GDI,
          CI,
          DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES
          >
        > Base;

    public:
      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update() after
       * all the children have been properly set up.
       *
       * \note This constructor must be present for ordering objects not at
       *       the leaf of the tree.
       */
      CompositeInterleavedOrdering(bool backend_blocked, const InterleavedOrderingTag& ordering_tag, DUNE_TYPETREE_COMPOSITENODE_STORAGE_CONSTRUCTOR_SIGNATURE)
        : Node(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES)
        , Base(*this,backend_blocked,ordering_tag)
      { }

      std::string name() const { return "CompositeInterleavedOrdering"; }

      void update()
      {
        TypeTree::applyToTree(*this,ordering::update_direct_children());
        Base::update();
      }
    };

#if HAVE_VARIADIC_TEMPLATES

    template<typename GFS, typename Transformation>
    struct composite_gfs_to_ordering_descriptor<GFS,Transformation,InterleavedOrderingTag>
    {

      static const bool recursive = true;

      template<typename... TC>
      struct result
      {

        typedef CompositeInterleavedOrdering<
          typename Transformation::DOFIndex,
          typename extract_first_child<TC...>::type::Traits::GlobalDOFIndex,
          typename Transformation::ContainerIndex,
          TC...
          > type;

        typedef shared_ptr<type> storage_type;

      };

      template<typename... TC>
      static typename result<TC...>::type transform(const GFS& gfs, const Transformation& t, shared_ptr<TC>... children)
      {
        return typename result<TC...>::type(gfs.backend().blocked(),gfs.orderingTag(),children...);
      }

      template<typename... TC>
      static typename result<TC...>::storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t, shared_ptr<TC>... children)
      {
        return make_shared<typename result<TC...>::type>(gfs->backend().blocked(),gfs.orderingTag(),children...);
      }

    };

#else // HAVE_VARIADIC_TEMPLATES

    //! Node transformation descriptor for CompositeGridFunctionSpace -> LexicographicOrdering (without variadic templates).
    template<typename GFS, typename Transformation>
    struct composite_gfs_to_ordering_descriptor<GFS,Transformation,InterleavedOrderingTag>
    {

      static const bool recursive = true;

      template<typename TC0,
               typename TC1,
               typename TC2,
               typename TC3,
               typename TC4,
               typename TC5,
               typename TC6,
               typename TC7,
               typename TC8,
               typename TC9>
      struct result
      {
        // TODO: FIXME - this has not been changed to new interface yet!
        typedef CompositeInterleavedOrdering<typename Transformation::GridFunctionSpace::Traits::SizeType,
                                             TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9> type;
        typedef shared_ptr<type> storage_type;
      };

      template<typename TC0,
               typename TC1,
               typename TC2,
               typename TC3,
               typename TC4,
               typename TC5,
               typename TC6,
               typename TC7,
               typename TC8,
               typename TC9>
      static typename result<TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9>::type
      transform(const GFSNode& gfs,
                const Transformation& t,
                shared_ptr<TC0> c0,
                shared_ptr<TC1> c1,
                shared_ptr<TC2> c2,
                shared_ptr<TC3> c3,
                shared_ptr<TC4> c4,
                shared_ptr<TC5> c5,
                shared_ptr<TC6> c6,
                shared_ptr<TC7> c7,
                shared_ptr<TC8> c8,
                shared_ptr<TC9> c9)
      {
        return typename result<TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9>::type(gfs.backend().blocked(),gfs.orderingTag(),c0,c1,c2,c3,c4,c5,c6,c7,c8,c9);
      }

      template<typename TC0,
               typename TC1,
               typename TC2,
               typename TC3,
               typename TC4,
               typename TC5,
               typename TC6,
               typename TC7,
               typename TC8,
               typename TC9>
      static typename result<TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9>::storage_type
      transform_storage(shared_ptr<const GFSNode> gfs,
                        const Transformation& t,
                        shared_ptr<TC0> c0,
                        shared_ptr<TC1> c1,
                        shared_ptr<TC2> c2,
                        shared_ptr<TC3> c3,
                        shared_ptr<TC4> c4,
                        shared_ptr<TC5> c5,
                        shared_ptr<TC6> c6,
                        shared_ptr<TC7> c7,
                        shared_ptr<TC8> c8,
                        shared_ptr<TC9> c9)
      {
        return make_shared<typename result<TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9>::type>(gfs->backend().blocked(),gfs->orderingTag(),c0,c1,c2,c3,c4,c5,c6,c7,c8,c9);
      }

    };

#endif // HAVE_VARIADIC_TEMPLATES

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_LEXICOGRAPHICORDERING_HH
