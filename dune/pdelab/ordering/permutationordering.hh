// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_PERMUTATIONORDERING_HH
#define DUNE_PDELAB_ORDERING_PERMUTATIONORDERING_HH

#include <cstddef>
#include <ostream>
#include <string>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/powernode.hh>
#include <dune/typetree/traversal.hh>
#include <dune/typetree/visitor.hh>

#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/ordering/lexicographicordering.hh>
#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/orderingbase.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    namespace permutation_ordering {

      //! Interface for merging index spaces
      template<typename DI, typename CI, typename Node>
      class Base
        : public lexicographic_ordering::Base<DI,CI,Node>,
          public VirtualOrderingBase<DI,CI>
      {

      public:

        typedef Base<DI,CI,Node> This;
        typedef lexicographic_ordering::Base<DI,CI,Node> BaseT;
        typedef typename BaseT::Traits Traits;

        typedef PermutationOrderingTag OrderingTag;

        static const bool consume_tree_index = true;

        //! Construct ordering object
        /**
         * In general, an ordering object is not properly setup after
         * construction.  This must be done by a separate call to update()
         * after all the children have been properly set up.
         */
        Base(Node& node, bool container_blocked, const OrderingTag& ordering_tag)
          : BaseT(node,container_blocked),
            _perm(ordering_tag.permutation())
        {
          // Make sure to use 'this' as a delegate in OrderingBase, so the virtual function call
          // to 'map_index_dynamic' is found!
          this->setDelegate(this);
        }

        //! Copy constructor
        Base(const This& other)
          : BaseT(other),
            _perm(other._perm)
        {
          // Make sure to use 'this' as a delegate in OrderingBase, so the virtual function call
          // to 'map_index_dynamic' is found!
          this->setDelegate(this);
        }

        Base(This&& other)
          : BaseT(std::move(other)),
            _perm(std::move(other._perm))
        {
          // Make sure to use 'this' as a delegate in OrderingBase, so the virtual function call
          // to 'map_index_dynamic' is found!
          this->setDelegate(this);
        }

        template<typename ItIn, typename ItOut>
        void map_lfs_indices(const ItIn begin, const ItIn end, ItOut out) const
        {
          BaseT::map_lfs_indices(begin,end,out);

          // Permute indices
          for (ItIn in = begin; in != end; ++in, ++out)
            out->back() = _perm[out->back()];
        }

        template<typename CIOutIterator, typename DIOutIterator = DummyDOFIndexIterator>
        typename Traits::SizeType
        extract_entity_indices(const typename Traits::DOFIndex::EntityIndex& ei,
                               typename Traits::SizeType child_index,
                               CIOutIterator ci_out, const CIOutIterator ci_end) const
        {
          BaseT::extract_entity_indices(ei,child_index,ci_out,ci_end);

          // Permute indices
          for (; ci_out != ci_end; ++ci_out)
          {
            ci_out->back() = _perm[ci_out->back()];
          }

          // The return value is not used for non-leaf orderings.
          return 0;
        }

        virtual void map_index_dynamic(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
        {
          this->_mapIndex(di,ci);

          // Permute index
          ci.back() = _perm[ci.back()];
        }

      private:
        const std::vector<std::size_t>& _perm;

      };
    }

    //! Interface for merging index spaces
    template<typename DI, typename CI, typename Child, std::size_t k>
    class PowerPermutationOrdering
      : public TypeTree::PowerNode<Child, k>
      , public permutation_ordering::Base<DI,
                                            CI,
                                            PowerLexicographicOrdering<DI,CI,Child,k>
                                            >
    {
      typedef TypeTree::PowerNode<Child, k> Node;

      typedef permutation_ordering::Base<DI,
                                           CI,
                                           PowerLexicographicOrdering<DI,CI,Child,k>
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
      PowerPermutationOrdering(bool container_blocked, const PermutationOrderingTag& ordering_tag, const typename Node::NodeStorage& children)
        : Node(children)
        , Base(*this,container_blocked,ordering_tag)
      { }

      void update()
      {
        for (std::size_t i = 0; i < Node::CHILDREN; ++i)
          {
            this->child(i).update();
          }
        Base::update();
      }

      std::string name() const { return "PowerPermutationOrdering"; }
    };


    template<typename GFS, typename Transformation>
    struct power_gfs_to_ordering_descriptor<GFS,Transformation,PermutationOrderingTag>
    {

      static const bool recursive = true;

      template<typename TC>
      struct result
      {

        typedef PowerPermutationOrdering<
          typename Transformation::DOFIndex,
          typename Transformation::ContainerIndex,
          TC,
          GFS::CHILDREN
          > type;

        typedef shared_ptr<type> storage_type;

      };

      template<typename TC>
      static typename result<TC>::type transform(const GFS& gfs, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
      {
        return typename result<TC>::type(gfs.backend().blocked(gfs),gfs->orderingTag(),children);
      }

      template<typename TC>
      static typename result<TC>::storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
      {
        return make_shared<typename result<TC>::type>(gfs->backend().blocked(*gfs),gfs->orderingTag(),children);
      }

    };

    // the generic registration for PowerGridFunctionSpace happens in transformations.hh


    //! Interface for merging index spaces
    template<typename DI, typename CI, typename... Children>
    class CompositePermutationOrdering :
      public TypeTree::CompositeNode<Children...>,
      public permutation_ordering::Base<DI,
                                          CI,
                                          CompositePermutationOrdering<
                                            DI,
                                            CI,
                                            Children...
                                            >
                                          >
    {
      typedef TypeTree::CompositeNode<Children...> Node;

      typedef permutation_ordering::Base<
        DI,
        CI,
        CompositePermutationOrdering<
          DI,
          CI,
          Children...
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
      CompositePermutationOrdering(bool backend_blocked, const PermutationOrderingTag& ordering_tag,
          Dune::shared_ptr<Children>... children)
        : Node(children...)
        , Base(*this,backend_blocked,ordering_tag)
      { }

      std::string name() const { return "CompositePermutationOrdering"; }

      void update()
      {
        TypeTree::applyToTree(*this,ordering::update_direct_children());
        Base::update();
      }
    };

    template<typename GFS, typename Transformation>
    struct composite_gfs_to_ordering_descriptor<GFS,Transformation,PermutationOrderingTag>
    {

      static const bool recursive = true;

      template<typename... TC>
      struct result
      {

        typedef CompositePermutationOrdering<
          typename Transformation::DOFIndex,
          typename Transformation::ContainerIndex,
          TC...
          > type;

        typedef shared_ptr<type> storage_type;

      };

      template<typename... TC>
      static typename result<TC...>::type transform(const GFS& gfs, const Transformation& t, shared_ptr<TC>... children)
      {
        return typename result<TC...>::type(gfs.backend().blocked(gfs),gfs.orderingTag(),children...);
      }

      template<typename... TC>
      static typename result<TC...>::storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t, shared_ptr<TC>... children)
      {
        return make_shared<typename result<TC...>::type>(gfs->backend().blocked(*gfs),gfs.orderingTag(),children...);
      }

    };

    // the generic registration for PowerGridFunctionSpace happens in transformations.hh

   //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_PERMUTATIONORDERING_HH
