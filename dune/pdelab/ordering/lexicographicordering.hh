// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_LEXICOGRAPHICORDERING_HH
#define DUNE_PDELAB_ORDERING_LEXICOGRAPHICORDERING_HH

#include <cstddef>
#include <ostream>
#include <string>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/powernode.hh>
#include <dune/typetree/traversal.hh>
#include <dune/typetree/visitor.hh>

#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/orderingbase.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    namespace lexicographic_ordering {

      template<typename DI, typename CI, typename Node>
      class Base
        : public OrderingBase<DI,CI>
      {

        typedef OrderingBase<DI,CI> BaseT;

      public:

        typedef typename OrderingBase<DI,CI>::Traits Traits;

        typedef LexicographicOrderingTag OrderingTag;

        static const bool consume_tree_index = true;

        //! Construct ordering object
        /**
         * In general, an ordering object is not properly setup after
         * construction.  This must be done by a seperate call to update()
         * after all the children have been properly set up.
         */
        Base(Node& node, bool container_blocked, typename BaseT::GFSData* gfs_data)
          : BaseT(node,container_blocked,gfs_data,nullptr)
        {
        }

        template<typename ItIn, typename ItOut>
        void map_lfs_indices(const ItIn begin, const ItIn end, ItOut out) const
        {
          if (this->_container_blocked)
            {
              for (ItIn in = begin; in != end; ++in, ++out)
                out->push_back(in->treeIndex().back());
            }
          else
            {
              for (ItIn in = begin; in != end; ++in, ++out)
                out->back() += (this->blockOffset(in->treeIndex().back()));
            }
        }

        template<typename CIOutIterator, typename DIOutIterator = DummyDOFIndexIterator>
        typename Traits::SizeType
        extract_entity_indices(const typename Traits::DOFIndex::EntityIndex& ei,
                               typename Traits::SizeType child_index,
                               CIOutIterator ci_out, const CIOutIterator ci_end) const
        {
          if (this->_container_blocked)
            {
              for (; ci_out != ci_end; ++ci_out)
                {
                  ci_out->push_back(child_index);
                }
            }
          else
            {
              for (; ci_out != ci_end; ++ci_out)
                {
                  ci_out->back() += (this->blockOffset(child_index));
                }
            }

          // The return value is not used for non-leaf orderings.
          return 0;
        }

      };
    }



    template<typename DI, typename CI, typename Child, std::size_t k>
    class PowerLexicographicOrdering
      : public TypeTree::PowerNode<Child, k>
      , public lexicographic_ordering::Base<DI,
                                            CI,
                                            PowerLexicographicOrdering<DI,CI,Child,k>
                                            >
    {
      typedef TypeTree::PowerNode<Child, k> Node;

      typedef lexicographic_ordering::Base<DI,
                                           CI,
                                           PowerLexicographicOrdering<DI,CI,Child,k>
                                           > Base;

    public:

      using Traits = typename Base::Traits;

      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update() after
       * all the children have been properly set up.
       *
       * \note This constructor must be present for ordering objects not at
       *       the leaf of the tree.
       */
      PowerLexicographicOrdering(bool container_blocked, const typename Node::NodeStorage& children, typename Base::GFSData* gfs_data)
        : Node(children)
        , Base(*this,container_blocked,gfs_data)
      { }

      void update()
      {
        for (std::size_t i = 0; i < k; ++i)
          {
            this->child(i).update();
          }
        Base::update();
      }

      std::string name() const { return "PowerLexicographicOrdering"; }

      using Base::size;

      /**
       * @brief Gives the size for a given prefix
       * @param prefix  MultiIndex with a partial path to a container
       * @return Traits::SizeType  The size required for such a path.
       */
      typename Traits::SizeType size(typename Traits::SizePrefix prefix) const {
        if (prefix.size() == Traits::SizePrefix::max_depth)
          return 0; // all indices in prefix were consumed, no more sizes to provide

        if (prefix.size() == 0)
          return this->blockCount();

        if (this->containerBlocked()) {
          auto child = prefix.back();
          assert(this->degree() > child);
          prefix.pop_back();
          return this->child(child).size(prefix);
        } else {
          auto it = std::upper_bound(this->_child_block_offsets.begin(), this->_child_block_offsets.end(), prefix.back());
          std::size_t child = *std::prev(it);
          return this->child(child).size(prefix);
        }
      }

    };


    template<typename GFS, typename Transformation>
    struct power_gfs_to_lexicographic_ordering_descriptor
    {

      static const bool recursive = true;

      template<typename TC>
      struct result
      {

        typedef PowerLexicographicOrdering<
          typename Transformation::DOFIndex,
          typename Transformation::ContainerIndex,
          TC,
          TypeTree::StaticDegree<GFS>::value
          > type;

        typedef std::shared_ptr<type> storage_type;

      };

      template<typename TC>
      static typename result<TC>::type transform(const GFS& gfs, const Transformation& t, const std::array<std::shared_ptr<TC>,TypeTree::StaticDegree<GFS>::value>& children)
      {
        return typename result<TC>::type(gfs.backend().blocked(gfs),children,const_cast<GFS*>(&gfs));
      }

      template<typename TC>
      static typename result<TC>::storage_type transform_storage(std::shared_ptr<const GFS> gfs, const Transformation& t, const std::array<std::shared_ptr<TC>,TypeTree::StaticDegree<GFS>::value>& children)
      {
        return std::make_shared<typename result<TC>::type>(gfs->backend().blocked(*gfs),children,const_cast<GFS*>(gfs.get()));
      }

    };

    template<typename GFS, typename Transformation>
    power_gfs_to_lexicographic_ordering_descriptor<GFS,Transformation>
    register_power_gfs_to_ordering_descriptor(GFS*,Transformation*,LexicographicOrderingTag*);

    // the generic registration for PowerGridFunctionSpace happens in transformations.hh


    //! Interface for merging index spaces
    template<typename DI, typename CI, typename... Children>
    class CompositeLexicographicOrdering :
      public TypeTree::CompositeNode<Children...>,
      public lexicographic_ordering::Base<DI,
                                          CI,
                                          CompositeLexicographicOrdering<
                                            DI,
                                            CI,
                                            Children...
                                            >
                                          >
    {
      typedef TypeTree::CompositeNode<Children...> Node;

      typedef lexicographic_ordering::Base<
        DI,
        CI,
        CompositeLexicographicOrdering<
          DI,
          CI,
          Children...
          >
        > Base;

    public:
      using Traits = typename Base::Traits;

      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update() after
       * all the children have been properly set up.
       *
       * \note This constructor must be present for ordering objects not at
       *       the leaf of the tree.
       */
      CompositeLexicographicOrdering(bool backend_blocked, typename Base::GFSData* gfs_data, std::shared_ptr<Children>... children)
        : Node(children...)
        , Base(*this,backend_blocked,gfs_data)
      { }

      std::string name() const { return "CompositeLexicographicOrdering"; }

      void update()
      {
        TypeTree::applyToTree(*this,ordering::update_direct_children());
        Base::update();
      }

      using Base::size;

      /**
       * @brief Gives the size for a given prefix
       * @param prefix  MultiIndex with a partial path to a container
       * @return Traits::SizeType  The size required for such a path.
       */
      typename Traits::SizeType size(typename Traits::SizePrefix prefix) const {
        if (prefix.size() == Traits::SizePrefix::max_depth)
          return 0; // all indices in prefix were consumed, no more sizes to provide

        if (prefix.size() == 0)
          return this->blockCount();

        auto indices = std::make_index_sequence<Node::degree()>{};
        typename Traits::SizeType _size;
        std::size_t _child;

        if (this->containerBlocked()) {
          _child = prefix.back();
          assert(this->degree() > _child);
          prefix.pop_back();
        } else {
          auto it = std::upper_bound(this->_child_block_offsets.begin(), this->_child_block_offsets.end(), prefix.back());
          _child = *std::prev(it);
        }

        Hybrid::forEach(indices, [&](auto i){
          if (i == _child)
            _size = this->template child<i>().size(prefix);
        });
        return _size;
      }
    };

    template<typename GFS, typename Transformation>
    struct composite_gfs_to_lexicographic_ordering_descriptor
    {

      static const bool recursive = true;

      template<typename... TC>
      struct result
      {

        typedef CompositeLexicographicOrdering<
          typename Transformation::DOFIndex,
          typename Transformation::ContainerIndex,
          TC...
          > type;

        typedef std::shared_ptr<type> storage_type;

      };

      template<typename... TC>
      static typename result<TC...>::type transform(const GFS& gfs, const Transformation& t, std::shared_ptr<TC>... children)
      {
        return typename result<TC...>::type(gfs.backend().blocked(gfs),const_cast<GFS*>(&gfs),children...);
      }

      template<typename... TC>
      static typename result<TC...>::storage_type transform_storage(std::shared_ptr<const GFS> gfs, const Transformation& t, std::shared_ptr<TC>... children)
      {
        return std::make_shared<typename result<TC...>::type>(gfs->backend().blocked(*gfs),const_cast<GFS*>(gfs.get()),children...);
      }

    };


    template<typename GFS, typename Transformation>
    composite_gfs_to_lexicographic_ordering_descriptor<GFS,Transformation>
    register_composite_gfs_to_ordering_descriptor(GFS*,Transformation*,LexicographicOrderingTag*);

   //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_LEXICOGRAPHICORDERING_HH
