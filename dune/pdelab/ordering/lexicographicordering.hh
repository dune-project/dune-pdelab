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

#include <dune/pdelab/common/typetree/compositenodemacros.hh>
#include <dune/pdelab/common/typetree/powernode.hh>
#include <dune/pdelab/common/typetree/traversal.hh>
#include <dune/pdelab/common/typetree/visitor.hh>

#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/orderingbase.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    namespace lexicographic_ordering {

      //! Interface for merging index spaces
      template<typename DI, typename CI, typename Node>
      class Base
        : public OrderingBase<DI,CI>
      {

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
        Base(Node& node, bool container_blocked)
        : OrderingBase<DI,CI>(node,container_blocked,nullptr)
        {
        }

        template<typename ItIn, typename ItOut>
        void map_indices(const ItIn begin, const ItIn end, ItOut out) const
        {
          if (this->_container_blocked)
            {
              for (ItIn in = begin; in != end; ++in, ++out)
                out->push_back(in->treeIndex().back());
            }
          else
            {
              for (ItIn in = begin; in != end; ++in, ++out)
                out->back() += (this->offset(in->treeIndex().back()));
            }
        }

        template<typename CIOutIterator, typename DIOutIterator = DummyDOFIndexIterator>
        typename Traits::SizeType
        containerIndices(const typename Traits::DOFIndex::EntityIndex& ei,
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
                  ci_out->back() += (this->offset(child_index));
                }
            }

          // The return value is not used for non-leaf orderings.
          return 0;
        }

      };
    }

    //! Interface for merging index spaces
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

      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update() after
       * all the children have been properly set up.
       *
       * \note This constructor must be present for ordering objects not at
       *       the leaf of the tree.
       */
      PowerLexicographicOrdering(bool container_blocked, const typename Node::NodeStorage& children)
        : Node(children)
        , Base(*this,container_blocked)
      { }

      void update()
      {
        for (std::size_t i = 0; i < Node::CHILDREN; ++i)
          {
            this->child(i).update();
          }
        Base::update();
      }

      std::string name() const { return "PowerLexicographicOrdering"; }
    };


    template<typename GFS, typename Transformation, typename OrderingTag>
    struct power_gfs_to_ordering_descriptor;

    template<typename GFS, typename Transformation>
    struct power_gfs_to_ordering_descriptor<GFS,Transformation,LexicographicOrderingTag>
    {

      static const bool recursive = true;

      template<typename TC>
      struct result
      {

        typedef PowerLexicographicOrdering<
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
        return typename result<TC>::type(gfs.backend().blocked(),children);
      }

      template<typename TC>
      static typename result<TC>::storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
      {
        return make_shared<typename result<TC>::type>(gfs->backend().blocked(),children);
      }

    };


    template<typename GridFunctionSpace, typename Params>
    power_gfs_to_ordering_descriptor<
      GridFunctionSpace,
      gfs_to_ordering<Params>,
      typename GridFunctionSpace::OrderingTag
      >
    lookupNodeTransformation(GridFunctionSpace* gfs, gfs_to_ordering<Params>* t, PowerGridFunctionSpaceTag tag);


    struct update_direct_children
      : public TypeTree::DirectChildrenVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename GFS, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(const GFS& gfs, Child& child, TreePath tp, ChildIndex childIndex) const
      {
        child.update();
      }

    };


    //! Interface for merging index spaces
    template<typename DI, typename CI, DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN>
    class CompositeLexicographicOrdering :
      public DUNE_TYPETREE_COMPOSITENODE_BASETYPE,
      public lexicographic_ordering::Base<DI,
                                          CI,
                                          CompositeLexicographicOrdering<
                                            DI,
                                            CI,
                                            DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES
                                            >
                                          >
    {
      typedef DUNE_TYPETREE_COMPOSITENODE_BASETYPE Node;

      typedef lexicographic_ordering::Base<
        DI,
        CI,
        CompositeLexicographicOrdering<
          DI,
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
      CompositeLexicographicOrdering(bool backend_blocked, DUNE_TYPETREE_COMPOSITENODE_STORAGE_CONSTRUCTOR_SIGNATURE)
        : Node(DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES)
        , Base(*this,backend_blocked)
      { }

      std::string name() const { return "CompositeLexicographicOrdering"; }

      void update()
      {
        TypeTree::applyToTree(*this,update_direct_children());
        Base::update();
      }
    };

#if HAVE_VARIADIC_TEMPLATES

    template<typename GFS, typename Transformation, typename OrderingTag>
    struct composite_gfs_to_ordering_descriptor;

    template<typename GFS, typename Transformation>
    struct composite_gfs_to_ordering_descriptor<GFS,Transformation,LexicographicOrderingTag>
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

        typedef shared_ptr<type> storage_type;

      };

      template<typename... TC>
      static typename result<TC...>::type transform(const GFS& gfs, const Transformation& t, shared_ptr<TC>... children)
      {
        return typename result<TC...>::type(gfs.backend().blocked(),children...);
      }

      template<typename... TC>
      static typename result<TC...>::storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t, shared_ptr<TC>... children)
      {
        return make_shared<typename result<TC...>::type>(gfs->backend().blocked(),children...);
      }

    };


    template<typename GridFunctionSpace, typename Params>
    composite_gfs_to_ordering_descriptor<
      GridFunctionSpace,
      gfs_to_ordering<Params>,
      typename GridFunctionSpace::OrderingTag
      >
    lookupNodeTransformation(GridFunctionSpace* gfs, gfs_to_ordering<Params>* t, CompositeGridFunctionSpaceTag tag);


#else // HAVE_VARIADIC_TEMPLATES

    //! Node transformation descriptor for CompositeGridFunctionSpace -> LexicographicOrdering (without variadic templates).
    template<typename GFSNode, typename Transformation>
    struct CompositeGFSToLexicographicOrderingTransformation
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
        typedef CompositeLexicographicOrdering<typename Transformation::GridFunctionSpace::Traits::SizeType,
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
      transform(const GFSNode& s,
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
        return typename result<TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9>::type(t.asGridFunctionSpace(s),c0,c1,c2,c3,c4,c5,c6,c7,c8,c9);
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
      transform_storage(shared_ptr<const GFSNode> s,
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
        return make_shared<typename result<TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9>::type>(t.asGridFunctionSpace(s),c0,c1,c2,c3,c4,c5,c6,c7,c8,c9);
      }

    };

    // Register transformation
    template<typename GFSNode, typename GFS>
    CompositeGFSToLexicographicOrderingTransformation<GFSNode,gfs_to_ordering<GFS,LexicographicOrderingTag> >
    lookupNodeTransformation(GFSNode*, gfs_to_ordering<GFS,LexicographicOrderingTag>*, CompositeGridFunctionSpaceBaseTag);

#endif // HAVE_VARIADIC_TEMPLATES

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_LEXICOGRAPHICORDERING_HH
