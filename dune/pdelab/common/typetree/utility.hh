// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_UTILITY_HH
#define DUNE_PDELAB_COMMON_TYPETREE_UTILITY_HH

#include <dune/common/shared_ptr.hh>
#include <dune/common/tuples.hh>
#include <dune/pdelab/common/typetree/nodetags.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup TypeTree
       *  \ingroup PDELab
       *  \{
       */

#ifndef DOXYGEN

      template<typename T>
      shared_ptr<T> convert_arg(const T& t)
      {
        return make_shared<T>(t);
      }

      template<typename T>
      shared_ptr<T> convert_arg(T& t)
      {
        return stackobject_to_shared_ptr(t);
      }

      template<typename BaseType, typename T>
      T& assertGridViewType(T& t)
      {
        dune_static_assert((is_same<typename BaseType::Traits::GridViewType,
                            typename T::Traits::GridViewType>::value),
                           "GridViewType must be equal in all components of composite type");
        return t;
      }

      // this partial specialization is required for the non-variadic case
      template<typename BaseType>
      TypeTree::EmptyNode assertGridViewType(TypeTree::EmptyNode e)
      {
        return e;
      }

#if HAVE_RVALUE_REFERENCES

      // only bind to real rvalues
      template<typename T>
      typename enable_if<!std::is_lvalue_reference<T>::value,shared_ptr<T> >::type convert_arg(T&& t)
      {
      return make_shared<T>(std::forward<T>(t));
    }

#endif

#endif // DOXYGEN


      namespace {

        //! Pointer to an empty node that is used for all empty slots
        /**
         * TODO: move into a library!
         */
        shared_ptr<EmptyNode> emptyNodePtr(make_shared<EmptyNode>());

      }


      //! Struct for obtaining some basic structural information about a TypeTree.
      /**
       * This struct extracts basic information about the passed TypeTree and
       * presents them in a static way suitable for use as compile-time constants.
       *
       * \tparam Tree  The TypeTree to examine.
       * \tparam Tag   Internal parameter, leave at default value.
       */
      template<typename Tree, typename Tag = StartTag>
      struct TreeInfo
      {

      private:
        // Start the tree traversal
        typedef TreeInfo<Tree,typename Tree::NodeTag> NodeInfo;

      public:

        //! The depth of the TypeTree.
        static const std::size_t depth = NodeInfo::depth;

        //! The total number of nodes in the TypeTree.
        static const std::size_t nodeCount = NodeInfo::nodeCount;

        //! The number of leaf nodes in the TypeTree.
        static const std::size_t leafCount = NodeInfo::leafCount;

      };


#ifndef DOXYGEN

      // ********************************************************************************
      // TreeInfo specializations for the different node types
      // ********************************************************************************


      // leaf node
      template<typename Node>
      struct TreeInfo<Node,LeafNodeTag>
      {

        static const std::size_t depth = 1;

        static const std::size_t nodeCount = 1;

        static const std::size_t leafCount = 1;

      };


      // power node - exploit the fact that all children are identical
      template<typename Node>
      struct TreeInfo<Node,PowerNodeTag>
      {

        typedef TreeInfo<typename Node::ChildType,typename Node::ChildType::NodeTag> ChildInfo;

        static const std::size_t depth = 1 + ChildInfo::depth;

        static const std::size_t nodeCount = 1 + Node::CHILDREN * ChildInfo::nodeCount;

        static const std::size_t leafCount = Node::CHILDREN * ChildInfo::leafCount;

      };


      namespace {

        // TMP for iterating over the children of a composite node
        // identical for both composite node implementations
        template<typename Node, std::size_t k, std::size_t n>
        struct generic_compositenode_children_info
        {

          typedef generic_compositenode_children_info<Node,k+1,n> NextChild;

          // extract child info
          typedef typename Node::template Child<k>::Type Child;
          typedef typename Child::NodeTag ChildTag;
          typedef TreeInfo<Child,ChildTag> ChildInfo;

          // combine information of current child with info about following children
          static const std::size_t maxDepth = ChildInfo::depth > NextChild::maxDepth ? ChildInfo::depth : NextChild::maxDepth;

          static const std::size_t nodeCount = ChildInfo::nodeCount + NextChild::nodeCount;

          static const std::size_t leafCount = ChildInfo::leafCount + NextChild::leafCount;

        };

        // End of recursion
        template<typename Node, std::size_t n>
        struct generic_compositenode_children_info<Node,n,n>
        {
          static const std::size_t maxDepth = 0;

          static const std::size_t nodeCount = 0;

          static const std::size_t leafCount = 0;
        };

      } // anonymous namespace


      // Struct for building information about composite node
      template<typename Node>
      struct GenericCompositeNodeInfo
      {

        typedef generic_compositenode_children_info<Node,0,Node::CHILDREN> Children;

        static const std::size_t depth = 1 + Children::maxDepth;

        static const std::size_t nodeCount = 1 + Children::nodeCount;

        static const std::size_t leafCount = Children::leafCount;

      };


      // CompositeNode: delegate to GenericCompositeNodeInfo
      template<typename Node>
      struct TreeInfo<Node,CompositeNodeTag>
        : public GenericCompositeNodeInfo<Node>
      {};


      // VariadicCompositeNode: delegate to GenericCompositeNodeInfo
      template<typename Node>
      struct TreeInfo<Node,VariadicCompositeNodeTag>
        : public GenericCompositeNodeInfo<Node>
      {};

#endif // DOXYGEN


#if HAVE_VARIADIC_TEMPLATES

      //! Simple holder class for a template argument pack of indices.
      /**
       * The main use of index_pack is to unpack variadically templated
       * data structures like this:
       *
       * \code
       * template<typename T, typename F, std::size_t... i>
       * void apply_to_tuple(const T& t, F f, index_pack<i...> indices)
       * {
       *   f(get<i>(t))...;
       * }
       *
       * tuple<int,double,...,char> t;
       * apply_to_tuple(t,foo,tuple_indices(t));
       * \endcode
       *
       * \sa tuple_indices()
       */
      template<std::size_t... i>
      struct index_pack {};

      //! TMP to build an index_pack containing the sequence 0,...,n-1.
      template<std::size_t n, std::size_t... i>
      struct build_index_pack
        : public build_index_pack<n-1,n-1,i...>
      {

#ifdef DOXYGEN
        //! Result.
        typedef index_pack<0,1,...,n-1> type;
#endif // DOXYGEN

      };

#ifndef DOXYGEN

      // end of recursion
      template<std::size_t... i>
      struct build_index_pack<0,i...>
      {
        typedef index_pack<i...> type;
      };

#endif // DOXYGEN

      //! TMP to build an index_pack for all elements in the tuple.
      template<typename tuple>
      struct build_tuple_index_pack
        : public build_index_pack<tuple_size<tuple>::value>
      {};

      //! Generate an index_pack for the tuple t.
      template<typename tuple>
      typename build_tuple_index_pack<tuple>::type tuple_indices(const tuple& t)
      {
        return typename build_tuple_index_pack<tuple>::type();
      }

#endif // HAVE_VARIADIC_TEMPLATES

      //! \} group TypeTree

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_UTILITY_HH
