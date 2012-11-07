// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_FILTEREDCOMPOSITENODE_HH
#define DUNE_PDELAB_COMMON_TYPETREE_FILTEREDCOMPOSITENODE_HH

#if !(HAVE_VARIADIC_TEMPLATES || DOXYGEN || HEADERCHECK)
#error The class FilteredCompositeNode requires compiler support for variadic templates, which your compiler lacks.
#endif

#if (HAVE_VARIADIC_TEMPLATES || DOXYGEN)
#include <dune/pdelab/common/typetree/nodetags.hh>
#include <dune/pdelab/common/typetree/filters.hh>
#include <dune/common/tuples.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Nodes
       *  \ingroup TypeTree
       *  \{
       */

#ifndef DOXYGEN
      namespace {

        // ********************************************************************************
        // Utility structs for filter construction and application
        // ********************************************************************************

        // Gets the filter and wraps it in case of a SimpleFilter.
        template<typename Filter, typename Tag>
        struct get_filter;

        // Helper struct to extract the child template parameter pack from the ChildTypes tuple.
        template<typename Filter, typename Node, typename ChildTypes>
        struct apply_filter_wrapper;

        template<typename Filter, typename Node, typename... Children>
        struct apply_filter_wrapper<Filter,Node,tuple<Children...> >
          : public Filter::template apply<Node,Children...>
        {};

        // specialization for SimpleFilter
        template<typename Filter>
        struct get_filter<Filter,SimpleFilterTag>
        {
          struct type
          {
            template<typename Node, typename ChildTypes>
            struct apply
              : public apply_filter_wrapper<filter<Filter>,Node,ChildTypes>
            {};
          };
        };

        // specialization for AdvancedFilter
        template<typename Filter>
        struct get_filter<Filter,AdvancedFilterTag>
        {
          struct type
          {
            template<typename Node, typename ChildTypes>
            struct apply
              : public apply_filter_wrapper<Filter,Node,ChildTypes>
            {};
          };
        };

      } // anonymous namespace
#endif // DOXYGEN


      //! Base class for composite nodes representing a filtered view on an underlying composite node.
      template<typename Node, typename Filter>
      class FilteredCompositeNode
      {

        typedef typename get_filter<Filter,typename Filter::FilterTag>::type filter;
        typedef typename filter::template apply<Node,typename Node::ChildTypes>::type filter_result;
        typedef typename filter_result::template apply<Node> mapped_children;

        static const bool nodeIsConst = IsConst<typename remove_reference<Node>::type>::value;

        template<std::size_t k>
        struct lazy_enable
        {
          static const bool value = !nodeIsConst;
        };

      public:

        //! The type tag that describes a VariadicCompositeNode.
        typedef VariadicCompositeNodeTag NodeTag;

        //! The type used for storing the children.
        typedef typename mapped_children::NodeStorage NodeStorage;

        //! A tuple storing the types of all children.
        typedef typename mapped_children::ChildTypes ChildTypes;

        //! Mark this class as non leaf in the \ref TypeTree.
        static const bool isLeaf = false;

        //! Mark this class as a non power in the \ref TypeTree.
        static const bool isPower = false;

        //! Mark this class as a composite in the \ref TypeTree.
        static const bool isComposite = true;

        //! The number of children.
        static const std::size_t CHILDREN = filter_result::size;

        //! Access to the type and storage type of the i-th child.
        template<std::size_t k>
        struct Child {

#ifndef DOXYGEN

          typedef typename tuple_element<k,typename mapped_children::Children>::type OriginalChild;

          static const std::size_t mapped_index = tuple_element<k,typename filter_result::IndexMap>::type::original_index;

#endif // DOXYGEN

          //! The type of the child.
          typedef typename OriginalChild::Type Type;

          //! The type of the child.
          typedef typename OriginalChild::type type;

          //! The storage type of the child.
          typedef typename OriginalChild::Storage Storage;

          //! The const storage type of the child.
          typedef typename OriginalChild::ConstStorage ConstStorage;
        };

        //! @name Child Access
        //! @{

        //! Returns the i-th child.
        /**
         * \returns a reference to the i-th child.
         */
        template<std::size_t k>
        typename enable_if<lazy_enable<k>::value,typename Child<k>::Type&>::type
        child()
        {
          return _node->template child<Child<k>::mapped_index>();
        }

        //! Returns the i-th child (const version).
        /**
         * \returns a const reference to the i-th child.
         */
        template<std::size_t k>
        const typename Child<k>::Type& child() const
        {
          return _node->template child<Child<k>::mapped_index>();
        }

        //! Returns the storage of the i-th child.
        /**
         * \returns a copy of the object storing the i-th child.
         */
        template<std::size_t k>
        typename enable_if<lazy_enable<k>::value,typename Child<k>::Storage>::type
        childStorage()
        {
          return _node->template childStorage<Child<k>::mapped_index>();
        }

        //! Returns the storage of the i-th child (const version).
        /**
         * This method is only important if the child is stored as
         * some kind of pointer, as this allows the pointee type to
         * become const.
         * \returns a copy of the object storing the i-th child.
         */
        template<std::size_t k>
        typename Child<k>::ConstStorage childStorage() const
        {
          return _node->template childStorage<Child<k>::mapped_index>();
        }

        //! Sets the i-th child to the passed-in value.
        template<std::size_t k>
        void setChild(typename Child<k>::type& child, typename enable_if<lazy_enable<k>::value,void*>::type = 0)
        {
          _node->template childStorage<Child<k>::mapped_index>() = stackobject_to_shared_ptr(child);
        }

        //! Sets the storage of the i-th child to the passed-in value.
        template<std::size_t k>
        void setChild(typename Child<k>::storage_type child, typename enable_if<lazy_enable<k>::value,void*>::type = 0)
        {
          _node->template childStorage<Child<k>::mapped_index>() = child;
        }

        //! @}

        //! @name Access to unfiltered node
        //! @{

      protected:

        //! Returns the unfiltered node.
        /**
         * \returns A reference to the original, unfiltered node.
         */
        template<bool enabled = !nodeIsConst>
        typename enable_if<enabled,Node&>::type
        unfiltered()
        {
          return *_node;
        }

        //! Returns the unfiltered node (const version).
        /**
         * \returns A const reference to the original, unfiltered node.
         */
        const Node& unfiltered() const
        {
          return *_node;
        }

        //! Returns the storage object of the unfiltered node.
        /**
         * \returns A shared_ptr to the original, unfiltered node.
         */
        template<bool enabled = !nodeIsConst>
        typename enable_if<enabled,shared_ptr<Node> >::type
        unfilteredStorage()
        {
          return _node;
        }

        //! Returns the storage object of the unfiltered node (const version).
        /**
         * \returns A shared_ptr to the original, unfiltered node.
         */
        shared_ptr<const Node> unfilteredStorage() const
        {
          return _node;
        }

        //! @}

      public:

        //! @name Constructors
        //! @{

        //! Initialize the VariadicCompositeNode with copies of the passed in Storage objects.
        FilteredCompositeNode(shared_ptr<Node> node)
          : _node(node)
        {}

        //! Initialize the VariadicCompositeNode with a copy of the passed-in storage type.
        FilteredCompositeNode(Node& node)
          : _node(stackobject_to_shared_ptr(node))
        {}

        //! @}

      private:
        shared_ptr<Node> _node;
      };

      //! \} group Nodes

    } // namespace TypeTree

  } // namespace PDELab
} //namespace Dune

#endif // (HAVE_VARIADIC_TEMPLATES || DOXYGEN)

#endif // DUNE_PDELAB_COMMON_TYPETREE_FILTEREDCOMPOSITENODE_HH
