// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_FILTEREDCOMPOSITENODE_HH
#define DUNE_PDELAB_COMMON_TYPETREE_FILTEREDCOMPOSITENODE_HH

#if !(HAVE_VARIADIC_TEMPLATES || DOXYGEN || HEADERCHECK)
#error The class FilteredCompositeNode requires compiler support for variadic templates, which your compiler lacks.
#endif

#include <dune/pdelab/common/typetree/nodetags.hh>
#include <dune/common/tuples.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Nodes
       *  \ingroup TypeTree
       *  \{
       */

      namespace {

        template<typename T, std::size_t new_k, std::size_t old_k>
        struct filter_entry
        {
          typedef T type;
          typedef shared_ptr<T> storage_type;
          typedef shared_ptr<const T> const_storage_type;
          static const std::size_t filtered_index = new_k;
          static const std::size_t original_index = old_k;
        };

        template<typename Filter, std::size_t new_k, std::size_t old_k, typename... tail>
        struct filter_helper
        {
          template<typename... head>
          struct apply
          {
            typedef tuple<typename head::type...> types;
            typedef tuple<typename head::storage_type... > storage_types;
            typedef tuple<head...> filter_entry_types;

            static const std::size_t size = sizeof...(head);
          };
        };

        template<typename Filter, std::size_t new_k, std::size_t old_k, typename child, typename... tail>
        struct filter_helper<Filter,new_k,old_k,child,tail...>
        {

          template<typename... head>
          struct apply
          {
            typedef typename SelectType<Filter::template apply<child,new_k,old_k>::value,
                                        typename filter_helper<Filter,new_k+1,old_k+1,tail...>::template apply<head...,filter_entry<child,new_k,old_k> >,
                                        typename filter_helper<Filter,new_k,old_k+1,tail...>::template apply<head...>
                                        >::Type result;

            typedef typename result::types types;
            typedef typename result::storage_types storage_types;
            typedef typename result::filter_entry_types filter_entry_types;

            static const std::size_t size = result::size;
          };

        };

        template<typename Filter, typename ChildContainer>
        struct filter;

        template<typename Filter, typename... Children>
        struct filter<Filter,tuple<Children...> >
        {
          typedef typename filter_helper<Filter,0,0,Children...>::template apply<> result;

          typedef typename result::types types;
          typedef typename result::storage_types storage_types;
          typedef typename result::filter_entry_types filter_entry_types;

          static const std::size_t size = result::size;
        };

      }

      //! Base class for composite nodes representing a filtered view on an underlying composite node.
      template<typename Node, typename Filter>
      class FilteredCompositeNode
      {

        typedef filter<Filter,typename Node::ChildTypes> filter_result;
        typedef typename filter_result::filter_entry_types index_map;

      public:

        //! The type tag that describes a VariadicCompositeNode.
        typedef VariadicCompositeNodeTag NodeTag;

        //! The type used for storing the children.
        typedef typename filter_result::storage_types NodeStorage;

        //! A tuple storing the types of all children.
        typedef typename filter_result::types ChildTypes;

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

          typedef typename tuple_element<k,index_map>::type map_entry;

          static const std::size_t mapped_index = map_entry::original_index;

#endif // DOXYGEN

          //! The type of the child.
          typedef typename map_entry::type Type;

          //! The storage type of the child.
          typedef typename map_entry::storage_type Storage;

          //! The const storage type of the child.
          typedef typename map_entry::const_storage_type ConstStorage;
        };

        //! @name Child Access
        //! @{

        //! Returns the i-th child.
        /**
         * \returns a reference to the i-th child.
         */
        template<std::size_t k>
        typename Child<k>::Type& child()
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
        typename Child<k>::Storage childStorage()
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
        void setChild(typename Child<k>::type& child)
        {
          _node->template childStorage<Child<k>::mapped_index>() = stackobject_to_shared_ptr(child);
        }

        //! Sets the storage of the i-th child to the passed-in value.
        template<std::size_t k>
        void setChild(typename Child<k>::storage_type child)
        {
          _node->template childStorage<Child<k>::mapped_index>() = child;
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

#endif // DUNE_PDELAB_COMMON_TYPETREE_FILTEREDCOMPOSITENODE_HH
