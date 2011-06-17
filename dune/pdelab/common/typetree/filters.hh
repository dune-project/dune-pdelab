// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_FILTERS_HH
#define DUNE_PDELAB_COMMON_TYPETREE_FILTERS_HH

#if !(HAVE_VARIADIC_TEMPLATES || DOXYGEN || HEADERCHECK)
#error The class FilteredCompositeNode requires compiler support for variadic templates, which your compiler lacks.
#endif

#if (HAVE_VARIADIC_TEMPLATES || DOXYGEN)
#include <dune/common/tuples.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Nodes
       *  \ingroup TypeTree
       *  \{
       */

      //! A filter entry describing the mapping of one child in the filtered node.
      template<std::size_t new_k, std::size_t old_k>
      struct FilterEntry
      {

#ifndef DOXYGEN

        // The precise contents of this class is an implementation detail.

        static const std::size_t filtered_index = new_k;
        static const std::size_t original_index = old_k;

#endif // DOXYGEN

      };

      //! The result of a filter.
      template<typename... FilterEntries>
      struct FilterResult
      {

        static const std::size_t size = sizeof...(FilterEntries);

        typedef tuple<FilterEntries...> IndexMap;

        template<typename Node>
        struct apply
        {
          typedef tuple<typename Node::template Child<FilterEntries::original_index>...> Children;
          typedef tuple<typename Node::template Child<FilterEntries::original_index>::Type...> ChildTypes;
          typedef tuple<typename Node::template Child<FilterEntries::original_index>::Storage...> NodeStorage;
        };

      };

      //! Tag describing a simple filter that can only decide whether or not to include a single given child.
      struct SimpleFilterTag {};

      //! Tag describing an advanced filter that has full control over the construction of the list of FilterEntries.
      struct AdvancedFilterTag {};


      //! Base class for advanced filters.
      struct AdvancedFilter
      {

        //! Filter tag for deciding on filter application mechanism.
        typedef AdvancedFilterTag FilterTag;

#ifdef DOXYGEN

        //! Apply this filter to the given node and children
        template<typename Node, typename... Children>
        struct apply
        {
          //! The result of the filtering process.
          /**
           * This type must be a model of FilterResult.
           */
          typedef implementation-defined type;
        };

#endif // DOXYGEN

      };

      //! Default simple filter that accepts any node and leaves its child structure unchanged.
      /**
       * This default filter causes the filtered node to be exactly identical to the original node.
       * It is useful as a base class for documentation purposes and if you do not need to validate
       * the filter, as it saves having to implement the validate template struct.
       */
      struct SimpleFilter
      {

        //! Filter tag for deciding on filter application mechanism.
        typedef SimpleFilterTag FilterTag;


        //! Validates the combination of filter and node.
        template<typename Node>
        struct validate
        {
          //! True if the combination of filter and node is valid.
          static const bool value = true;
        };

        //! Applies the filter to the given child node.
        /**
         * This struct applies the filter to the given child to decide whether or not it will be
         * included in the filtered node.
         *
         * \tparam Child     The type of the child node.
         * \tparam new_index The index this child would receive in the filtered node.
         * \tparam old_index The index of this child in the unfiltered node.
         */
        template<typename Child, std::size_t new_index, std::size_t old_index>
        struct apply
        {
          //! True if the child will be included in the filtered node.
          static const bool value = true;
        };

      };

      namespace {

        // ********************************************************************************
        // IndexFilter helpers
        // ********************************************************************************

        template<typename Node, std::size_t new_index, std::size_t... indices>
        struct index_filter_helper
        {
          template<typename... FilterEntries>
          struct apply
          {
            typedef FilterResult<FilterEntries...> type;
          };
        };

        template<typename Node, std::size_t new_index, std::size_t old_index, std::size_t... indices>
        struct index_filter_helper<Node,new_index,old_index,indices...>
        {
          template<typename... FilterEntries>
          struct apply
            : public index_filter_helper<Node,new_index+1,indices...>::template apply<FilterEntries...,
                                                                                      FilterEntry<new_index,
                                                                                                  old_index>
                                                                                      >
          {};
        };

      } // anonymous namespace


      //! Filter class for FilteredCompositeNode that selects the children with the given indices.
      template<std::size_t... indices>
      struct IndexFilter
        : public AdvancedFilter
      {

#ifndef DOXYGEN

        template<typename Node, typename... Children>
        struct apply
        {
          typedef typename index_filter_helper<Node,0,indices...>::template apply<>::type type;
        };

#endif // DOXYGEN

      };


      // ********************************************************************************
      // filter: Wrapper class for turning a simple filter into an advanced filter
      //         usable by FilteredCompositeNode
      // ********************************************************************************

      namespace {

        template<typename Filter, std::size_t new_k, std::size_t old_k, typename... tail>
        struct filter_helper
        {
          template<typename... FilterDescriptors>
          struct apply
          {
            typedef FilterResult<FilterDescriptors...> type;
          };
        };

        template<typename Filter, std::size_t new_k, std::size_t old_k, typename child, typename... tail>
        struct filter_helper<Filter,new_k,old_k,child,tail...>
        {

          template<typename... FilterDescriptors>
          struct apply
            : public SelectType<Filter::template apply<child,new_k,old_k>::value,
                                typename filter_helper<Filter,new_k+1,old_k+1,tail...>::template apply<FilterDescriptors...,FilterEntry<new_k,old_k> >,
                                typename filter_helper<Filter,new_k,old_k+1,tail...>::template apply<FilterDescriptors...>
                                >::Type
          {};

        };

      } // anonymous namespace

      //! Adapter class that takes a SimpleFilter, validated it and turns it into an AdvancedFilter.
      template<typename Filter>
      struct filter
      {

        //! Apply the filter.
        template<typename Node, typename... Children>
        struct apply
        {

          dune_static_assert((Filter::template validate<Node>::value),"Invalid simple filter");

          typedef typename filter_helper<Filter,0,0,Children...>::template apply<>::type type;

        };

      };

      //! \} group Nodes

    } // namespace TypeTree

  } // namespace PDELab
} //namespace Dune

#endif // (HAVE_VARIADIC_TEMPLATES || DOXYGEN)

#endif // DUNE_PDELAB_COMMON_TYPETREE_FILTERS_HH
