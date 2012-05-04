// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_CHILDEXTRACTION_HH
#define DUNE_PDELAB_COMMON_TYPETREE_CHILDEXTRACTION_HH

#include <dune/common/documentation.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/pdelab/common/typetree/treepath.hh>


namespace Dune {
  namespace PDELab {

  namespace TypeTree {


    //! \addtogroup TypeTreeChildExtraction Child Extraction
    //! Utility functions and metafunctions for extracting children from a TypeTree.
    //! \ingroup TypeTree
    //! \{

    //! Extract the type of the child of Node at position TreePath.
    template<typename Node, typename TreePath>
    struct extract_child_type
    {

      //! The type of the child.
      typedef typename extract_child_type<
        typename Node::template Child<TypeTree::TreePathFront<TreePath>::value>::Type,
        typename TypeTree::TreePathPopFront<TreePath>::type
        >::type type;

      //! The storage type of the child.
      typedef typename extract_child_type<
        typename Node::template Child<TypeTree::TreePathFront<TreePath>::value>::Type,
        typename TypeTree::TreePathPopFront<TreePath>::type
        >::storage_type storage_type;

      //! The const storage type of the child.
      typedef typename extract_child_type<
        typename Node::template Child<TypeTree::TreePathFront<TreePath>::value>::Type,
        typename TypeTree::TreePathPopFront<TreePath>::type
        >::const_storage_type const_storage_type;


    };

#ifndef DOXYGEN

#if HAVE_VARIADIC_TEMPLATES

    // end of recursion
    template<typename Node>
    struct extract_child_type<Node,TypeTree::TreePath<> >
    {
      typedef Node type;
      typedef shared_ptr<Node> storage_type;
      typedef shared_ptr<const Node> const_storage_type;
    };

#else // HAVE_VARIADIC_TEMPLATES


    // end of recursion
    template<typename Node>
    struct extract_child_type<Node,
                              TypeTree::TreePath<
                                noChildIndex,noChildIndex,
                                noChildIndex,noChildIndex,
                                noChildIndex,noChildIndex,
                                noChildIndex,noChildIndex,
                                noChildIndex,noChildIndex
                                >
                              >
    {
      typedef Node type;
      typedef shared_ptr<Node> storage_type;
      typedef shared_ptr<const Node> const_storage_type;
    };

#endif // HAVE_VARIADIC_TEMPLATES

#endif // DOXYGEN



#ifdef DOXYGEN

    //! Extract the child of a node located at tp (non-const version).
    /**
     * Use this function to extract a (possibly indirect) child of
     * a TypeTree node.
     *
     * Example:
     *
     * extract_child(node,Dune::PDELab::TypeTree::TreePath<2,3,0>())
     *
     * returns the first child of the fourth child of the third child
     * of node.
     *
     * \sa Use extract_child_type to determine the type of the return
     *     value.
     *
     * \param node      The node from which to extract the child.
     * \param tp        The path into the tree leading to the child.
     *                  Note that the actual instance is not used
     *                  at all by this function, only the type of
     *                  the parameter.
     * \tparam TreePath A TreePath instantiation which statically
     *                  encodes the path to the child.
     * \return          A reference to the child.
     */
    template<typename Node, typename TreePath>
    ImplementationDefined& extract_child(Node& node, Treepath tp)
    {}

    //! Extract the child of a node located at tp (const version).
    /**
     * Use this function to extract a (possibly indirect) child of
     * a TypeTree node.
     *
     * Example:
     *
     * extract_child(node,Dune::PDELab::TypeTree::TreePath<2,3,0>())
     *
     * returns the first child of the fourth child of the third child
     * of node.
     *
     * \sa Use extract_child_type to determine the type of the return
     *     value.
     *
     * \param node      The node from which to extract the child.
     * \param tp        The path into the tree leading to the child.
     *                  Note that the actual instance is not used
     *                  at all by this function, only the type of
     *                  the parameter.
     * \tparam TreePath A TreePath instantiation which statically
     *                  encodes the path to the child.
     * \return          A reference to the child.
     */
    template<typename Node, typename TreePath>
    const ImplementationDefined& extract_child(const Node& node, Treepath tp)
    {}

#else // DOXYGEN

    // ********************************************************************************
    // non-const implementation
    // ********************************************************************************

    template<typename Node, typename TreePath>
    typename enable_if<
      (TypeTree::TreePathSize<TreePath>::value > 1),
      typename extract_child_type<Node,TreePath>::type&
      >::type
    extract_child(Node& node, TreePath tp)
    {
      return extract_child(node.template child<TypeTree::TreePathFront<TreePath>::value>(),
                           typename TypeTree::TreePathPopFront<TreePath>::type());
    }

    template<typename Node, typename TreePath>
    typename enable_if<
      TypeTree::TreePathSize<TreePath>::value == 1,
      typename Node::template Child<TypeTree::TreePathFront<TreePath>::value>::Type&
      >::type
    extract_child(Node& node, TreePath tp)
    {
      return node.template child<TypeTree::TreePathFront<TreePath>::value>();
    }

    template<typename Node, typename TreePath>
    typename enable_if<
      TypeTree::TreePathSize<TreePath>::value == 0,
      Node&
      >::type
    extract_child(Node& node, TreePath tp)
    {
      return node;
    }

    // ********************************************************************************
    // const implementation
    // ********************************************************************************

    template<typename Node, typename TreePath>
    typename enable_if<
      (TypeTree::TreePathSize<TreePath>::value > 1),
      const typename extract_child_type<Node,TreePath>::type&
      >::type
    extract_child(const Node& node, TreePath tp)
    {
      return extract_child(node.template child<TypeTree::TreePathFront<TreePath>::value>(),
                           typename TypeTree::TreePathPopFront<TreePath>::type());
    }

    template<typename Node, typename TreePath>
    typename enable_if<
      TypeTree::TreePathSize<TreePath>::value == 1,
      const typename Node::template Child<TypeTree::TreePathFront<TreePath>::value>::Type&
      >::type
    extract_child(const Node& node, TreePath tp)
    {
      return node.template child<TypeTree::TreePathFront<TreePath>::value>();
    }

    template<typename Node, typename TreePath>
    typename enable_if<
      TypeTree::TreePathSize<TreePath>::value == 0,
      const Node&
      >::type
    extract_child(const Node& node, TreePath tp)
    {
      return node;
    }


#endif // DOXYGEN



#ifdef DOXYGEN

    //! Extract the storage for the child of a node located at tp
    //! (non-const version).
    /**
     * Use this function to extract the storage (usually a shared_ptr)
     * of a (possibly indirect) child of a TypeTree node.
     *
     * Example:
     *
     * extract_child_storage(node,Dune::PDELab::TypeTree::TreePath<2,3,0>())
     *
     * returns the first child of the fourth child of the third child
     * of node.
     *
     * \sa Use extract_child_type to determine the type of the return
     *     value.
     *
     * \param node      The node from which to extract the child.
     * \param tp        The path into the tree leading to the child.
     *                  Note that the actual instance is not used
     *                  at all by this function, only the type of
     *                  the parameter.
     * \tparam TreePath A TreePath instantiation which statically
     *                  encodes the path to the child.
     * \return          A reference to the child.
     */
    template<typename Node, typename TreePath>
    ImplementationDefined extract_child_storage(Node& node, Treepath tp)
    {}

    //! Extract the storage for the child of a node located at tp
    //! (const version).
    /**
     * Use this function to extract the const storage (usually a shared_ptr)
     * of a (possibly indirect) child of a TypeTree node.
     *
     * Example:
     *
     * extract_child_storage(node,Dune::PDELab::TypeTree::TreePath<2,3,0>())
     *
     * returns the first child of the foruth child of the third child
     * of node.
     *
     * \sa Use extract_child_type to determine the type of the return
     *     value.
     *
     * \param node      The node from which to extract the child.
     * \param tp        The path into the tree leading to the child.
     *                  Note that the actual instance is not used
     *                  at all by this function, only the type of
     *                  the parameter.
     * \tparam TreePath A TreePath instantiation which statically
     *                  encodes the path to the child.
     * \return          A reference to the child.
     */
    template<typename Node, typename TreePath>
    ImplementationDefined extract_child_storage(const Node& node, Treepath tp)
    {}

#else // DOXYGEN

    // ********************************************************************************
    // non-const implementation
    // ********************************************************************************

    template<typename Node, typename TreePath>
    typename enable_if<
      (TypeTree::TreePathSize<TreePath>::value > 1),
      typename extract_child_type<Node,TreePath>::storage_type
      >::type
    extract_child_storage(Node& node, TreePath tp)
    {
      return extract_child_storage(node.template child<TypeTree::TreePathFront<TreePath>::value>(),
                                   typename TypeTree::TreePathPopFront<TreePath>::type());
    }

    template<typename Node, typename TreePath>
    typename enable_if<
      TypeTree::TreePathSize<TreePath>::value == 1,
      typename Node::template Child<TypeTree::TreePathFront<TreePath>::value>::Storage&
      >::type
    extract_child_storage(Node& node, TreePath tp)
    {
      return node.template childStorage<TypeTree::TreePathFront<TreePath>::value>();
    }

    template<typename Node, typename TreePath>
    typename enable_if<
      TypeTree::TreePathSize<TreePath>::value == 0
      >::type
    extract_child_storage(Node& node, TreePath tp)
    {
      dune_static_assert((Dune::AlwaysFalse<Node>::value),
                          "extract_child_storage only works for real children, not the node itself.");
    }

    // ********************************************************************************
    // const implementation
    // ********************************************************************************

    template<typename Node, typename TreePath>
    typename enable_if<
      (TypeTree::TreePathSize<TreePath>::value > 1),
      typename extract_child_type<Node,TreePath>::const_storage_type
      >::type
    extract_child_storage(const Node& node, TreePath tp)
    {
      return extract_child_storage(node.template child<TypeTree::TreePathFront<TreePath>::value>(),
                                   typename TypeTree::TreePathPopFront<TreePath>::type());
    }

    template<typename Node, typename TreePath>
    typename enable_if<
      TypeTree::TreePathSize<TreePath>::value == 1,
      const typename Node::template Child<TypeTree::TreePathFront<TreePath>::value>::ConstStorage&
      >::type
    extract_child_storage(const Node& node, TreePath tp)
    {
      return node.template childStorage<TypeTree::TreePathFront<TreePath>::value>();
    }

    template<typename Node, typename TreePath>
    typename enable_if<
      TypeTree::TreePathSize<TreePath>::value == 0
      >::type
    extract_child_storage(const Node& node, TreePath tp)
    {
      dune_static_assert((Dune::AlwaysFalse<Node>::value),
                          "extract_child_storage only works for real children, not the node itself.");
    }


#endif // DOXYGEN


    //! \} group TypeTree

  } // namespace TypeTree

  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_CHILDEXTRACTION_HH
