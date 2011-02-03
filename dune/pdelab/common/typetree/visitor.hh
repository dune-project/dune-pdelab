// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_VISITOR_HH
#define DUNE_PDELAB_COMMON_TYPETREE_VISITOR_HH

#include <dune/pdelab/common/typetree/treepath.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Tree Traversal
       *  \ingroup TypeTree
       *  \{
       */

      //! Visitor interface and base class for TypeTree visitors.
      /**
       * DefaultVisitor defines the interface for visitors that can be applied to a TypeTree
       * using applyToTree(). Each method of the visitor is passed a node of the tree (either as
       * a mutable or a const reference, depending on the constness of the tree applyToTree() was
       * called with). The second argument is of type TreePath and denotes the exact position of the
       * node within the TypeTree, encoded as child indices starting at the root node.
       *
       * \note This class can also be used as a convenient base class if the implemented visitor
       * only needs to act on some of the possible callback sites, avoiding a lot of boilerplate code.
       */
      struct DefaultVisitor
      {

#if HAVE_RVALUE_REFERENCES || DOXYGEN

        //! Method for prefix tree traversal.
        /**
         * This method gets called when first encountering a non-leaf node and
         * before visiting any of its children.
         *
         * \param t        The node to visit.
         * \param treePath The position of the node within the TypeTree.
         */
        template<typename T, typename TreePath>
        void pre(T&& t, TreePath treePath) const {}

        //! Method for infix tree traversal.
        /**
         * This method gets called BETWEEN visits of children of a non-leaf node.
         * That definition implies that this method will only be called for nodes
         * with at least two children.
         *
         * \param t        The node to visit.
         * \param treePath The position of the node within the TypeTree.
         */
        template<typename T, typename TreePath>
        void in(T&& t, TreePath treePath) const {}

        //! Method for postfix tree traversal.
        /**
         * This method gets called after all children of a non-leaf node have
         * been visited.
         *
         * \param t        The node to visit.
         * \param treePath The position of the node within the TypeTree.
         */
        template<typename T, typename TreePath>
        void post(T&& t, TreePath treePath) const {}

        //! Method for leaf traversal.
        /**
         * This method gets called when encountering a leaf-node within the TypeTree.
         *
         * \param t        The node to visit.
         * \param treePath The position of the node within the TypeTree.
         */
        template<typename T, typename TreePath>
        void leaf(T&& t, TreePath treePath) const {}

        //! Method for parent-child traversal.
        /**
         * This method gets called before visiting a child node.
         *
         * \note This method gets called even if the visitor decides not to visit the child in question.
         *
         * \param t          The parent node.
         * \param Child      The child node that will (potentially) be visited next.
         * \param treePath   The position of the parent node within the TypeTree.
         * \param childIndex The index of the child node in relation to the parent node.
         */
        template<typename T, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(T&& t, Child&& child, TreePath treePath, ChildIndex childIndex) const {}

        //! Method for parent-child traversal.
        /**
         * This method gets called after visiting a child node.
         *
         * \note This method gets called even if the child node was not visited because the visitor
         *       chose not to do so.
         *
         * \param t          The parent node.
         * \param Child      The child node that was visited last (if the visitor did not reject it).
         * \param treePath   The position of the parent node within the TypeTree.
         * \param childIndex The index of the child node in relation to the parent node.
         */
        template<typename T, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(T&& t, Child&& child, TreePath treePath, ChildIndex childIndex) const {}

#else // HAVE_RVALUE_REFERENCES

        // These are just a repeat of the above if the compiler does not support
        // rvalue references. Since our methods are all empty, we only need methods
        // which take the node as const reference.

        // Method for prefix traversal
        template<typename T, typename TreePath>
        void pre(const T& t, TreePath treePath) const {}

        // Method for infix traversal
        template<typename T, typename TreePath>
        void in(const T& t, TreePath treePath) const {}

        // Method for postfix traversal
        template<typename T, typename TreePath>
        void post(const T& t, TreePath treePath) const {}

        // Method for leaf traversal
        template<typename T, typename TreePath>
        void leaf(const T& t, TreePath treePath) const {}

        // Method for parent-child traversal
        template<typename T, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(const T& t, const Child& child, TreePath treePath, ChildIndex childIndex) const {}

        // Method for child-parent traversal
        template<typename T, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(const T& t, const Child& child, TreePath treePath, ChildIndex childIndex) const {}

#endif // HAVE_RVALUE_REFERENCES || DOXYGEN

      };


      //! Visitor interface and base class for visitors of pairs of TypeTrees.
      /**
       * DefaultPairVisitor defines the interface for visitors that can be applied to a TypeTree
       * using applyToTree(). Each method of the visitor is passed a node of the tree (either as
       * a mutable or a const reference, depending on the constness of the tree applyToTree() was
       * called with). The second argument is of type TreePath and denotes the exact position of the
       * node within the TypeTree, encoded as child indices starting at the root node.
       *
       * \note This class can also be used as a convenient base class if the implemented visitor
       * only needs to act on some of the possible callback sites, avoiding a lot of boilerplate code.
       */
      struct DefaultPairVisitor
      {

#if HAVE_RVALUE_REFERENCES || DOXYGEN

        //! Method for prefix tree traversal.
        /**
         * This method gets called when first encountering a non-leaf node and
         * before visiting any of its children.
         *
         * \param t        The node to visit.
         * \param treePath The position of the node within the TypeTree.
         */
        template<typename T1, typename T2, typename TreePath>
        void pre(T1&& t1, T2&& t2, TreePath treePath) const {}

        //! Method for infix tree traversal.
        /**
         * This method gets called BETWEEN visits of children of a non-leaf node.
         * That definition implies that this method will only be called for nodes
         * with at least two children.
         *
         * \param t        The node to visit.
         * \param treePath The position of the node within the TypeTree.
         */
        template<typename T1, typename T2, typename TreePath>
        void in(T1&& t1, T2&& t2, TreePath treePath) const {}

        //! Method for postfix traversal.
        /**
         * This method gets called after all children of a non-leaf node have
         * been visited.
         *
         * \param t        The node to visit.
         * \param treePath The position of the node within the TypeTree.
         */
        template<typename T1, typename T2, typename TreePath>
        void post(T1&& t1, T2&& t2, TreePath treePath) const {}

        //! Method for leaf traversal.
        /**
         * This method gets called when encountering a leaf-node within the TypeTree.
         *
         * \param t        The node to visit.
         * \param treePath The position of the node within the TypeTree.
         */
        template<typename T1, typename T2, typename TreePath>
        void leaf(T1&& t1, T2&& t2, TreePath treePath) const {}

        //! Method for parent-child traversal.
        /**
         * This method gets called before visiting a child node.
         *
         * \note This method gets called even if the visitor decides not to visit the child in question.
         *
         * \param t        The node to visit.
         * \param treePath The position of the node within the TypeTree.
         */
        template<typename T1, typename Child1, typename T2, typename Child2, typename TreePath, typename ChildIndex>
        void beforeChild(T1&& t1, Child1&& child1, T2&& t2, Child2&& child2, TreePath treePath, ChildIndex childIndex) const {}

        template<typename T1, typename Child1, typename T2, typename Child2, typename TreePath, typename ChildIndex>
        void afterChild(T1&& t1, Child1&& child1, T2&& t2, Child2&& child2, TreePath treePath, ChildIndex childIndex) const {}

#else // HAVE_RVALUE_REFERENCES

        // These are just a repeat of the above if the compiler does not support
        // rvalue references. In this case, we need variants for const and non-const
        // nodes.

        // Method for prefix traversal
        template<typename T1, typename T2, typename TreePath>
        void pre(T1& t1, T2& t2, TreePath treePath) const {}

        // Method for infix traversal
        template<typename T1, typename T2, typename TreePath>
        void in(T1& t1, T2& t2, TreePath treePath) const {}

        // Method for postfix traversal
        template<typename T1, typename T2, typename TreePath>
        void post(T1& t1, T2& t2, TreePath treePath) const {}

        // Method for leaf traversal
        template<typename T1, typename T2, typename TreePath>
        void leaf(T1& t1, T2& t2, TreePath treePath) const {}

        template<typename T1, typename Child1, typename T2, typename Child2, typename TreePath, typename ChildIndex>
        void beforeChild(T1& t1, Child1& child1,
                         T2& t2, Child2& child2,
                         TreePath treePath,
                         ChildIndex childIndex) const {}

        template<typename T1, typename Child1, typename T2, typename Child2, typename TreePath, typename ChildIndex>
        void afterChild(T1& t1, Child1& child1,
                        T2& t2, Child2& child2,
                        TreePath treePath,
                        ChildIndex childIndex) const {}

        // Method for prefix traversal
        template<typename T1, typename T2, typename TreePath>
        void pre(const T1& t1, const T2& t2, TreePath treePath) const {}

        // Method for infix traversal
        template<typename T1, typename T2, typename TreePath>
        void in(const T1& t1, const T2& t2, TreePath treePath) const {}

        // Method for postfix traversal
        template<typename T1, typename T2, typename TreePath>
        void post(const T1& t1, const T2& t2, TreePath treePath) const {}

        // Method for leaf traversal
        template<typename T1, typename T2, typename TreePath>
        void leaf(const T1& t1, const T2& t2, TreePath treePath) const {}

        template<typename T1, typename Child1, typename T2, typename Child2, typename TreePath, typename ChildIndex>
        void beforeChild(const T1& t1, const Child1& child1,
                         const T2& t2, const Child2& child2,
                         TreePath treePath,
                         ChildIndex childIndex) const {}

        template<typename T1, typename Child1, typename T2, typename Child2, typename TreePath, typename ChildIndex>
        void afterChild(const T1& t1, const Child1& child1,
                        const T2& t2, const Child2& child2,
                        TreePath treePath,
                        ChildIndex childIndex) const {}

#endif // HAVE_RVALUE_REFERENCES || DOXYGEN

      };

      //! Mixin base class for visitors that only want to visit the direct children of a node.
      /**
       * This mixin class will reject all children presented to it, causing the algorithm to
       * only visit the root node and call DefaultVisitor::beforeChild() and DefaultVisitor::afterChild()
       * for its direct children.
       */
      struct VisitDirectChildren
      {

        // the little trick with the default template arguments
        // makes the class usable for both single-tree visitors
        // and visitors for pairs of trees
        //! Template struct for determining whether or not to visit a given child.
        template<typename Node1,
                 typename Child1,
                 typename Node2,
                 typename Child2 = void,
                 typename TreePath = void>
        struct VisitChild
        {
          //! Do not visit any child.
          static const bool value = false;
        };

      };


      //! Mixin base class for visitors that want to visit the complete tree.
      /**
       * This mixin class will accept all children presented to it and thus make the iterator
       * traverse the entire tree.
       */
      struct VisitTree
      {

        // the little trick with the default template arguments
        // makes the class usable for both single-tree visitors
        // and visitors for pairs of trees
        //! Template struct for determining whether or not to visit a given child.
        template<typename Node1,
                 typename Child1,
                 typename Node2,
                 typename Child2 = void,
                 typename TreePath = void>
        struct VisitChild
        {
          //! Visit any child.
          static const bool value = true;
        };

      };

      //! Mixin base class for visitors that require a static TreePath during traversal.
      /**
       * \warning Static traversal should only be used if absolutely necessary, as it tends
       *          to increase compilation times and object sizes (especially if compiling
       *          with debug information)!
       *
       * \sa DynamicTraversal
       */
      struct StaticTraversal
      {
        //! Use the static tree traversal algorithm.
        static const TreePathType::Type treePathType = TreePathType::fullyStatic;
      };

      //! Mixin base class for visitors that only need a dynamic TreePath during traversal.
      /**
       * \note Dynamic traversal is preferable to static traversal, as it causes fewer
       *       template instantiations, which improves compile time and reduces object
       *       size (especially if compiling with debug information).
       *
       * \sa StaticTraversal
       */
      struct DynamicTraversal
      {
        //! Use the dynamic tree traversal algorithm.
        static const TreePathType::Type treePathType = TreePathType::dynamic;
      };

      //! Convenience base class for visiting the entire tree.
      struct TreeVisitor
        : public DefaultVisitor
        , public VisitTree
      {};

      //! Convenience base class for visiting the direct children of a node.
      struct DirectChildrenVisitor
        : public DefaultVisitor
        , public VisitDirectChildren
      {};

      //! Convenience base class for visiting an entire tree pair.
      struct TreePairVisitor
        : public DefaultPairVisitor
        , public VisitTree
      {};

      //! Convenience base class for visiting the direct children of a node pair.
      struct DirectChildrenPairVisitor
        : public DefaultPairVisitor
        , public VisitDirectChildren
      {};

      //! \} group Tree Traversal

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_VISITOR_HH
