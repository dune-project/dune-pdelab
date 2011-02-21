// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_PAIRTRAVERSAL_HH
#define DUNE_PDELAB_COMMON_TYPETREE_PAIRTRAVERSAL_HH

#include <dune/pdelab/common/typetree/nodetags.hh>
#include <dune/pdelab/common/typetree/treepath.hh>
#include <dune/pdelab/common/typetree/visitor.hh>

#if HAVE_RVALUE_REFERENCES
#include <utility>
#endif

#ifndef DOXYGEN

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      // forward declaration for included implementation parts
      template<TreePathType::Type tpType, typename tag1 = StartTag, typename tag2 = StartTag, bool doApply = true>
      struct ApplyToTreePair;

    }
  }
}

#endif // DOXYGEN

#include <dune/pdelab/common/typetree/applytochildrentreepair.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Tree Traversal
       *  \ingroup TypeTree
       *  \{
       */

#ifndef DOXYGEN // these are all internals and not public API. Only access is using applyToTree().

      template<TreePathType::Type tpType>
      struct ApplyToTreePair<tpType,StartTag,StartTag,true>
      {

#if HAVE_RVALUE_REFERENCES

        template<typename Node1, typename Node2, typename Visitor>
        static void apply(Node1&& node1, Node2&& node2, Visitor&& visitor)
        {
          ApplyToTreePair<tpType,
                          typename remove_reference<Node1>::type::NodeTag,
                          typename remove_reference<Node2>::type::NodeTag
                          >::apply(std::forward<Node1>(node1),
                                   std::forward<Node2>(node2),
                                   std::forward<Visitor>(visitor),
                                   TreePathFactory<tpType>::create().mutablePath());
        }

#else

        // The next two methods do some nasty trickery to make sure that both trees
        // are either const or non-const and that no mixed case can occur. For this
        // purpose, the enable_if on the first method makes sure that it will only
        // ever match if both trees are non-const, and the second method casts both
        // trees to const before passing them on.
        template<typename Node1, typename Node2, typename Visitor>
        static
        typename enable_if<!(IsConst<Node1>::Value || IsConst<Node2>::Value)>::type
        apply(Node1& node1, Node2& node2, Visitor& visitor)
        {
          ApplyToTreePair<tpType,
                          typename Node1::NodeTag,
                          typename Node2::NodeTag
                          >::apply(node1,
                                   node2,
                                   visitor,
                                   TreePathFactory<tpType>::create().mutablePath());
        }

        template<typename Node1, typename Node2, typename Visitor>
        static void apply(const Node1& node1, const Node2& node2, Visitor& visitor)
        {
          ApplyToTreePair<tpType,
                          typename Node1::NodeTag,
                          typename Node2::NodeTag
                          >::apply(const_cast<const Node1&>(node1), // see previous method
                                   const_cast<const Node2&>(node2), // for explanation
                                   visitor,
                                   TreePathFactory<tpType>::create().mutablePath());
        }


        // The next two methods do some nasty trickery to make sure that both trees
        // are either const or non-const and that no mixed case can occur. For this
        // purpose, the enable_if on the first method makes sure that it will only
        // ever match if both trees are non-const, and the second method casts both
        // trees to const before passing them on.
        template<typename Node1, typename Node2, typename Visitor>
        static
        typename enable_if<!(IsConst<Node1>::Value || IsConst<Node2>::Value)>::type
        apply(Node1& node1, Node2& node2, const Visitor& visitor)
        {
          ApplyToTreePair<tpType,
                          typename Node1::NodeTag,
                          typename Node2::NodeTag
                          >::apply(node1,
                                   node2,
                                   visitor,
                                   TreePathFactory<tpType>::create().mutablePath());
        }

        template<typename Node1, typename Node2, typename Visitor>
        static void apply(const Node1& node1, const Node2& node2, const Visitor& visitor)
        {
          ApplyToTreePair<tpType,
                          typename Node1::NodeTag,
                          typename Node2::NodeTag
                          >::apply(const_cast<const Node1&>(node1), // see previous method
                                   const_cast<const Node2&>(node2), // for explanation
                                   visitor,
                                   TreePathFactory<tpType>::create().mutablePath());
        }

#endif // HAVE_RVALUE_REFERENCES

      };


      // Do not visit nodes the visitor is not interested in
      template<TreePathType::Type tpType, typename Tag1, typename Tag2>
      struct ApplyToTreePair<tpType,Tag1,Tag2,false>
      {
        template<typename Node1, typename Node2, typename Visitor, typename TreePath>
        static void apply(const Node1& node1, const Node2& node2, const Visitor& visitor, TreePath treePath)
        {}
      };


      /*

    // LeafNode - again, this is easy: just do all three visits
    template<TreePathType::Type tpType>
    struct ApplyToTree<tpType,LeafNodeTag,LeafNodeTag,true>
    {

#if HAVE_RVALUE_REFERENCES

      template<typename N1, typename N2, typename V, typename TreePath>
      static void apply(N1&& n1, N2&& n2, V&& v, TreePath tp)
      {
        v.leaf(std::forward<N1>(n1),std::forward<N2>(n2),tp.view());
      }

#else

      template<typename N1, typename N2, typename V, typename TreePath>
      static void apply(N1& n1, N2& n2, V& v, TreePath tp)
      {
        v.leaf(n1,n2,tp.view());
      }

      template<typename N1, typename N2, typename V, typename TreePath>
      static void apply(const N1& n1, const N2& n2, V& v, TreePath tp)
      {
        v.leaf(n1,n2,tp.view());
      }

      template<typename N1, typename N2, typename V, typename TreePath>
      static void apply(N1& n1, N2& n2, const V& v, TreePath tp)
      {
        v.leaf(n1,n2,tp.view());
      }

      template<typename N1, typename N2, typename V, typename TreePath>
      static void apply(const N1& n1, const N2& n2, const V& v, TreePath tp)
      {
        v.leaf(n1,n2,tp.view());
      }

#endif // HAVE_RVALUE_REFERENCES

    };
      */


      // Automatically pick the correct traversal algorithm for the two nodes
      template<TreePathType::Type treePathType,typename FirstTag, typename SecondTag>
      struct ApplyToTreePair<treePathType,FirstTag,SecondTag,true>
        : public ApplyToGenericCompositeNodePair<treePathType>
      {
      };



      // ********************************************************************************
      // Specialization for dynamic traversal and two PowerNodes -> use runtime iteration
      // ********************************************************************************

      template<>
      struct ApplyToTreePair<TreePathType::dynamic,PowerNodeTag,PowerNodeTag,true>
      {

#if HAVE_RVALUE_REFERENCES

        template<typename N1, typename N2, typename V, typename TreePath>
        static void apply(N1&& n1, N2&& n2, V&& v, TreePath tp)
        {
          v.pre(std::forward<N1>(n1),std::forward<N2>(n2),tp.view());
          typedef typename remove_reference<N1>::type Node1;
          typedef typename remove_reference<N2>::type Node2;
          typedef typename Node1::template Child<0>::Type C1;
          typedef typename Node2::template Child<0>::Type C2;
          dune_static_assert(Node1::CHILDREN == Node2::CHILDREN,
                             "non-leaf nodes with different numbers of children " \
                             "are not allowed during simultaneous grid traversal");
          const bool visit = std::remove_reference<V>::type
            ::template VisitChild<Node1,C1,Node2,C2,typename TreePath::ViewType>::value;
          for (std::size_t k = 0; k < Node1::CHILDREN; ++k)
            {
              v.beforeChild(std::forward<N1>(n1),n1.child(k),std::forward<N2>(n2),n2.child(k),tp.view(),k);
              tp.push_back(k);
              ApplyToTreePair<TreePathType::dynamic, // we know that due to the specialization
                              typename C1::NodeTag,
                              typename C2::NodeTag,
                              visit>::apply(n1.child(k),
                                            n2.child(k),
                                            std::forward<V>(v),
                                            tp);
              tp.pop_back();
              v.afterChild(std::forward<N1>(n1),n1.child(k),std::forward<N2>(n2),n2.child(k),tp.view(),k);
              if (k < Node1::CHILDREN-1)
                v.in(std::forward<N1>(n1),std::forward<N2>(n2),tp.view());
            }
          v.post(std::forward<N1>(n1),std::forward<N2>(n2),tp.view());
        }

#else

        template<typename N1, typename N2, typename V, typename TreePath>
        static void apply(N1& n1, N2& n2, V& v, TreePath tp)
        {
          v.pre(n1,n2,tp.view());
          typedef typename N1::template Child<0>::Type C1;
          typedef typename N2::template Child<0>::Type C2;
          dune_static_assert(N1::CHILDREN == N2::CHILDREN,
                             "non-leaf nodes with different numbers of children " \
                             "are not allowed during simultaneous grid traversal");
          const bool visit = V::template VisitChild<N1,C1,N2,C2,typename TreePath::ViewType>::value;
          for (std::size_t k = 0; k < N1::CHILDREN; ++k)
            {
              v.beforeChild(n1,n1.child(k),n2,n2.child(k),tp.view(),k);
              tp.push_back(k);
              ApplyToTreePair<TreePathType::dynamic, // we know that due to the specialization
                              typename C1::NodeTag,
                              typename C2::NodeTag,
                              visit>::apply(n1.child(k),
                                            n2.child(k),
                                            v,
                                            tp);
              tp.pop_back();
              v.afterChild(n1,n1.child(k),n2,n2.child(k),tp.view(),k);
              if (k < N1::CHILDREN-1)
                v.in(n1,n2,tp.view());
            }
          v.post(n1,n2,tp.view());
        }

        template<typename N1, typename N2, typename V, typename TreePath>
        static void apply(const N1& n1, const N2& n2, V& v, TreePath tp)
        {
          v.pre(n1,n2,tp.view());
          typedef typename N1::template Child<0>::Type C1;
          typedef typename N2::template Child<0>::Type C2;
          dune_static_assert(N1::CHILDREN == N2::CHILDREN,
                             "non-leaf nodes with different numbers of children " \
                             "are not allowed during simultaneous grid traversal");
          const bool visit = V::template VisitChild<N1,C1,N2,C2,typename TreePath::ViewType>::value;
          for (std::size_t k = 0; k < N1::CHILDREN; ++k)
            {
              v.beforeChild(n1,n1.child(k),n2,n2.child(k),tp.view(),k);
              tp.push_back(k);
              ApplyToTreePair<TreePathType::dynamic, // we know that due to the specialization
                              typename C1::NodeTag,
                              typename C2::NodeTag,
                              visit>::apply(n1.child(k),
                                            n2.child(k),
                                            v,
                                            tp);
              tp.pop_back();
              v.afterChild(n1,n1.child(k),n2,n2.child(k),tp.view(),k);
              if (k < N1::CHILDREN-1)
                v.in(n1,n2,tp.view());
            }
          v.post(n1,n2,tp.view());
        }

        template<typename N1, typename N2, typename V, typename TreePath>
        static void apply(N1& n1, N2& n2, const V& v, TreePath tp)
        {
          v.pre(n1,n2,tp.view());
          typedef typename N1::template Child<0>::Type C1;
          typedef typename N2::template Child<0>::Type C2;
          dune_static_assert(N1::CHILDREN == N2::CHILDREN,
                             "non-leaf nodes with different numbers of children " \
                             "are not allowed during simultaneous grid traversal");
          const bool visit = V::template VisitChild<N1,C1,N2,C2,typename TreePath::ViewType>::value;
          for (std::size_t k = 0; k < N1::CHILDREN; ++k)
            {
              v.beforeChild(n1,n1.child(k),n2,n2.child(k),tp.view(),k);
              tp.push_back(k);
              ApplyToTreePair<TreePathType::dynamic, // we know that due to the specialization
                              typename C1::NodeTag,
                              typename C2::NodeTag,
                              visit>::apply(n1.child(k),
                                            n2.child(k),
                                            v,
                                            tp);
              tp.pop_back();
              v.afterChild(n1,n1.child(k),n2,n2.child(k),tp.view(),k);
              if (k < N1::CHILDREN-1)
                v.in(n1,n2,tp.view());
            }
          v.post(n1,n2,tp.view());
        }

        template<typename N1, typename N2, typename V, typename TreePath>
        static void apply(const N1& n1, const N2& n2, const V& v, TreePath tp)
        {
          v.pre(n1,n2,tp.view());
          typedef typename N1::template Child<0>::Type C1;
          typedef typename N2::template Child<0>::Type C2;
          dune_static_assert(N1::CHILDREN == N2::CHILDREN,
                             "non-leaf nodes with different numbers of children " \
                             "are not allowed during simultaneous grid traversal");
          const bool visit = V::template VisitChild<N1,C1,N2,C2,typename TreePath::ViewType>::value;
          for (std::size_t k = 0; k < N1::CHILDREN; ++k)
            {
              v.beforeChild(n1,n1.child(k),n2,n2.child(k),tp.view(),k);
              tp.push_back(k);
              ApplyToTreePair<TreePathType::dynamic, // we know that due to the specialization
                              typename C1::NodeTag,
                              typename C2::NodeTag,
                              visit>::apply(n1.child(k),
                                            n2.child(k),
                                            v,
                                            tp);
              tp.pop_back();
              v.afterChild(n1,n1.child(k),n2,n2.child(k),tp.view(),k);
              if (k < N1::CHILDREN-1)
                v.in(n1,n2,tp.view());
            }
          v.post(n1,n2,tp.view());
        }

#endif // HAVE_RVALUE_REFERENCES

      };

#endif // DOXYGEN

#if HAVE_RVALUE_REFERENCES || DOXYGEN

      //! Apply visitor to a pair of TypeTrees.
      /**
       * This function applies the given visitor to the given tree. Both visitor and tree may be const
       * or non-const. If the compiler supports rvalue references, they may even be a non-const temporary;
       * otherwise both trees must be either const or non-const. If they have different constness, both will
       * be promoted to const.
       *
       * \note The visitor must implement the interface laid out by DefaultPairVisitor (most easily achieved by
       *       inheriting from it) and specify the required type of tree traversal (static or dynamic) by
       *       inheriting from either StaticTraversal or DynamicTraversal.
       *
       * \param tree1   The first tree the visitor will be applied to.
       * \param tree2   The second tree the visitor will be applied to.
       * \param visitor The visitor to apply to the trees.
       */
      template<typename Tree1, typename Tree2, typename Visitor>
      void applyToTreePair(Tree1&& tree1, Tree2&& tree2, Visitor&& visitor)
      {
        ApplyToTreePair<std::remove_reference<Visitor>::type::treePathType>::apply(std::forward<Tree1>(tree1),
                                                                                   std::forward<Tree2>(tree2),
                                                                                   std::forward<Visitor>(visitor));
      }

#else

      template<typename Tree1, typename Tree2, typename Visitor>
      void applyToTreePair(Tree1& tree1, Tree2& tree2, Visitor& visitor)
      {
        ApplyToTreePair<Visitor::treePathType>::apply(tree1,tree2,visitor);
      }

      template<typename Tree1, typename Tree2, typename Visitor>
      void applyToTreePair(const Tree1& tree1, const Tree2& tree2, Visitor& visitor)
      {
        ApplyToTreePair<Visitor::treePathType>::apply(tree1,tree2,visitor);
      }

      template<typename Tree1, typename Tree2, typename Visitor>
      void applyToTreePair(Tree1& tree1, Tree2& tree2, const Visitor& visitor)
      {
        ApplyToTreePair<Visitor::treePathType>::apply(tree1,tree2,visitor);
      }

      template<typename Tree1, typename Tree2, typename Visitor>
      void applyToTreePair(const Tree1& tree1, const Tree2& tree2, const Visitor& visitor)
      {
        ApplyToTreePair<Visitor::treePathType>::apply(tree1,tree2,visitor);
      }

#endif // HAVE_RVALUE_REFERENCES || DOXYGEN

      //! \} group Tree Traversal

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_PAIRTRAVERSAL_HH
