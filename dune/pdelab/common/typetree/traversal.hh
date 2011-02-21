// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_TRAVERSAL_HH
#define DUNE_PDELAB_COMMON_TYPETREE_TRAVERSAL_HH

#if HAVE_RVALUE_REFERENCES
#include <utility>
#endif

#include <dune/pdelab/common/typetree/nodetags.hh>
#include <dune/pdelab/common/typetree/treepath.hh>
#include <dune/pdelab/common/typetree/visitor.hh>

#ifndef DOXYGEN

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      // forward declaration of TMP struct for child visits - required for included implementation parts
      template<TreePathType::Type tpType, typename tag = StartTag, bool doApply = true>
      struct ApplyToTree;

    }
  }
}

#endif // DOXYGEN

#include <dune/pdelab/common/typetree/applytochildrensingletree.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Tree Traversal
       *  \ingroup TypeTree
       *  \{
       */


#ifndef DOXYGEN // these are all internals and not public API. Only access is using applyToTree().


      // This struct is the core of the algorithm. While this specialization simply serves as the starting point
      // of the traversal and takes care of some setup work, the struct has to be specialized for each TreeType node type it
      // should support.
      // The first parameter specifies the kind of TreePath (dynamic/static) to use, the second one is the tag of the node type
      // and the third one must always be specialized as true, as a value of false means the node should in fact not be visited.
      // That case is already handled by a specialization of the struct.
      template<TreePathType::Type tpType, bool doApply>
      struct ApplyToTree<tpType,StartTag,doApply>
      {

#if HAVE_RVALUE_REFERENCES

        template<typename Node, typename Visitor>
        static void apply(Node&& node, Visitor&& visitor)
        {
          ApplyToTree<tpType,typename remove_reference<Node>::type::NodeTag>::apply(std::forward<Node>(node),
                                                                                    std::forward<Visitor>(visitor),
                                                                                    TreePathFactory<tpType>::create().mutablePath());
        }

#else

        template<typename Node, typename Visitor>
        static void apply(Node& node, Visitor& visitor)
        {
          ApplyToTree<tpType,typename Node::NodeTag>::apply(node,visitor,TreePathFactory<tpType>::create().mutablePath());
        }

        template<typename Node, typename Visitor>
        static void apply(const Node& node, Visitor& visitor)
        {
          ApplyToTree<tpType,typename Node::NodeTag>::apply(node,visitor,TreePathFactory<tpType>::create().mutablePath());
        }

        template<typename Node, typename Visitor>
        static void apply(Node& node, const Visitor& visitor)
        {
          ApplyToTree<tpType,typename Node::NodeTag>::apply(node,visitor,TreePathFactory<tpType>::create().mutablePath());
        }

        template<typename Node, typename Visitor>
        static void apply(const Node& node, const Visitor& visitor)
        {
          ApplyToTree<tpType,typename Node::NodeTag>::apply(node,visitor,TreePathFactory<tpType>::create().mutablePath());
        }

#endif // HAVE_RVALUE_REFERENCES

      };


      // Do not visit nodes the visitor is not interested in
      template<TreePathType::Type tpType, typename NodeTag>
      struct ApplyToTree<tpType,NodeTag,false>
      {

        // we won't do anything with the objects, so having them all const
        // works fine.
        template<typename Node, typename Visitor, typename TreePath>
        static void apply(const Node& node, const Visitor& visitor, TreePath treePath)
        {}

      };



      // ********************************************************************************
      // LeafNode
      // ********************************************************************************

      // LeafNode - just call the leaf() callback
      template<TreePathType::Type tpType>
      struct ApplyToTree<tpType,LeafNodeTag,true>
      {

#if HAVE_RVALUE_REFERENCES

        template<typename N, typename V, typename TreePath>
        static void apply(N&& n, V&& v, TreePath tp)
        {
          v.leaf(std::forward<N>(n),tp.view());
        }

#else

        template<typename N, typename V, typename TreePath>
        static void apply(N& n, V& v, TreePath tp)
        {
          v.leaf(n,tp.view());
        }

        template<typename N, typename V, typename TreePath>
        static void apply(const N& n, V& v, TreePath tp)
        {
          v.leaf(n,tp.view());
        }

        template<typename N, typename V, typename TreePath>
        static void apply(N& n, const V& v, TreePath tp)
        {
          v.leaf(n,tp.view());
        }

        template<typename N, typename V, typename TreePath>
        static void apply(const N& n, const V& v, TreePath tp)
        {
          v.leaf(n,tp.view());
        }

#endif // HAVE_RVALUE_REFERENCES

      };



      // ********************************************************************************
      // PowerNode
      // ********************************************************************************

      // Traverse PowerNode statically - in this case, we simply use the
      // generic child traversal algorithm
      template<>
      struct ApplyToTree<TreePathType::fullyStatic,PowerNodeTag,true>
        : public ApplyToGenericCompositeNode
      {
      };

      // Traverse PowerNode dynamically. Here, we exploit the fact that is possible
      // to do the child traversal using runtime iteration, as that saves a lot of
      // template instantiations.
      template<>
      struct ApplyToTree<TreePathType::dynamic,PowerNodeTag,true>
      {

#if HAVE_RVALUE_REFERENCES

        template<typename N, typename V, typename TreePath>
        static void apply(N&& n, V&& v, TreePath tp)
        {
          // first encounter of this node
          v.pre(std::forward<N>(n),tp.view());

          // strip types of possible references
          typedef typename remove_reference<N>::type Node;
          typedef typename remove_reference<V>::type Visitor;

          // get child type
          typedef typename Node::template Child<0>::Type C;

          // Do we have to visit the children? As the TreePath is dynamic, it does not
          // contain any information that could be evaluated at compile time, so we only
          // have to query the visitor once.
          const bool visit = Visitor::template VisitChild<Node,C,typename TreePath::ViewType>::value;

          // iterate over children
          for (std::size_t k = 0; k < Node::CHILDREN; ++k)
            {
              // always call beforeChild(), regardless of the value of visit
              v.beforeChild(std::forward<N>(n),n.child(k),tp.view(),k);

              // update TreePath
              tp.push_back(k);

              // descend to child
              ApplyToTree<Visitor::treePathType,typename C::NodeTag,visit>::apply(n.child(k),std::forward<V>(v),tp);

              // restore TreePath
              tp.pop_back();

              // always call afterChild(), regardless of the value of visit
              v.afterChild(std::forward<N>(n),n.child(k),tp.view(),k);

              // if this is not the last child, call infix callback
              if (k < Node::CHILDREN-1)
                v.in(std::forward<N>(n),tp.view());
            }

          // node is done - call postfix callback
          v.post(std::forward<N>(n),tp.view());
        }

#else

        // non-const node, non-const visitor
        template<typename N, typename V, typename TreePath>
        static void apply(N& n, V& v, TreePath tp)
        {
          v.pre(n,tp.view());
          typedef typename N::template Child<0>::Type C;
          const bool visit = V::template VisitChild<N,C,typename TreePath::ViewType>::value;
          for (std::size_t k = 0; k < N::CHILDREN; ++k)
            {
              v.beforeChild(n,n.child(k),tp.view(),k);
              tp.push_back(k);
              ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.child(k),v,tp);
              tp.pop_back();
              v.afterChild(n,n.child(k),tp.view(),k);
              if (k < N::CHILDREN-1)
                v.in(n,tp.view());
            }
          v.post(n,tp.view());
        }

        // const node, non-const visitor
        template<typename N, typename V, typename TreePath>
        static void apply(const N& n, V& v, TreePath tp)
        {
          v.pre(n,tp.view());
          typedef typename N::template Child<0>::Type C;
          const bool visit = V::template VisitChild<N,C,typename TreePath::ViewType>::value;
          for (std::size_t k = 0; k < N::CHILDREN; ++k)
            {
              v.beforeChild(n,n.child(k),tp.view(),k);
              tp.push_back(k);
              ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.child(k),v,tp);
              tp.pop_back();
              v.afterChild(n,n.child(k),tp.view(),k);
              if (k < N::CHILDREN-1)
                v.in(n,tp.view());
            }
          v.post(n,tp.view());
        }

        // non-const node, const visitor
        template<typename N, typename V, typename TreePath>
        static void apply(N& n, const V& v, TreePath tp)
        {
          v.pre(n,tp.view());
          typedef typename N::template Child<0>::Type C;
          const bool visit = V::template VisitChild<N,C,typename TreePath::ViewType>::value;
          for (std::size_t k = 0; k < N::CHILDREN; ++k)
            {
              v.beforeChild(n,n.child(k),tp.view(),k);
              tp.push_back(k);
              ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.child(k),v,tp);
              tp.pop_back();
              v.afterChild(n,n.child(k),tp.view(),k);
              if (k < N::CHILDREN-1)
                v.in(n,tp.view());
            }
          v.post(n,tp.view());
        }

        // const node, const visitor
        template<typename N, typename V, typename TreePath>
        static void apply(const N& n, const V& v, TreePath tp)
        {
          v.pre(n,tp.view());
          typedef typename N::template Child<0>::Type C;
          const bool visit = V::template VisitChild<N,C,typename TreePath::ViewType>::value;
          for (std::size_t k = 0; k < N::CHILDREN; ++k)
            {
              v.beforeChild(n,n.child(k),tp.view(),k);
              tp.push_back(k);
              ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.child(k),v,tp);
              tp.pop_back();
              v.afterChild(n,n.child(k),tp.view(),k);
              if (k < N::CHILDREN-1)
                v.in(n,tp.view());
            }
          v.post(n,tp.view());
        }

#endif // HAVE_RVALUE_REFERENCES

      };



      // ********************************************************************************
      // CompositeNode, VariadicCompositeNode
      // ********************************************************************************

      // Traverse CompositeNode - just forward to the generic algorithm
      template<TreePathType::Type treePathType>
      struct ApplyToTree<treePathType,CompositeNodeTag,true>
        : public ApplyToGenericCompositeNode
      {
      };


      // Traverse VariadicCompositeNode - just forward to the generic algorithm
      template<TreePathType::Type treePathType>
      struct ApplyToTree<treePathType,VariadicCompositeNodeTag,true>
        : public ApplyToGenericCompositeNode
      {
      };

#endif // DOXYGEN



      // ********************************************************************************
      // Public Interface
      // ********************************************************************************

#if HAVE_RVALUE_REFERENCES || DOXYGEN

      //! Apply visitor to TypeTree.
      /**
       * This function applies the given visitor to the given tree. Both visitor and tree may be const
       * or non-const (if the compiler supports rvalue references, they may even be a non-const temporary).
       *
       * \note The visitor must implement the interface laid out by DefaultVisitor (most easily achieved by
       *       inheriting from it) and specify the required type of tree traversal (static or dynamic) by
       *       inheriting from either StaticTraversal or DynamicTraversal.
       *
       * \param tree    The tree the visitor will be applied to.
       * \param visitor The visitor to apply to the tree.
       */
      template<typename Tree, typename Visitor>
      void applyToTree(Tree&& tree, Visitor&& visitor)
      {
        ApplyToTree<std::remove_reference<Visitor>::type::treePathType>::apply(std::forward<Tree>(tree),
                                                                               std::forward<Visitor>(visitor));
      }

#else

      template<typename Tree, typename Visitor>
      void applyToTree(Tree& tree, Visitor& visitor)
      {
        ApplyToTree<Visitor::treePathType>::apply(tree,visitor);
      }

      template<typename Tree, typename Visitor>
      void applyToTree(const Tree& tree, Visitor& visitor)
      {
        ApplyToTree<Visitor::treePathType>::apply(tree,visitor);
      }

      template<typename Tree, typename Visitor>
      void applyToTree(Tree& tree, const Visitor& visitor)
      {
        ApplyToTree<Visitor::treePathType>::apply(tree,visitor);
      }

      template<typename Tree, typename Visitor>
      void applyToTree(const Tree& tree, const Visitor& visitor)
      {
        ApplyToTree<Visitor::treePathType>::apply(tree,visitor);
      }

#endif // HAVE_RVALUE_REFERENCES || DOXYGEN

      //! \} group Tree Traversal

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_TRAVERSAL_HH
