// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_TRAVERSAL_HH
#define DUNE_PDELAB_COMMON_TYPETREE_TRAVERSAL_HH

#include <dune/pdelab/common/typetree/nodetags.hh>
#include <dune/pdelab/common/typetree/treepath.hh>

#if HAVE_RVALUE_REFERENCES
#include <utility>
#endif

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup TypeTree
       *  \ingroup PDELab
       *  \{
       */

      //! Visitor interface and base class for TypeTree visitors.
      /**
       * TypeTreeVisitor defines the interface for visitors that can be applied to a TypeTree
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

        //! Method for postfix traversal.
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

        template<typename T, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(T&& t, Child&& child, TreePath treePath, ChildIndex childIndex) const {}

        template<typename T, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(T&& t, Child&& child, TreePath treePath, ChildIndex childIndex) const {}

#else // HAVE_RVALUE_REFERENCES

        // These are just a repeat of the above if the compiler does not support
        // rvalue references. In this case, we need variants for const and non-const
        // nodes.

        // Method for prefix traversal
        template<typename T, typename TreePath>
        void pre(T& t, TreePath treePath) const {}

        // Method for infix traversal
        template<typename T, typename TreePath>
        void in(T& t, TreePath treePath) const {}

        // Method for postfix traversal
        template<typename T, typename TreePath>
        void post(T& t, TreePath treePath) const {}

        // Method for leaf traversal
        template<typename T, typename TreePath>
        void leaf(T& t, TreePath treePath) const {}

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

        template<typename T, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(const T& t, const Child& child, TreePath treePath, ChildIndex childIndex) const {}

        template<typename T, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(T& t, Child& child, TreePath treePath, ChildIndex childIndex) const {}

        template<typename T, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(const T& t, const Child& child, TreePath treePath, ChildIndex childIndex) const {}

        template<typename T, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(T& t, Child& child, TreePath treePath, ChildIndex childIndex) const {}

#endif // HAVE_RVALUE_REFERENCES || DOXYGEN

      };

      struct VisitDirectChildren
      {
        template<typename Node, typename Child, typename TreePath>
        struct VisitChild
        {
          static const bool value = false;
        };
      };

      struct VisitTree
      {
        template<typename Node, typename Child, typename TreePath>
        struct VisitChild
        {
          static const bool value = true;
        };
      };

      struct TreeVisitor
        : public DefaultVisitor
        , public VisitTree
      {};

      struct DirectChildrenVisitor
        : public DefaultVisitor
        , public VisitDirectChildren
      {};


#ifndef DOXYGEN // these are all internals and not public API. Only access is using applyToTree().

      template<TreePathType::Type tpType, typename tag = StartTag, bool doApply = true>
      struct ApplyToTree
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
          ApplyToTree<typename Node::NodeTag,tpType>::apply(node,visitor,TreePathFactory<tpType>::create().mutablePath());
        }

        template<typename Node, typename Visitor>
        static void apply(const Node& node, Visitor& visitor)
        {
          ApplyToTree<typename Node::NodeTag,tpType>::apply(node,visitor,TreePathFactory<tpType>::create().mutablePath());
        }

        template<typename Node, typename Visitor>
        static void apply(Node& node, const Visitor& visitor)
        {
          ApplyToTree<typename Node::NodeTag,tpType>::apply(node,visitor,TreePathFactory<tpType>::create().mutablePath());
        }

        template<typename Node, typename Visitor>
        static void apply(const Node& node, const Visitor& visitor)
        {
          ApplyToTree<typename Node::NodeTag,tpType>::apply(node,visitor,TreePathFactory<tpType>::create().mutablePath());
        }

#endif // HAVE_RVALUE_REFERENCES

      };


      // Do not visit nodes the visitor is not interested in
      template<TreePathType::Type tpType, typename tag>
      struct ApplyToTree<tpType,tag,false>
      {
        template<typename Node, typename Visitor, typename TreePath>
        static void apply(const Node& node, const Visitor& visitor, TreePath treePath)
        {}
      };


      // LeafNode
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

      namespace {

        // For the CompositeNode, we do again need a TMP for iterating over the
        // children. Note that we use an index that counts down instead of up.
        // This allows us to specialize the TMP for the last child, where we
        // do not want to invoke the infix visitor on the CompositeNode.
        template<std::size_t inverse_k, std::size_t count>
        struct apply_to_children_fully_static
        {

#if HAVE_RVALUE_REFERENCES

          template<typename N, typename V, typename TreePath>
          static void apply(N&& n, V&& v, TreePath tp)
          {
            typedef typename remove_reference<N>::type Node;
            typedef typename remove_reference<V>::type Visitor;
            typedef typename Node::template Child<count-inverse_k>::Type C;
            typedef typename TreePathPushBack<TreePath,count-inverse_k>::type ChildTreePath;
            const bool visit = Visitor::template VisitChild<Node,C,ChildTreePath>::value;
            v.beforeChild(std::forward<N>(n),n.template child<count-inverse_k>(),tp,std::integral_constant<std::size_t,count-inverse_k>());
            ApplyToTree<Visitor::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-inverse_k>(),
                                                                                std::forward<V>(v),
                                                                                ChildTreePath());
            v.afterChild(std::forward<N>(n),n.template child<count-inverse_k>(),tp,std::integral_constant<std::size_t,count-inverse_k>());
            v.in(std::forward<N>(n),tp);
            apply_to_children_fully_static<inverse_k-1,count>::apply(std::forward<N>(n),
                                                                     std::forward<V>(v),
                                                                     tp);
          }

#else

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, V& v, TreePath tp)
          {
            typedef typename N::template Child<count-inverse_k>::Type C;
            typedef typename TreePathPushBack<TreePath,count-inverse_k>::type ChildTreePath;
            const bool visit = V::template VisitChild<N,C,ChildTreePath>::value;
            v.beforeChild(n,n.template child<count-inverse_k>(),tp,std::integral_constant<std::size_t,count-inverse_k>());
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-inverse_k>(),v,ChildTreePath());
            v.afterChild(n,n.template child<count-inverse_k>(),tp,std::integral_constant<std::size_t,count-inverse_k>());
            v.in(n,tp);
            apply_to_children_fully_static<inverse_k-1,count>::apply(n,v,tp);
          }

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, V& v, TreePath tp)
          {
            typedef typename N::template Child<count-inverse_k>::Type C;
            typedef typename TreePathPushBack<TreePath,count-inverse_k>::type ChildTreePath;
            const bool visit = V::template VisitChild<N,C,ChildTreePath>::value;
            v.beforeChild(n,n.template child<count-inverse_k>(),tp,std::integral_constant<std::size_t,count-inverse_k>());
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-inverse_k>(),v,ChildTreePath());
            v.afterChild(n,n.template child<count-inverse_k>(),tp,std::integral_constant<std::size_t,count-inverse_k>());
            v.in(n,tp);
            apply_to_children_fully_static<inverse_k-1,count>::apply(n,v,tp);
          }

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, const V& v, TreePath tp)
          {
            typedef typename N::template Child<count-inverse_k>::Type C;
            typedef typename TreePathPushBack<TreePath,count-inverse_k>::type ChildTreePath;
            const bool visit = V::template VisitChild<N,C,ChildTreePath>::value;
            v.beforeChild(n,n.template child<count-inverse_k>(),tp,std::integral_constant<std::size_t,count-inverse_k>());
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-inverse_k>(),v,ChildTreePath());
            v.afterChild(n,n.template child<count-inverse_k>(),tp,std::integral_constant<std::size_t,count-inverse_k>());
            v.in(n,tp);
            apply_to_children_fully_static<inverse_k-1,count>::apply(n,v,tp);
          }

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, const V& v, TreePath tp)
          {
            typedef typename N::template Child<count-inverse_k>::Type C;
            typedef typename TreePathPushBack<TreePath,count-inverse_k>::type ChildTreePath;
            const bool visit = V::template VisitChild<N,C,ChildTreePath>::value;
            v.beforeChild(n,n.template child<count-inverse_k>(),tp,std::integral_constant<std::size_t,count-inverse_k>());
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-inverse_k>(),v,ChildTreePath());
            v.afterChild(n,n.template child<count-inverse_k>(),tp,std::integral_constant<std::size_t,count-inverse_k>());
            v.in(n,tp);
            apply_to_children_fully_static<inverse_k-1,count>::apply(n,v,tp);
          }

#endif // HAVE_RVALUE_REFERENCES

        };

        // Specialization for last child. This specialization stops the recursion and
        // does not call the infix visitor on the CompositeNode.
        template<std::size_t count>
        struct apply_to_children_fully_static<1,count>
        {

#if HAVE_RVALUE_REFERENCES

          template<typename N, typename V, typename TreePath>
          static void apply(N&& n, V&& v, TreePath tp)
          {
            typedef typename remove_reference<N>::type Node;
            typedef typename remove_reference<V>::type Visitor;
            typedef typename Node::template Child<count-1>::Type C;
            typedef typename TreePathPushBack<TreePath,count-1>::type ChildTreePath;
            const bool visit = Visitor::template VisitChild<Node,C,ChildTreePath>::value;
            v.beforeChild(std::forward<N>(n),n.template child<count-1>(),tp,std::integral_constant<std::size_t,count-1>());
            ApplyToTree<Visitor::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-1>(),
                                                                                std::forward<V>(v),
                                                                                ChildTreePath());
            v.afterChild(std::forward<N>(n),n.template child<count-1>(),tp,std::integral_constant<std::size_t,count-1>());
          }

#else

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, V& v, TreePath tp)
          {
            typedef typename N::template Child<count-1>::Type C;
            typedef typename TreePathPushBack<TreePath,count-1>::type ChildTreePath;
            const bool visit = V::template VisitChild<N,C,ChildTreePath>::value;
            v.beforeChild(n,n.template child<count-1>(),tp,std::integral_constant<std::size_t,count-1>());
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-1>(),v,ChildTreePath());
            v.afterChild(n,n.template child<count-1>(),tp,std::integral_constant<std::size_t,count-1>());
          }

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, V& v, TreePath tp)
          {
            typedef typename N::template Child<count-1>::Type C;
            typedef typename TreePathPushBack<TreePath,count-1>::type ChildTreePath;
            const bool visit = V::template VisitChild<N,C,ChildTreePath>::value;
            v.beforeChild(n,n.template child<count-1>(),tp,std::integral_constant<std::size_t,count-1>());
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-1>(),v,ChildTreePath());
            v.afterChild(n,n.template child<count-1>(),tp,std::integral_constant<std::size_t,count-1>());
          }

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, const V& v, TreePath tp)
          {
            typedef typename N::template Child<count-1>::Type C;
            typedef typename TreePathPushBack<TreePath,count-1>::type ChildTreePath;
            const bool visit = V::template VisitChild<N,C,ChildTreePath>::value;
            v.beforeChild(n,n.template child<count-1>(),tp,std::integral_constant<std::size_t,count-1>());
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-1>(),v,ChildTreePath());
            v.afterChild(n,n.template child<count-1>(),tp,std::integral_constant<std::size_t,count-1>());
          }

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, const V& v, TreePath tp)
          {
            typedef typename N::template Child<count-1>::Type C;
            typedef typename TreePathPushBack<TreePath,count-1>::type ChildTreePath;
            const bool visit = V::template VisitChild<N,C,ChildTreePath>::value;
            v.beforeChild(n,n.template child<count-1>(),tp,std::integral_constant<std::size_t,count-1>());
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-1>(),v,ChildTreePath());
            v.afterChild(n,n.template child<count-1>(),tp,std::integral_constant<std::size_t,count-1>());
          }

#endif // HAVE_RVALUE_REFERENCES
        };

        // Specialization for CompositeNode without any children.
        template<>
        struct apply_to_children_fully_static<0,0>
        {

#if HAVE_RVALUE_REFERENCES

          template<typename N, typename V, typename TreePath>
          static void apply(N&& n, V&& v, TreePath tp) {}

#else

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, V& v, TreePath tp) {}

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, V& v, TreePath tp) {}

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, const V& v, TreePath tp) {}

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, const V& v, TreePath tp) {}

#endif // HAVE_RVALUE_REFERENCES

        };

        // For the CompositeNode, we do again need a TMP for iterating over the
        // children. Note that we use an index that counts down instead of up.
        // This allows us to specialize the TMP for the last child, where we
        // do not want to invoke the infix visitor on the CompositeNode.
        template<std::size_t inverse_k, std::size_t count>
        struct apply_to_children_dynamic
        {

#if HAVE_RVALUE_REFERENCES

          template<typename N, typename V, typename TreePath>
          static void apply(N&& n, V&& v, TreePath tp)
          {
            typedef typename remove_reference<N>::type Node;
            typedef typename remove_reference<V>::type Visitor;
            typedef typename Node::template Child<count-inverse_k>::Type C;
            const bool visit = Visitor::template VisitChild<Node,C,typename TreePath::ViewType>::value;
            v.beforeChild(std::forward<N>(n),n.template child<count-inverse_k>(),tp.view(),count-inverse_k);
            tp.push_back(count-inverse_k);
            ApplyToTree<Visitor::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-inverse_k>(),
                                                                                std::forward<V>(v),
                                                                                tp);
            tp.pop_back();
            v.afterChild(std::forward<N>(n),n.template child<count-inverse_k>(),tp.view(),count-inverse_k);
            v.in(std::forward<N>(n),tp.view());
            apply_to_children_dynamic<inverse_k-1,count>::apply(std::forward<N>(n),
                                                                std::forward<V>(v),
                                                                tp);
          }

#else

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, V& v, TreePath tp)
          {
            typedef typename N::template Child<count-inverse_k>::Type C;
            const bool visit = V::template VisitChild<N,C,typename TreePath::ViewType>::value;
            v.beforeChild(n,n.template child<count-inverse_k>(),tp.view(),count-inverse_k);
            tp.push_back(count-inverse_k);
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-inverse_k>(),
                                                                          v,
                                                                          tp);
            tp.pop_back();
            v.afterChild(n,n.template child<count-inverse_k>(),tp.view(),count-inverse_k);
            v.in(n,tp.view());
            apply_to_children_dynamic<inverse_k-1,count>::apply(n,v,tp);
          }

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, V& v, TreePath tp)
          {
            typedef typename N::template Child<count-inverse_k>::Type C;
            const bool visit = V::template VisitChild<N,C,typename TreePath::ViewType>::value;
            v.beforeChild(n,n.template child<count-inverse_k>(),tp.view(),count-inverse_k);
            tp.push_back(count-inverse_k);
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-inverse_k>(),
                                                                          v,
                                                                          tp);
            tp.pop_back();
            v.afterChild(n,n.template child<count-inverse_k>(),tp.view(),count-inverse_k);
            v.in(n,tp.view());
            apply_to_children_dynamic<inverse_k-1,count>::apply(n,v,tp);
          }

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, const V& v, TreePath tp)
          {
            typedef typename N::template Child<count-inverse_k>::Type C;
            const bool visit = V::template VisitChild<N,C,typename TreePath::ViewType>::value;
            v.beforeChild(n,n.template child<count-inverse_k>(),tp.view(),count-inverse_k);
            tp.push_back(count-inverse_k);
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-inverse_k>(),
                                                                          v,
                                                                          tp);
            tp.pop_back();
            v.afterChild(n,n.template child<count-inverse_k>(),tp.view(),count-inverse_k);
            v.in(n,tp.view());
            apply_to_children_dynamic<inverse_k-1,count>::apply(n,v,tp);
          }

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, const V& v, TreePath tp)
          {
            typedef typename N::template Child<count-inverse_k>::Type C;
            const bool visit = V::template VisitChild<N,C,typename TreePath::ViewType>::value;
            v.beforeChild(n,n.template child<count-inverse_k>(),tp.view(),count-inverse_k);
            tp.push_back(count-inverse_k);
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-inverse_k>(),
                                                                          v,
                                                                          tp);
            tp.pop_back();
            v.afterChild(n,n.template child<count-inverse_k>(),tp.view(),count-inverse_k);
            v.in(n,tp.view());
            apply_to_children_dynamic<inverse_k-1,count>::apply(n,v,tp);
          }

#endif // HAVE_RVALUE_REFERENCES

        };

        // Specialization for last child. This specialization stops the recursion and
        // does not call the infix visitor on the CompositeNode.
        template<std::size_t count>
        struct apply_to_children_dynamic<1,count>
        {

#if HAVE_RVALUE_REFERENCES

          template<typename N, typename V, typename TreePath>
          static void apply(N&& n, V&& v, TreePath tp)
          {
            typedef typename remove_reference<N>::type Node;
            typedef typename remove_reference<V>::type Visitor;
            typedef typename Node::template Child<count-1>::Type C;
            const bool visit = Visitor::template VisitChild<Node,C,typename TreePath::ViewType>::value;
            v.beforeChild(std::forward<N>(n),n.template child<count-1>(),tp.view(),count-1);
            tp.push_back(count-1);
            ApplyToTree<Visitor::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-1>(),
                                                                                std::forward<V>(v),
                                                                                tp);
            tp.pop_back();
            v.afterChild(std::forward<N>(n),n.template child<count-1>(),tp.view(),count-1);
          }

#else

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, V& v, TreePath tp)
          {
            typedef typename N::template Child<count-1>::Type C;
            const bool visit = V::template VisitChild<N,C,typename TreePath::ViewType>::value;
            v.beforeChild(n,n.template child<count-1>(),tp.view(),count-1);
            tp.push_back(count-1);
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-1>(),
                                                                          v,
                                                                          tp);
            tp.pop_back();
            v.afterChild(n,n.template child<count-1>(),tp.view(),count-1);
          }

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, V& v, TreePath tp)
          {
            typedef typename N::template Child<count-1>::Type C;
            const bool visit = V::template VisitChild<N,C,typename TreePath::ViewType>::value;
            v.beforeChild(n,n.template child<count-1>(),tp.view(),count-1);
            tp.push_back(count-1);
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-1>(),
                                                                          v,
                                                                          tp);
            tp.pop_back();
            v.afterChild(n,n.template child<count-1>(),tp.view(),count-1);
          }

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, const V& v, TreePath tp)
          {
            typedef typename N::template Child<count-1>::Type C;
            const bool visit = V::template VisitChild<N,C,typename TreePath::ViewType>::value;
            v.beforeChild(n,n.template child<count-1>(),tp.view(),count-1);
            tp.push_back(count-1);
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-1>(),
                                                                          v,
                                                                          tp);
            tp.pop_back();
            v.afterChild(n,n.template child<count-1>(),tp.view(),count-1);
          }

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, const V& v, TreePath tp)
          {
            typedef typename N::template Child<count-1>::Type C;
            const bool visit = V::template VisitChild<N,C,typename TreePath::ViewType>::value;
            v.beforeChild(n,n.template child<count-1>(),tp.view(),count-1);
            tp.push_back(count-1);
            ApplyToTree<V::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-1>(),
                                                                          v,
                                                                          tp);
            tp.pop_back();
            v.afterChild(n,n.template child<count-1>(),tp.view(),count-1);
          }

#endif // HAVE_RVALUE_REFERENCES
        };

        // Specialization for CompositeNode without any children.
        template<>
        struct apply_to_children_dynamic<0,0>
        {

#if HAVE_RVALUE_REFERENCES

          template<typename N, typename V, typename TreePath>
          static void apply(N&& n, V&& v, TreePath tp) {}

#else

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, V& v, TreePath tp) {}

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, V& v, TreePath tp) {}

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, const V& v, TreePath tp) {}

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, const V& v, TreePath tp) {}

#endif // HAVE_RVALUE_REFERENCES

        };

        template<TreePathType::Type treePathType, std::size_t CHILDREN>
        struct apply_to_children;

        template<std::size_t CHILDREN>
        struct apply_to_children<TreePathType::fullyStatic,CHILDREN>
          : public apply_to_children_fully_static<CHILDREN,CHILDREN>
        {};

        template<std::size_t CHILDREN>
        struct apply_to_children<TreePathType::dynamic,CHILDREN>
          : public apply_to_children_dynamic<CHILDREN,CHILDREN>
        {};


      } // anonymous namespace

        // Base class for composite node traversal
      struct ApplyToGenericCompositeNode
      {

#if HAVE_RVALUE_REFERENCES

        template<typename N, typename V, typename TreePath>
        static void apply(N&& n, V&& v, TreePath tp)
        {
          v.pre(std::forward<N>(n),tp);
          typedef typename remove_reference<N>::type Node;
          typedef typename remove_reference<V>::type Visitor;
          apply_to_children<Visitor::treePathType,Node::CHILDREN>::apply(std::forward<N>(n),
                                                                         std::forward<V>(v),
                                                                         tp);
          v.post(std::forward<N>(n),tp);
        }

#else

        template<typename N, typename V, typename TreePath>
        static void apply(N& n, V& v, TreePath tp)
        {
          v.pre(n,tp);
          apply_to_children<V::treePathType,N::CHILDREN>::apply(n,v,tp);
          v.post(n,tp);
        }

        template<typename N, typename V, typename TreePath>
        static void apply(const N& n, V& v, TreePath tp)
        {
          v.pre(n,tp);
          apply_to_children<V::treePathType,N::CHILDREN>::apply(n,v,tp);
          v.post(n,tp);
        }

        template<typename N, typename V, typename TreePath>
        static void apply(N& n, const V& v, TreePath tp)
        {
          v.pre(n,tp);
          apply_to_children<V::treePathType,N::CHILDREN>::apply(n,v,tp);
          v.post(n,tp);
        }

        template<typename N, typename V, typename TreePath>
        static void apply(const N& n, const V& v, TreePath tp)
        {
          v.pre(n,tp);
          apply_to_children<V::treePathType,N::CHILDREN>::apply(n,v,tp);
          v.post(n,tp);
        }

#endif // HAVE_RVALUE_REFERENCES

      };

      // Traverse PowerNode
      template<>
      struct ApplyToTree<TreePathType::fullyStatic,PowerNodeTag,true>
        : public ApplyToGenericCompositeNode
      {
      };


      template<>
      struct ApplyToTree<TreePathType::dynamic,PowerNodeTag,true>
      {

#if HAVE_RVALUE_REFERENCES

        template<typename N, typename V, typename TreePath>
        static void apply(N&& n, V&& v, TreePath tp)
        {
          v.pre(std::forward<N>(n),tp.view());
          typedef typename remove_reference<N>::type Node;
          typedef typename remove_reference<V>::type Visitor;
          typedef typename Node::template Child<0>::Type C;
          const bool visit = Visitor::template VisitChild<Node,C,typename TreePath::ViewType>::value;
          for (std::size_t k = 0; k < Node::CHILDREN; ++k)
            {
              v.beforeChild(std::forward<N>(n),n.child(k),tp.view(),k);
              tp.push_back(k);
              ApplyToTree<Visitor::treePathType,typename C::NodeTag,visit>::apply(n.child(k),std::forward<V>(v),tp);
              tp.pop_back();
              v.afterChild(std::forward<N>(n),n.child(k),tp.view(),k);
              if (k < Node::CHILDREN-1)
                v.in(std::forward<N>(n),tp.view());
            }
          v.post(std::forward<N>(n),tp.view());
        }

#else

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


      // Traverse CompositeNode
      template<TreePathType::Type treePathType>
      struct ApplyToTree<treePathType,CompositeNodeTag,true>
        : public ApplyToGenericCompositeNode
      {
      };


      // Traverse VariadicCompositeNode
      template<TreePathType::Type treePathType>
      struct ApplyToTree<treePathType,VariadicCompositeNodeTag,true>
        : public ApplyToGenericCompositeNode
      {
      };

#endif // DOXYGEN

#if HAVE_RVALUE_REFERENCES

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

#endif // HAVE_RVALUE_REFERENCES

      //! \} group TypeTree

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_TRAVERSAL_HH
