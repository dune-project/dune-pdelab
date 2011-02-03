// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_PAIRTRAVERSAL_HH
#define DUNE_PDELAB_COMMON_TYPETREE_PAIRTRAVERSAL_HH

#include <dune/pdelab/common/typetree/nodetags.hh>
#include <dune/pdelab/common/typetree/treepath.hh>

#if HAVE_RVALUE_REFERENCES
#include <utility>
#endif

namespace Dune {
  namespace PDELab {
    namespace TypeTree {
      namespace PairTraversal {

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

        struct VisitDirectChildren
        {
          template<typename Node1, typename Child1, typename Node2, typename Child2, typename TreePath>
          struct VisitChild
          {
            static const bool value = false;
          };
        };

        struct VisitTree
        {
          template<typename Node1, typename Child1, typename Node2, typename Child2, typename TreePath>
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


        template<TreePathType::Type tpType, typename tag1 = StartTag, typename tag2 = StartTag, bool doApply = true>
        struct ApplyToTree;

        template<TreePathType::Type tpType>
        struct ApplyToTree<tpType,StartTag,StartTag,true>
        {

#if HAVE_RVALUE_REFERENCES

          template<typename Node1, typename Node2, typename Visitor>
          static void apply(Node1&& node1, Node2&& node2, Visitor&& visitor)
          {
            ApplyToTree<tpType,
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
            ApplyToTree<tpType,
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
            ApplyToTree<tpType,
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
            ApplyToTree<tpType,
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
            ApplyToTree<tpType,
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
        struct ApplyToTree<tpType,Tag1,Tag2,false>
        {
          template<typename Node1, typename Node2, typename Visitor, typename TreePath>
          static void apply(const Node1& node1, const Node2& node2, const Visitor& visitor, TreePath treePath)
          {}
        };


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

        namespace {

          // TMP for traversing the children of non-leaf nodes with a static TreePath.
          // Due to the static TreePath, we have to use this TMP for both CompositeNode
          // and PowerNode.
          template<std::size_t inverse_k, std::size_t count>
          struct apply_to_children_fully_static
          {

#if HAVE_RVALUE_REFERENCES

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1&& n1, N2&& n2, V&& v, TreePath tp)
            {
              typedef typename remove_reference<N1>::type::template Child<count-inverse_k>::Type C1;
              typedef typename remove_reference<N2>::type::template Child<count-inverse_k>::Type C2;
              typedef typename TreePathPushBack<TreePath,count-inverse_k>::type ChildTreePath;
              const bool visit = std::remove_reference<V>::type
                ::template VisitChild<typename remove_reference<N1>::type,
                                      C1,
                                      typename remove_reference<N2>::type,
                                      C2,
                                      ChildTreePath>::value;
              v.beforeChild(std::forward<N1>(n1),n1.template child<count-inverse_k>(),
                            std::forward<N2>(n2),n2.template child<count-inverse_k>(),
                            tp,std::integral_constant<std::size_t,count-inverse_k>());
              ApplyToTree<std::remove_reference<V>::type::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-inverse_k>(),
                                        n2.template child<count-inverse_k>(),
                                        std::forward<V>(v),
                                        ChildTreePath());
              v.afterChild(std::forward<N1>(n1),n1.template child<count-inverse_k>(),
                           std::forward<N2>(n2),n2.template child<count-inverse_k>(),
                           tp,std::integral_constant<std::size_t,count-inverse_k>());
              v.in(std::forward<N1>(n1),std::forward<N2>(n2),tp);
              apply_to_children_fully_static<inverse_k-1,count>::apply(std::forward<N1>(n1),
                                                                       std::forward<N2>(n2),
                                                                       std::forward<V>(v),
                                                                       tp);
            }

#else

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1& n1, N2& n2, V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-inverse_k>::Type C1;
              typedef typename N2::template Child<count-inverse_k>::Type C2;
              typedef typename TreePathPushBack<TreePath,count-inverse_k>::type ChildTreePath;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,ChildTreePath>::value;
              v.beforeChild(n1,n1.template child<count-inverse_k>(),
                            n2,n2.template child<count-inverse_k>(),
                            tp,
                            std::integral_constant<std::size_t,count-inverse_k>());
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-inverse_k>(),
                                        n2.template child<count-inverse_k>(),
                                        v,
                                        ChildTreePath());
              v.afterChild(n1,n1.template child<count-inverse_k>(),
                           n2,n2.template child<count-inverse_k>(),
                           tp,
                           std::integral_constant<std::size_t,count-inverse_k>());
              v.in(n1,n2,tp);
              apply_to_children_fully_static<inverse_k-1,count>::apply(n1,n2,v,tp);
            }

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(const N1& n1, const N2& n2, V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-inverse_k>::Type C1;
              typedef typename N2::template Child<count-inverse_k>::Type C2;
              typedef typename TreePathPushBack<TreePath,count-inverse_k>::type ChildTreePath;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,ChildTreePath>::value;
              v.beforeChild(n1,n1.template child<count-inverse_k>(),
                            n2,n2.template child<count-inverse_k>(),
                            tp,
                            std::integral_constant<std::size_t,count-inverse_k>());
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-inverse_k>(),
                                        n2.template child<count-inverse_k>(),
                                        v,
                                        ChildTreePath());
              v.afterChild(n1,n1.template child<count-inverse_k>(),
                           n2,n2.template child<count-inverse_k>(),
                           tp,
                           std::integral_constant<std::size_t,count-inverse_k>());
              v.in(n1,n2,tp);
              apply_to_children_fully_static<inverse_k-1,count>::apply(n1,n2,v,tp);
            }

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1& n1, N2& n2, const V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-inverse_k>::Type C1;
              typedef typename N2::template Child<count-inverse_k>::Type C2;
              typedef typename TreePathPushBack<TreePath,count-inverse_k>::type ChildTreePath;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,ChildTreePath>::value;
              v.beforeChild(n1,n1.template child<count-inverse_k>(),
                            n2,n2.template child<count-inverse_k>(),
                            tp,
                            std::integral_constant<std::size_t,count-inverse_k>());
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-inverse_k>(),
                                        n2.template child<count-inverse_k>(),
                                        v,
                                        ChildTreePath());
              v.afterChild(n1,n1.template child<count-inverse_k>(),
                           n2,n2.template child<count-inverse_k>(),
                           tp,
                           std::integral_constant<std::size_t,count-inverse_k>());
              v.in(n1,n2,tp);
              apply_to_children_fully_static<inverse_k-1,count>::apply(n1,n2,v,tp);
            }

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(const N1& n1, const N2& n2, const V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-inverse_k>::Type C1;
              typedef typename N2::template Child<count-inverse_k>::Type C2;
              typedef typename TreePathPushBack<TreePath,count-inverse_k>::type ChildTreePath;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,ChildTreePath>::value;
              v.beforeChild(n1,n1.template child<count-inverse_k>(),
                            n2,n2.template child<count-inverse_k>(),
                            tp,
                            std::integral_constant<std::size_t,count-inverse_k>());
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-inverse_k>(),
                                        n2.template child<count-inverse_k>(),
                                        v,
                                        ChildTreePath());
              v.afterChild(n1,n1.template child<count-inverse_k>(),
                           n2,n2.template child<count-inverse_k>(),
                           tp,
                           std::integral_constant<std::size_t,count-inverse_k>());
              v.in(n1,n2,tp);
              apply_to_children_fully_static<inverse_k-1,count>::apply(n1,n2,v,tp);
            }

#endif // HAVE_RVALUE_REFERENCES

          };

          // Specialization for last child. This specialization stops the recursion and
          // does not call the infix visitor.
          template<std::size_t count>
          struct apply_to_children_fully_static<1,count>
          {

#if HAVE_RVALUE_REFERENCES

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1&& n1, N2&& n2, V&& v, TreePath tp)
            {
              typedef typename remove_reference<N1>::type::template Child<count-1>::Type C1;
              typedef typename remove_reference<N2>::type::template Child<count-1>::Type C2;
              typedef typename TreePathPushBack<TreePath,count-1>::type ChildTreePath;
              const bool visit = std::remove_reference<V>::type
                ::template VisitChild<typename remove_reference<N1>::type,
                                      C1,
                                      typename remove_reference<N2>::type,
                                      C2,
                                      ChildTreePath>::value;
              v.beforeChild(std::forward<N1>(n1),n1.template child<count-1>(),
                            std::forward<N2>(n2),n2.template child<count-1>(),
                            tp,std::integral_constant<std::size_t,count-1>());
              ApplyToTree<std::remove_reference<V>::type::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-1>(),
                                        n2.template child<count-1>(),
                                        std::forward<V>(v),
                                        ChildTreePath());
              v.afterChild(std::forward<N1>(n1),n1.template child<count-1>(),
                           std::forward<N2>(n2),n2.template child<count-1>(),
                           tp,std::integral_constant<std::size_t,count-1>());
            }

#else

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1& n1, N2& n2, V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-1>::Type C1;
              typedef typename N2::template Child<count-1>::Type C2;
              typedef typename TreePathPushBack<TreePath,count-1>::type ChildTreePath;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,ChildTreePath>::value;
              v.beforeChild(n1,n1.template child<count-1>(),
                            n2,n2.template child<count-1>(),
                            tp,
                            std::integral_constant<std::size_t,count-1>());
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-1>(),
                                        n2.template child<count-1>(),
                                        v,
                                        ChildTreePath());
              v.afterChild(n1,n1.template child<count-1>(),
                           n2,n2.template child<count-1>(),
                           tp,
                           std::integral_constant<std::size_t,count-1>());
            }

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(const N1& n1, const N2& n2, V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-1>::Type C1;
              typedef typename N2::template Child<count-1>::Type C2;
              typedef typename TreePathPushBack<TreePath,count-1>::type ChildTreePath;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,ChildTreePath>::value;
              v.beforeChild(n1,n1.template child<count-1>(),
                            n2,n2.template child<count-1>(),
                            tp,
                            std::integral_constant<std::size_t,count-1>());
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-1>(),
                                        n2.template child<count-1>(),
                                        v,
                                        ChildTreePath());
              v.afterChild(n1,n1.template child<count-1>(),
                           n2,n2.template child<count-1>(),
                           tp,
                           std::integral_constant<std::size_t,count-1>());
            }

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1& n1, N2& n2, const V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-1>::Type C1;
              typedef typename N2::template Child<count-1>::Type C2;
              typedef typename TreePathPushBack<TreePath,count-1>::type ChildTreePath;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,ChildTreePath>::value;
              v.beforeChild(n1,n1.template child<count-1>(),
                            n2,n2.template child<count-1>(),
                            tp,
                            std::integral_constant<std::size_t,count-1>());
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-1>(),
                                        n2.template child<count-1>(),
                                        v,
                                        ChildTreePath());
              v.afterChild(n1,n1.template child<count-1>(),
                           n2,n2.template child<count-1>(),
                           tp,
                           std::integral_constant<std::size_t,count-1>());
            }

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(const N1& n1, const N2& n2, const V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-1>::Type C1;
              typedef typename N2::template Child<count-1>::Type C2;
              typedef typename TreePathPushBack<TreePath,count-1>::type ChildTreePath;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,ChildTreePath>::value;
              v.beforeChild(n1,n1.template child<count-1>(),
                            n2,n2.template child<count-1>(),
                            tp,
                            std::integral_constant<std::size_t,count-1>());
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-1>(),
                                        n2.template child<count-1>(),
                                        v,
                                        ChildTreePath());
              v.afterChild(n1,n1.template child<count-1>(),
                           n2,n2.template child<count-1>(),
                           tp,
                           std::integral_constant<std::size_t,count-1>());
            }

#endif // HAVE_RVALUE_REFERENCES

          };

          // Specialization for CompositeNode without any children.
          template<>
          struct apply_to_children_fully_static<0,0>
          {

#if HAVE_RVALUE_REFERENCES

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1&& n1, N2&& n2, V&& v, TreePath tp) {}

#else

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1& n1, N2& n2, V& v, TreePath tp) {}

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(const N1& n1, const N2& n2, V& v, TreePath tp) {}

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1& n1, N2& n2, const V& v, TreePath tp) {}

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(const N1& n1, const N2& n2, const V& v, TreePath tp) {}

#endif // HAVE_RVALUE_REFERENCES

          };


          // TMP for traversing the children of non-leaf nodes with a dynamic TreePath.
          // This is used if at least one of the nodes is a CompositeNode.
          // If both nodes are PowerNodes, this will not be instantiated and we simply
          // use a runtime loop to iterate over the children.
          template<std::size_t inverse_k, std::size_t count>
          struct apply_to_children_dynamic
          {

#if HAVE_RVALUE_REFERENCES

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1&& n1, N2&& n2, V&& v, TreePath tp)
            {
              typedef typename remove_reference<N1>::type::template Child<count-inverse_k>::Type C1;
              typedef typename remove_reference<N2>::type::template Child<count-inverse_k>::Type C2;
              const bool visit = std::remove_reference<V>::type
                ::template VisitChild<typename remove_reference<N1>::type,
                                      C1,
                                      typename remove_reference<N2>::type,
                                      C2,
                                      typename TreePath::ViewType>::value;
              v.beforeChild(std::forward<N1>(n1),n1.template child<count-inverse_k>(),
                            std::forward<N2>(n2),n2.template child<count-inverse_k>(),
                            tp.view(),count-inverse_k);
              tp.push_back(count-inverse_k);
              ApplyToTree<std::remove_reference<V>::type::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-inverse_k>(),
                                        n2.template child<count-inverse_k>(),
                                        std::forward<V>(v),
                                        tp);
              tp.pop_back();
              v.afterChild(std::forward<N1>(n1),n1.template child<count-inverse_k>(),
                           std::forward<N2>(n2),n2.template child<count-inverse_k>(),
                           tp.view(),count-inverse_k);
              v.in(std::forward<N1>(n1),std::forward<N2>(n2),tp.view());
              apply_to_children_dynamic<inverse_k-1,count>::apply(std::forward<N1>(n1),
                                                                  std::forward<N2>(n2),
                                                                  std::forward<V>(v),
                                                                  tp);
            }

#else

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1& n1, N2& n2, V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-inverse_k>::Type C1;
              typedef typename N2::template Child<count-inverse_k>::Type C2;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,typename TreePath::ViewType>::value;
              v.beforeChild(n1,n1.template child<count-inverse_k>(),
                            n2,n2.template child<count-inverse_k>(),
                            tp.view(),count-inverse_k);
              tp.push_back(count-inverse_k);
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-inverse_k>(),
                                        n2.template child<count-inverse_k>(),
                                        v,
                                        tp);
              tp.pop_back();
              v.afterChild(n1,n1.template child<count-inverse_k>(),
                           n2,n2.template child<count-inverse_k>(),
                           tp.view(),count-inverse_k);
              v.in(n1,n2,tp.view());
              apply_to_children_dynamic<inverse_k-1,count>::apply(n1,n2,v,tp);
            }

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(const N1& n1, const N2& n2, V& v, TreePath tp)            {
              typedef typename N1::template Child<count-inverse_k>::Type C1;
              typedef typename N2::template Child<count-inverse_k>::Type C2;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,typename TreePath::ViewType>::value;
              v.beforeChild(n1,n1.template child<count-inverse_k>(),
                            n2,n2.template child<count-inverse_k>(),
                            tp.view(),count-inverse_k);
              tp.push_back(count-inverse_k);
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-inverse_k>(),
                                        n2.template child<count-inverse_k>(),
                                        v,
                                        tp);
              tp.pop_back();
              v.afterChild(n1,n1.template child<count-inverse_k>(),
                           n2,n2.template child<count-inverse_k>(),
                           tp.view(),count-inverse_k);
              v.in(n1,n2,tp.view());
              apply_to_children_dynamic<inverse_k-1,count>::apply(n1,n2,v,tp);
            }

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1& n1, N2& n2, const V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-inverse_k>::Type C1;
              typedef typename N2::template Child<count-inverse_k>::Type C2;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,typename TreePath::ViewType>::value;
              v.beforeChild(n1,n1.template child<count-inverse_k>(),
                            n2,n2.template child<count-inverse_k>(),
                            tp.view(),count-inverse_k);
              tp.push_back(count-inverse_k);
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-inverse_k>(),
                                        n2.template child<count-inverse_k>(),
                                        v,
                                        tp);
              tp.pop_back();
              v.afterChild(n1,n1.template child<count-inverse_k>(),
                           n2,n2.template child<count-inverse_k>(),
                           tp.view(),count-inverse_k);
              v.in(n1,n2,tp.view());
              apply_to_children_dynamic<inverse_k-1,count>::apply(n1,n2,v,tp);
            }

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(const N1& n1, const N2& n2, const V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-inverse_k>::Type C1;
              typedef typename N2::template Child<count-inverse_k>::Type C2;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,typename TreePath::ViewType>::value;
              v.beforeChild(n1,n1.template child<count-inverse_k>(),
                            n2,n2.template child<count-inverse_k>(),
                            tp.view(),count-inverse_k);
              tp.push_back(count-inverse_k);
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-inverse_k>(),
                                        n2.template child<count-inverse_k>(),
                                        v,
                                        tp);
              tp.pop_back();
              v.afterChild(n1,n1.template child<count-inverse_k>(),
                           n2,n2.template child<count-inverse_k>(),
                           tp.view(),count-inverse_k);
              v.in(n1,n2,tp.view());
              apply_to_children_dynamic<inverse_k-1,count>::apply(n1,n2,v,tp);
            }

#endif // HAVE_RVALUE_REFERENCES

          };

          // Specialization for last child. This specialization stops the recursion and
          // does not call the infix visitor on the CompositeNode.
          template<std::size_t count>
          struct apply_to_children_dynamic<1,count>
          {

#if HAVE_RVALUE_REFERENCES

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1&& n1, N2&& n2, V&& v, TreePath tp)
            {
              typedef typename remove_reference<N1>::type::template Child<count-1>::Type C1;
              typedef typename remove_reference<N2>::type::template Child<count-1>::Type C2;
              const bool visit = std::remove_reference<V>::type
                ::template VisitChild<typename remove_reference<N1>::type,
                                      C1,
                                      typename remove_reference<N2>::type,
                                      C2,
                                      typename TreePath::ViewType>::value;
              v.beforeChild(std::forward<N1>(n1),n1.template child<count-1>(),
                            std::forward<N2>(n2),n2.template child<count-1>(),
                            tp.view(),count-1);
              tp.push_back(count-1);
              ApplyToTree<std::remove_reference<V>::type::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-1>(),
                                        n2.template child<count-1>(),
                                        std::forward<V>(v),
                                        tp);
              tp.pop_back();
              v.afterChild(std::forward<N1>(n1),n1.template child<count-1>(),
                           std::forward<N2>(n2),n2.template child<count-1>(),
                           tp.view(),count-1);
            }

#else

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1& n1, N2& n2, V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-1>::Type C1;
              typedef typename N2::template Child<count-1>::Type C2;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,typename TreePath::ViewType>::value;
              v.beforeChild(n1,n1.template child<count-1>(),
                            n2,n2.template child<count-1>(),
                            tp.view(),count-1);
              tp.push_back(count-1);
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-1>(),
                                        n2.template child<count-1>(),
                                        v,
                                        tp);
              tp.pop_back();
              v.afterChild(n1,n1.template child<count-1>(),
                           n2,n2.template child<count-1>(),
                           tp.view(),count-1);
            }

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(const N1& n1, const N2& n2, V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-1>::Type C1;
              typedef typename N2::template Child<count-1>::Type C2;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,typename TreePath::ViewType>::value;
              v.beforeChild(n1,n1.template child<count-1>(),
                            n2,n2.template child<count-1>(),
                            tp.view(),count-1);
              tp.push_back(count-1);
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-1>(),
                                        n2.template child<count-1>(),
                                        v,
                                        tp);
              tp.pop_back();
              v.afterChild(n1,n1.template child<count-1>(),
                           n2,n2.template child<count-1>(),
                           tp.view(),count-1);
            }

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1& n1, N2& n2, const V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-1>::Type C1;
              typedef typename N2::template Child<count-1>::Type C2;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,typename TreePath::ViewType>::value;
              v.beforeChild(n1,n1.template child<count-1>(),
                            n2,n2.template child<count-1>(),
                            tp.view(),count-1);
              tp.push_back(count-1);
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-1>(),
                                        n2.template child<count-1>(),
                                        v,
                                        tp);
              tp.pop_back();
              v.afterChild(n1,n1.template child<count-1>(),
                           n2,n2.template child<count-1>(),
                           tp.view(),count-1);
            }

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(const N1& n1, const N2& n2, const V& v, TreePath tp)
            {
              typedef typename N1::template Child<count-1>::Type C1;
              typedef typename N2::template Child<count-1>::Type C2;
              const bool visit = V::template VisitChild<N1,C1,N2,C2,typename TreePath::ViewType>::value;
              v.beforeChild(n1,n1.template child<count-1>(),
                            n2,n2.template child<count-1>(),
                            tp.view(),count-1);
              tp.push_back(count-1);
              ApplyToTree<V::treePathType,
                          typename C1::NodeTag,
                          typename C2::NodeTag,
                          visit>::apply(n1.template child<count-1>(),
                                        n2.template child<count-1>(),
                                        v,
                                        tp);
              tp.pop_back();
              v.afterChild(n1,n1.template child<count-1>(),
                           n2,n2.template child<count-1>(),
                           tp.view(),count-1);
            }

#endif // HAVE_RVALUE_REFERENCES

          };

          // Specialization for CompositeNode without any children.
          template<>
          struct apply_to_children_dynamic<0,0>
          {

#if HAVE_RVALUE_REFERENCES

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1&& n1, N2&& n2, V&& v, TreePath tp) {}

#else

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1& n1, N2& n2, V& v, TreePath tp) {}

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(const N1& n1, const N2& n2, V& v, TreePath tp) {}

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(N1& n1, N2& n2, const V& v, TreePath tp) {}

            template<typename N1, typename N2, typename V, typename TreePath>
            static void apply(const N1& n1, const N2& n2, const V& v, TreePath tp) {}

#endif // HAVE_RVALUE_REFERENCES

          };

        } // anonymous namespace

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


        // Base class for composite node traversal
        template<TreePathType::Type treePathType>
        struct ApplyToGenericCompositeNode
        {

#if HAVE_RVALUE_REFERENCES

          template<typename N1, typename N2, typename V, typename TreePath>
          static typename enable_if<(remove_reference<N1>::type::isLeaf || remove_reference<N2>::type::isLeaf)>::type
          apply(N1&& n1, N2&& n2, V&& v, TreePath tp)
          {
            v.leaf(std::forward<N1>(n1),std::forward<N2>(n2),tp.view());
          }


          template<typename N1, typename N2, typename V, typename TreePath>
          static typename enable_if<!(remove_reference<N1>::type::isLeaf || remove_reference<N2>::type::isLeaf)>::type
          apply(N1&& n1, N2&& n2, V&& v, TreePath tp)
          {
            v.pre(std::forward<N1>(n1),std::forward<N2>(n2),tp.view());
            typedef typename remove_reference<N1>::type Node1;
            typedef typename remove_reference<N2>::type Node2;
            dune_static_assert(Node1::CHILDREN == Node2::CHILDREN,
                               "non-leaf nodes with different numbers of children " \
                               "are not allowed during simultaneous grid traversal");
            apply_to_children<treePathType,Node1::CHILDREN>::apply(std::forward<N1>(n1),
                                                                   std::forward<N2>(n2),
                                                                   std::forward<V>(v),
                                                                   tp);
            v.post(std::forward<N1>(n1),std::forward<N2>(n2),tp.view());
          }

#else

          template<typename N1, typename N2, typename V, typename TreePath>
          static typename enable_if<N1::isLeaf || N2::isLeaf>::type
          apply(N1& n1, N2& n2, V& v, TreePath tp)
          {
            v.leaf(n1,n2,tp.view());
          }

          template<typename N1, typename N2, typename V, typename TreePath>
          static typename enable_if<N1::isLeaf || N2::isLeaf>::type
          apply(const N1& n1, const N2& n2, V& v, TreePath tp)
          {
            v.leaf(n1,n2,tp.view());
          }

          template<typename N1, typename N2, typename V, typename TreePath>
          static typename enable_if<N1::isLeaf || N2::isLeaf>::type
          apply(N1& n1, N2& n2, const V& v, TreePath tp)
          {
            v.leaf(n1,n2,tp.view());
          }

          template<typename N1, typename N2, typename V, typename TreePath>
          static typename enable_if<N1::isLeaf || N2::isLeaf>::type
          apply(const N1& n1, const N2& n2, const V& v, TreePath tp)
          {
            v.leaf(n1,n2,tp.view());
          }


          template<typename N1, typename N2, typename V, typename TreePath>
          static typename enable_if<!(N1::isLeaf || N2::isLeaf)>::type
          apply(N1& n1, N2& n2, V& v, TreePath tp)
          {
            v.pre(n1,n2,tp.view());
            dune_static_assert(N1::CHILDREN == N2::CHILDREN,
                               "non-leaf nodes with different numbers of children " \
                               "are not allowed during simultaneous grid traversal");
            apply_to_children<treePathType,N1::CHILDREN>::apply(n1,n2,v,tp);
            v.post(n1,n2,tp.view());
          }

          template<typename N1, typename N2, typename V, typename TreePath>
          static typename enable_if<!(N1::isLeaf || N2::isLeaf)>::type
          apply(const N1& n1, const N2& n2, V& v, TreePath tp)
          {
            v.pre(n1,n2,tp.view());
            dune_static_assert(N1::CHILDREN == N2::CHILDREN,
                               "non-leaf nodes with different numbers of children " \
                               "are not allowed during simultaneous grid traversal");
            apply_to_children<treePathType,N1::CHILDREN>::apply(n1,n2,v,tp);
            v.post(n1,n2,tp.view());
          }

          template<typename N1, typename N2, typename V, typename TreePath>
          static typename enable_if<!(N1::isLeaf || N2::isLeaf)>::type
          apply(N1& n1, N2& n2, const V& v, TreePath tp)
          {
            v.pre(n1,n2,tp.view());
            dune_static_assert(N1::CHILDREN == N2::CHILDREN,
                               "non-leaf nodes with different numbers of children " \
                               "are not allowed during simultaneous grid traversal");
            apply_to_children<treePathType,N1::CHILDREN>::apply(n1,n2,v,tp);
            v.post(n1,n2,tp.view());
          }

          template<typename N1, typename N2, typename V, typename TreePath>
          static typename enable_if<!(N1::isLeaf || N2::isLeaf)>::type
          apply(const N1& n1, const N2& n2, const V& v, TreePath tp)
          {
            v.pre(n1,n2,tp.view());
            dune_static_assert(N1::CHILDREN == N2::CHILDREN,
                               "non-leaf nodes with different numbers of children " \
                               "are not allowed during simultaneous grid traversal");
            apply_to_children<treePathType,N1::CHILDREN>::apply(n1,n2,v,tp);
            v.post(n1,n2,tp.view());
          }

#endif // HAVE_RVALUE_REFERENCES

        };

        // Automatically pick the correct traversal algorithm
        template<TreePathType::Type treePathType,typename FirstTag, typename SecondTag>
        struct ApplyToTree<treePathType,FirstTag,SecondTag,true>
          : public ApplyToGenericCompositeNode<treePathType>
        {
        };

        // Specialization for dynamic traversal and two PowerNodes -> use runtime iteration
        template<>
        struct ApplyToTree<TreePathType::dynamic,PowerNodeTag,PowerNodeTag,true>
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
                ApplyToTree<TreePathType::dynamic, // we know that due to the specialization
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
                ApplyToTree<TreePathType::dynamic, // we know that due to the specialization
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
                ApplyToTree<TreePathType::dynamic, // we know that due to the specialization
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
                ApplyToTree<TreePathType::dynamic, // we know that due to the specialization
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
                ApplyToTree<TreePathType::dynamic, // we know that due to the specialization
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

#if HAVE_RVALUE_REFERENCES

        template<typename Tree1, typename Tree2, typename Visitor>
        void applyToTree(Tree1&& tree1, Tree2&& tree2, Visitor&& visitor)
        {
          ApplyToTree<std::remove_reference<Visitor>::type::treePathType>::apply(std::forward<Tree1>(tree1),
                                                                                 std::forward<Tree2>(tree2),
                                                                                 std::forward<Visitor>(visitor));
        }

#else

        template<typename Tree1, typename Tree2, typename Visitor>
        void applyToTree(Tree1& tree1, Tree2& tree2, Visitor& visitor)
        {
          ApplyToTree<Visitor::treePathType>::apply(tree1,tree2,visitor);
        }

        template<typename Tree1, typename Tree2, typename Visitor>
        void applyToTree(const Tree1& tree1, const Tree2& tree2, Visitor& visitor)
        {
          ApplyToTree<Visitor::treePathType>::apply(tree1,tree2,visitor);
        }

        template<typename Tree1, typename Tree2, typename Visitor>
        void applyToTree(Tree1& tree1, Tree2& tree2, const Visitor& visitor)
        {
          ApplyToTree<Visitor::treePathType>::apply(tree1,tree2,visitor);
        }

        template<typename Tree1, typename Tree2, typename Visitor>
        void applyToTree(const Tree1& tree1, const Tree2& tree2, const Visitor& visitor)
        {
          ApplyToTree<Visitor::treePathType>::apply(tree1,tree2,visitor);
        }

#endif // HAVE_RVALUE_REFERENCES

        //! \} group TypeTree

      } // namespace PairTraversal
    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_PAIRTRAVERSAL_HH
