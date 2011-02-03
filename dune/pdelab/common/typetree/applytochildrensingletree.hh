// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_APPLYTOCHILDRENSINGLETREE_HH
#define DUNE_PDELAB_COMMON_TYPETREE_APPLYTOCHILDRENSINGLETREE_HH

#include <dune/pdelab/common/typetree/nodetags.hh>
#include <dune/pdelab/common/typetree/treepath.hh>
#include <dune/pdelab/common/typetree/visitor.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Tree Traversal
       *  \ingroup TypeTree
       *  \{
       */

#ifndef DOXYGEN // these are all internals and not public API.

      namespace {

        // For the CompositeNode, we do need a TMP for iterating over the
        // children. Note that we use an index that counts down instead of up.
        // This allows us to specialize the TMP for the last child, where we
        // do not want to invoke the infix visitor on the CompositeNode.

        // There are two versions of this TMP, one for iteration with a static TreePath, and one
        // for iteration with a dynamic TreePath.




        // ********************************************************************************
        // Static Version
        // ********************************************************************************

        template<std::size_t inverse_k, std::size_t count>
        struct apply_to_children_fully_static
        {

#if HAVE_RVALUE_REFERENCES

          template<typename N, typename V, typename TreePath>
          static void apply(N&& n, V&& v, TreePath tp)
          {
            // make sure we do not try to work with references to the actual types
            typedef typename remove_reference<N>::type Node;
            typedef typename remove_reference<V>::type Visitor;

            // get child type
            typedef typename Node::template Child<count-inverse_k>::Type C;

            // extend TreePath by child index
            typedef typename TreePathPushBack<TreePath,count-inverse_k>::type ChildTreePath;

            // is the visitor interested in this child?
            const bool visit = Visitor::template VisitChild<Node,C,ChildTreePath>::value;

            // beforeChild() gets called regardless of the value of visit
            v.beforeChild(std::forward<N>(n),n.template child<count-inverse_k>(),tp,std::integral_constant<std::size_t,count-inverse_k>());

            // traverse to child
            ApplyToTree<Visitor::treePathType,typename C::NodeTag,visit>::apply(n.template child<count-inverse_k>(),
                                                                                std::forward<V>(v),
                                                                                ChildTreePath());

            // afterChild() gets called regardless of the value of visit
            v.afterChild(std::forward<N>(n),n.template child<count-inverse_k>(),tp,std::integral_constant<std::size_t,count-inverse_k>());

            // we are not at the last child (that is specialized), so call infix visitor callback
            v.in(std::forward<N>(n),tp);

            // continue with next child
            apply_to_children_fully_static<inverse_k-1,count>::apply(std::forward<N>(n),
                                                                     std::forward<V>(v),
                                                                     tp);
          }

#else

          // non-const tree, non-const visitor
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

          // const tree, non-const visitor
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

          // non-const tree, const visitor
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

          // const tree, const visitor
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

          // non-const tree, non-const visitor
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

          // const tree, non-const visitor
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

          // non-const tree, const visitor
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

          // const tree, const visitor
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




        // ********************************************************************************
        // Dynamic Version
        // ********************************************************************************

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


        // helper struct for automatically picking the right traversal
        // algorithm variant
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

      // The traversal algorithm is identical for CompositeNode and VariadicCompositeNode
      // (and even PowerNode for static traversal), so the implementation can be bundled
      // in a single base class.
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

#endif // DOXYGEN

      //! \} group Tree Traversal

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_APPLYTOCHILDRENSINGLETREE_HH
