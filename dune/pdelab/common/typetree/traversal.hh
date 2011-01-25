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
      struct TypeTreeVisitor
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

#endif // HAVE_RVALUE_REFERENCES || DOXYGEN

      };


#ifndef DOXYGEN // these are all internals and not public API. Only access is using applyToTree().

      template<typename tag = StartTag>
      struct ApplyToTree
      {

#if HAVE_RVALUE_REFERENCES

        template<typename Node, typename Visitor>
        static void apply(Node&& node, Visitor&& visitor)
        {
          ApplyToTree<typename remove_reference<Node>::type::NodeTag>::apply(std::forward<Node>(node),
                                                                             std::forward<Visitor>(visitor),
                                                                             TreePath<>());
        }

#else

        template<typename Node, typename Visitor>
        static void apply(Node& node, Visitor& visitor)
        {
          ApplyToTree<typename Node::NodeTag>::apply(node,visitor,TreePath<>());
        }

        template<typename Node, typename Visitor>
        static void apply(const Node& node, Visitor& visitor)
        {
          ApplyToTree<typename Node::NodeTag>::apply(node,visitor,TreePath<>());
        }

        template<typename Node, typename Visitor>
        static void apply(Node& node, const Visitor& visitor)
        {
          ApplyToTree<typename Node::NodeTag>::apply(node,visitor,TreePath<>());
        }

        template<typename Node, typename Visitor>
        static void apply(const Node& node, const Visitor& visitor)
        {
          ApplyToTree<typename Node::NodeTag>::apply(node,visitor,TreePath<>());
        }

#endif // HAVE_RVALUE_REFERENCES

      };


      // LeafNode - again, this is easy: just do all three visits
      template<>
      struct ApplyToTree<LeafNodeTag>
      {

#if HAVE_RVALUE_REFERENCES

        template<typename N, typename V, typename TreePath>
        static void apply(N&& n, V&& v, TreePath tp)
        {
          v.leaf(std::forward<N>(n),tp);
        }

#else

        template<typename N, typename V, typename TreePath>
        static void apply(N& n, V& v, TreePath tp)
        {
          v.leaf(n,tp);
        }

        template<typename N, typename V, typename TreePath>
        static void apply(const N& n, V& v, TreePath tp)
        {
          v.leaf(n,tp);
        }

        template<typename N, typename V, typename TreePath>
        static void apply(N& n, const V& v, TreePath tp)
        {
          v.leaf(n,tp);
        }

        template<typename N, typename V, typename TreePath>
        static void apply(const N& n, const V& v, TreePath tp)
        {
          v.leaf(n,tp);
        }

#endif // HAVE_RVALUE_REFERENCES

      };

      namespace {

        // For the CompositeNode, we do again need a TMP for iterating over the
        // children. Note that we use an index that counts down instead of up.
        // This allows us to specialize the TMP for the last child, where we
        // do not want to invoke the infix visitor on the CompositeNode.
        template<std::size_t inverse_k, std::size_t count>
        struct apply_to_children
        {

#if HAVE_RVALUE_REFERENCES

          template<typename N, typename V, typename TreePath>
          static void apply(N&& n, V&& v, TreePath tp)
          {
            typedef typename remove_reference<N>::type::template Child<count-inverse_k>::Type C;
            ApplyToTree<typename C::NodeTag>::apply(n.template child<count-inverse_k>(),
                                                    std::forward<V>(v),
                                                    typename TreePathPushBack<TreePath,count-inverse_k>::type());
            v.in(std::forward<N>(n),tp);
            apply_to_children<inverse_k-1,count>::apply(std::forward<N>(n),
                                                        std::forward<V>(v),
                                                        tp);
          }

#else

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, V& v, TreePath tp)
          {
            typedef typename N::template Child<count-inverse_k>::Type C;
            ApplyToTree<typename C::NodeTag>::apply(n.template child<count-inverse_k>(),v,typename TreePathPushBack<TreePath,count-inverse_k>::type());
            v.in(n,tp);
            apply_to_children<inverse_k-1,count>::apply(n,v,tp);
          }

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, V& v, TreePath tp)
          {
            typedef typename N::template Child<count-inverse_k>::Type C;
            ApplyToTree<typename C::NodeTag>::apply(n.template child<count-inverse_k>(),v,typename TreePathPushBack<TreePath,count-inverse_k>::type());
            v.in(n,tp);
            apply_to_children<inverse_k-1,count>::apply(n,v,tp);
          }

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, const V& v, TreePath tp)
          {
            typedef typename N::template Child<count-inverse_k>::Type C;
            ApplyToTree<typename C::NodeTag>::apply(n.template child<count-inverse_k>(),v,typename TreePathPushBack<TreePath,count-inverse_k>::type());
            v.in(n,tp);
            apply_to_children<inverse_k-1,count>::apply(n,v,tp);
          }

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, const V& v, TreePath tp)
          {
            typedef typename N::template Child<count-inverse_k>::Type C;
            ApplyToTree<typename C::NodeTag>::apply(n.template child<count-inverse_k>(),v,typename TreePathPushBack<TreePath,count-inverse_k>::type());
            v.in(n,tp);
            apply_to_children<inverse_k-1,count>::apply(n,v,tp);
          }

#endif // HAVE_RVALUE_REFERENCES

        };

        // Specialization for last child. This specialization stops the recursion and
        // does not call the infix visitor on the CompositeNode.
        template<std::size_t count>
        struct apply_to_children<1,count>
        {

#if HAVE_RVALUE_REFERENCES

          template<typename N, typename V, typename TreePath>
          static void apply(N&& n, V&& v, TreePath tp)
          {
            typedef typename remove_reference<N>::type::template Child<count-1>::Type C;
            ApplyToTree<typename C::NodeTag>::apply(n.template child<count-1>(),
                                                      std::forward<V>(v),
                                                      typename TreePathPushBack<TreePath,count-1>::type());
          }

#else

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, V& v, TreePath tp)
          {
            typedef typename N::template Child<count-1>::Type C;
            ApplyToTree<typename C::NodeTag>::apply(n.template child<count-1>(),v,typename TreePathPushBack<TreePath,count-1>::type());
          }

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, V& v, TreePath tp)
          {
            typedef typename N::template Child<count-1>::Type C;
            ApplyToTree<typename C::NodeTag>::apply(n.template child<count-1>(),v,typename TreePathPushBack<TreePath,count-1>::type());
          }

          template<typename N, typename V, typename TreePath>
          static void apply(N& n, const V& v, TreePath tp)
          {
            typedef typename N::template Child<count-1>::Type C;
            ApplyToTree<typename C::NodeTag>::apply(n.template child<count-1>(),v,typename TreePathPushBack<TreePath,count-1>::type());
          }

          template<typename N, typename V, typename TreePath>
          static void apply(const N& n, const V& v, TreePath tp)
          {
            typedef typename N::template Child<count-1>::Type C;
            ApplyToTree<typename C::NodeTag>::apply(n.template child<count-1>(),v,typename TreePathPushBack<TreePath,count-1>::type());
          }

#endif // HAVE_RVALUE_REFERENCES
        };

        // Specialization for CompositeNode without any children.
        template<std::size_t count>
        struct apply_to_children<0,count>
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

      } // anonymous namespace

        // Base class for composite node traversal
      struct ApplyToGenericCompositeNode
      {

#if HAVE_RVALUE_REFERENCES

        template<typename N, typename V, typename TreePath>
        static void apply(N&& n, V&& v, TreePath tp)
        {
          v.pre(n,tp);
          typedef typename remove_reference<N>::type Node;
          apply_to_children<Node::CHILDREN,Node::CHILDREN>::apply(std::forward<N>(n),
                                                                  std::forward<V>(v),
                                                                  tp);
          v.post(n,tp);
        }

#else

        template<typename N, typename V, typename TreePath>
        static void apply(N& n, V& v, TreePath tp)
        {
          v.pre(n,tp);
          apply_to_children<N::CHILDREN,N::CHILDREN>::apply(n,v,tp);
          v.post(n,tp);
        }

        template<typename N, typename V, typename TreePath>
        static void apply(const N& n, V& v, TreePath tp)
        {
          v.pre(n,tp);
          apply_to_children<N::CHILDREN,N::CHILDREN>::apply(n,v,tp);
          v.post(n,tp);
        }

        template<typename N, typename V, typename TreePath>
        static void apply(N& n, const V& v, TreePath tp)
        {
          v.pre(n,tp);
          apply_to_children<N::CHILDREN,N::CHILDREN>::apply(n,v,tp);
          v.post(n,tp);
        }

        template<typename N, typename V, typename TreePath>
        static void apply(const N& n, const V& v, TreePath tp)
        {
          v.pre(n,tp);
          apply_to_children<N::CHILDREN,N::CHILDREN>::apply(n,v,tp);
          v.post(n,tp);
        }

#endif // HAVE_RVALUE_REFERENCES

      };

      // Traverse PowerNode
      template<>
      struct ApplyToTree<PowerNodeTag>
        : public ApplyToGenericCompositeNode
      {

        using ApplyToGenericCompositeNode::apply;

      };


      // Traverse CompositeNode
      template<>
      struct ApplyToTree<CompositeNodeTag>
        : public ApplyToGenericCompositeNode
      {

        using ApplyToGenericCompositeNode::apply;

      };


      // Traverse VariadicCompositeNode
      template<>
      struct ApplyToTree<VariadicCompositeNodeTag>
        : public ApplyToGenericCompositeNode
      {

        using ApplyToGenericCompositeNode::apply;

      };

#endif // DOXYGEN

#if HAVE_RVALUE_REFERENCES

      template<typename Tree, typename Visitor>
      void applyToTree(Tree&& tree, Visitor&& visitor)
      {
        ApplyToTree<>::apply(std::forward<Tree>(tree),
                             std::forward<Visitor>(visitor));
      }

#else

      template<typename Tree, typename Visitor>
      void applyToTree(Tree& tree, Visitor& visitor)
      {
        ApplyToTree<>::apply(tree,visitor);
      }

      template<typename Tree, typename Visitor>
      void applyToTree(const Tree& tree, Visitor& visitor)
      {
        ApplyToTree<>::apply(tree,visitor);
      }

      template<typename Tree, typename Visitor>
      void applyToTree(Tree& tree, const Visitor& visitor)
      {
        ApplyToTree<>::apply(tree,visitor);
      }

      template<typename Tree, typename Visitor>
      void applyToTree(const Tree& tree, const Visitor& visitor)
      {
        ApplyToTree<>::apply(tree,visitor);
      }

#endif // HAVE_RVALUE_REFERENCES


    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_TRAVERSAL_HH
