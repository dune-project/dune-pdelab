// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_PROXYNODE_HH
#define DUNE_PDELAB_COMMON_TYPETREE_PROXYNODE_HH

#include <dune/pdelab/common/typetree/nodetags.hh>
#include <dune/common/shared_ptr.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Nodes
       *  \ingroup TypeTree
       *  \{
       */

      template<typename Node>
      class ProxyNode;

      //! Mixin class providing methods for child access with compile-time parameter.
      template<typename ProxiedNode>
      class StaticChildAccessors
      {

        static const bool proxiedNodeIsConst = IsConst<typename remove_reference<ProxiedNode>::type>::value;

        template<std::size_t k>
        struct lazy_enabled
        {
          static const bool value = !proxiedNodeIsConst;
        };

        typedef ProxyNode<ProxiedNode> Node;

        template<bool enabled = !proxiedNodeIsConst>
        typename enable_if<enabled,Node&>::type
        node()
        {
          return static_cast<Node&>(*this);
        }

        const Node& node() const
        {
          return static_cast<const Node&>(*this);
        }

      public:

        //! Access to the type and storage type of the i-th child.
        template<std::size_t k>
        struct Child
          : public ProxiedNode::template Child<k>
        {};

        //! @name Child Access
        //! @{

        //! Returns the i-th child.
        /**
         * \returns a reference to the i-th child.
         */
        template<std::size_t k>
        typename enable_if<lazy_enabled<k>::value,typename Child<k>::Type&>::type
        child()
        {
          return node().proxiedNode().template child<k>();
        }

        //! Returns the i-th child (const version).
        /**
         * \returns a const reference to the i-th child.
         */
        template<std::size_t k>
        const typename Child<k>::Type& child() const
        {
          return node().proxiedNode().template child<k>();
        }

        //! Returns the storage of the i-th child.
        /**
         * \returns a copy of the object storing the i-th child.
         */
        template<std::size_t k>
        typename enable_if<lazy_enabled<k>::value,typename Child<k>::Storage>::type
        childStorage()
        {
          return node().proxiedNode().template childStorage<k>();
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
          return node().proxiedNode().template childStorage<k>();
        }

        //! Sets the i-th child to the passed-in value.
        template<std::size_t k>
        void setChild(typename Child<k>::type& child, typename enable_if<lazy_enabled<k>::value,void*>::type = 0)
        {
          node().proxiedNode().template childStorage<k>() = stackobject_to_shared_ptr(child);
        }

        //! Sets the storage of the i-th child to the passed-in value.
        template<std::size_t k>
        void setChild(typename Child<k>::storage_type child, typename enable_if<lazy_enabled<k>::value,void*>::type = 0)
        {
           node().proxiedNode().template childStorage<k>() = child;
        }

        const typename ProxiedNode::NodeStorage& nodeStorage() const
        {
          return node().proxiedNode().nodeStorage();
        }

      };

      //! Mixin class providing methods for child access with run-time parameter.
      /**
       * This class also provides the compile-time parameter based methods, as
       * multiple inheritance from both DynamicChildAccessors and StaticChildAccessors
       * creates ambigous method lookups.
       */
      template<typename ProxiedNode>
      class DynamicChildAccessors
        : public StaticChildAccessors<ProxiedNode>
      {

        typedef ProxyNode<ProxiedNode> Node;

        static const bool proxiedNodeIsConst = IsConst<typename remove_reference<ProxiedNode>::type>::value;

        template<bool enabled = !proxiedNodeIsConst>
        typename enable_if<enabled,Node&>::type
        node()
        {
          return static_cast<Node&>(*this);
        }

        const Node& node() const
        {
          return static_cast<const Node&>(*this);
        }

      public:

        //! @name Child Access (Dynamic methods)
        //! @{

        //! Returns the i-th child.
        /**
         * \returns a reference to the i-th child.
         */
        template<bool enabled = !proxiedNodeIsConst>
        typename enable_if<enabled,typename ProxiedNode::ChildType&>::type
        child (std::size_t i)
        {
          return node().proxiedNode().child(i);
        }

        //! Returns the i-th child (const version).
        /**
         * \returns a const reference to the i-th child.
         */
        const typename ProxiedNode::ChildType& child (std::size_t i) const
        {
          return node().proxiedNode().child(i);
        }

        //! Returns the storage of the i-th child.
        /**
         * \returns a copy of the object storing the i-th child.
         */
        template<bool enabled = !proxiedNodeIsConst>
        typename enable_if<enabled,typename ProxiedNode::ChildStorageType>::type
        childStorage(std::size_t i)
        {
          return node().proxiedNode().childStorage(i);
        }

        //! Returns the storage of the i-th child (const version).
        /**
         * This method is only important if the child is stored as
         * some kind of pointer, as this allows the pointee type to
         * become const.
         * \returns a copy of the object storing the i-th child.
         */
        typename ProxiedNode::ChildConstStorageType childStorage (std::size_t i) const
        {
          return node().proxiedNode().childStorage(i);
        }

        //! Sets the i-th child to the passed-in value.
        template<bool enabled = !proxiedNodeIsConst>
        void setChild (std::size_t i, typename ProxiedNode::ChildType& t, typename enable_if<enabled,void*>::type = 0)
        {
          node().proxiedNode().childStorage(i) = stackobject_to_shared_ptr(t);
        }

        //! Sets the stored value representing the i-th child to the passed-in value.
        template<bool enabled = !proxiedNodeIsConst>
        void setChild (std::size_t i, typename ProxiedNode::ChildStorageType st, typename enable_if<enabled,void*>::type = 0)
        {
          node().proxiedNode().childStorage(i) = st;
        }

      };

      //! Tag-based dispatch to appropiate base class that provides necessary functionality.
      template<typename Node, typename NodeTag>
      struct ProxyNodeBase;

      //! ProxyNode base class for LeafNode.
      template<typename Node>
      struct ProxyNodeBase<Node,LeafNodeTag>
      {
      };

      //! ProxyNode base class for CompositeNode.
      template<typename Node>
      struct ProxyNodeBase<Node,CompositeNodeTag>
        : public StaticChildAccessors<Node>
      {
        typedef typename Node::ChildTypes ChildTypes;
        typedef typename Node::NodeStorage NodeStorage;
      };

      //! ProxyNode base class for VariadicCompositeNode.
      template<typename Node>
      struct ProxyNodeBase<Node,VariadicCompositeNodeTag>
        : public StaticChildAccessors<Node>
      {
        typedef typename Node::ChildTypes ChildTypes;
        typedef typename Node::NodeStorage NodeStorage;
      };

      //! ProxyNode base class for PowerNode.
      template<typename Node>
      struct ProxyNodeBase<Node,PowerNodeTag>
        : public DynamicChildAccessors<Node>
      {
        typedef typename Node::ChildType ChildType;
        typedef typename Node::NodeStorage NodeStorage;
      };


      //! Base class for nodes acting as a proxy for an existing node.
      /**
       * ProxyNode is a utility class for implementing proxy classes
       * that need to provide the TypeTree node functionality of the
       * proxied class. It exactly mirrors the TypeTree node characteristics
       * of the proxied node.
       */
      template<typename Node>
      class ProxyNode
        : public ProxyNodeBase<Node,typename Node::NodeTag>
      {

        static const bool proxiedNodeIsConst = IsConst<typename remove_reference<Node>::type>::value;

      public:

        typedef Node ProxiedNode;

        typedef typename Node::NodeTag NodeTag;

        //! Mark this class as non leaf in the \ref TypeTree.
        static const bool isLeaf = Node::isLeaf;

        //! Mark this class as a non power in the \ref TypeTree.
        static const bool isPower = Node::isPower;

        //! Mark this class as a composite in the \ref TypeTree.
        static const bool isComposite = Node::isComposite;

        //! The number of children.
        static const std::size_t CHILDREN = Node::CHILDREN;

        //! @name Access to the proxied node
        //! @{

        //! Returns the proxied node.
        template<bool enabled = !proxiedNodeIsConst>
        typename enable_if<enabled,Node&>::type
        proxiedNode()
        {
          return *_node;
        }

        //! Returns the proxied node (const version).
        const Node& proxiedNode() const
        {
          return *_node;
        }

        //! Returns the storage of the proxied node.
        template<bool enabled = !proxiedNodeIsConst>
        typename enable_if<enabled,shared_ptr<Node> >::type
        proxiedNodeStorage()
        {
          return _node;
        }

        //! Returns the storage of the proxied node (const version).
        shared_ptr<const Node> proxiedNodeStorage() const
        {
          return _node;
        }

        //! @}

      protected:

        //! @name Constructors
        //! @{

        ProxyNode(Node& node)
          : _node(stackobject_to_shared_ptr(node))
        {}

        ProxyNode(shared_ptr<Node> node)
          : _node(node)
        {}

        //! @}

      private:

        shared_ptr<Node> _node;
      };

      //! \} group Nodes

    } // namespace TypeTree

  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_VARIADICCOMPOSITENODE_HH
