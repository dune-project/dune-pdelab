// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_COMPOSITENODE_HH
#define DUNE_PDELAB_COMMON_TYPETREE_COMPOSITENODE_HH

#include <dune/pdelab/common/nodetags.hh>
#include <dune/pdelab/common/utility.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/tuples.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Nodes
       *  \ingroup TypeTree
       *  \{
       */

      namespace {

        //! TMP for counting the actual number of children of the composite node.
        /**
         * This TMP counts the number of children that are not of type EmptyNode.
         * Moreover, it makes sure that a child of type EmptyNode is not followed by
         * a non-empty child.
         */
        template<typename Children, std::size_t i, std::size_t n, bool atEnd = false>
        struct count_children
        {

          static const bool emptyNode = is_same<typename tuple_element<i,Children>::type,EmptyNode>::value;

          dune_static_assert(atEnd ? emptyNode : true,"invalid child structure (EmptyNode followed by real node)");

          static const std::size_t value = count_children<Children,i+1,n,emptyNode>::value + (emptyNode ? 0 : 1);

        };

        //! End of TMP recursion
        template<typename Children, std::size_t n, bool atEnd>
        struct count_children<Children,n,n,atEnd>
        {

          static const std::size_t value = 0;

        };

        //! Pointer to an empty node that is used for all empty slots
        /**
         * TODO: move into a library!
         */
        shared_ptr<EmptyNode> emptyNodePtr(make_shared<EmptyNode>());

      } // anonymous namespace


      //! Implementation Helper for constructors of composite nodes.
      /**
       * Using this struct for all but the first constructor argument in
       * a composite node implementation makes it possible to only have
       * a single constructor regardless of the number of actual children
       * of the node.
       *
       * It should be used like this:
       *
       * \code
       * template<typename C1, typename C2, ...>
       * class MyCompositeNode {
       *   ...
       *   MyCompositeNode(C1& c1,
       *                   typename OptionalChild<C2>::type c2 = OptionalChild<C2>::default_value(),
       *                   ...)
       *     : BaseT(c1,c2,...)
       *   {}
       * };
       * \endcode
       */
      template<typename T>
      struct OptionalChild
      {
        //! The correct child type.
        typedef T& type;

        //! Method providing a default value for empty children.
        static T default_value()
        {
          dune_static_assert((<AlwaysFalse<T>::value), "You must provide a constructor parameter for non-empty children!");
          DUNE_THROW(InvalidArgument,"You must provide a constructor parameter for non-empty children!");
        }
      };

#ifndef DOXYGEN

      //! Specialization for empty children.
      template<>
      struct OptionalChild<EmptyNode>
      {
        typedef EmptyNode type;

        static EmptyNode default_value()
        {
          return EmptyNode();
        }
      };

#endif // DOXYGEN


      /** \brief Base class for composite nodes combining children of different types within a TypeTree.
       *
       * A CompositeNode can tie together up to 10 children of different types.
       *
       * \note If you need more than 10 children in a composite node and can use a compiler that supports
       * the upcoming C++0x standard, consider using a VariadicCompositeNode instead.
       *
       * \tparam C0,...,C9 The types of the children.
       */
      template<typename C0, typename C1 = EmptyNode, typename C2 = EmptyNode, typename C3 = EmptyNode, typename C4 = EmptyNode,
               typename C5 = EmptyNode, typename C6 = EmptyNode, typename C7 = EmptyNode, typename C8 = EmptyNode, typename C9 = EmptyNode>
      class CompositeNode
      {

      public:

        //! The type used for storing the children.
        typedef tuple<shared_ptr<C0>,
                      shared_ptr<C1>,
                      shared_ptr<C2>,
                      shared_ptr<C3>,
                      shared_ptr<C4>,
                      shared_ptr<C5>,
                      shared_ptr<C6>,
                      shared_ptr<C7>,
                      shared_ptr<C8>,
                      shared_ptr<C9>
                      > NodeStorage;

        //! The types of all children.
        typedef tuple<C0,C1,C2,C3,C4,C5,C6,C7,C8,C9> ChildTypes;

        //! Mark this class as non leaf in the TypeTree.
        static const bool isLeaf = false;

        //! Mark this class as a composite in the TypeTree.
        static const bool isComposite = true;

        //! Mark this class as a non power in the typeTree.
        static const bool isPower = false;

        //! The type tag that describes a CompositeNode.
        typedef CompositeNodeTag NodeTag;

#ifdef DOXYGEN
        //! The number of children of the CompositeNode.
        static const std::size_t CHILDREN = implementation-defined;
#else
        static const std::size_t CHILDREN = count_children<ChildTypes,0,tuple_size<ChildTypes>::value>::value;
#endif

        //! Access to the type and storage type of the i-th child.
        template<std::size_t k>
        struct Child {
          //! The type of the child.
          typedef typename tuple_element<k,ChildTypes>::type Type;

          //! The storage type of the child.
          typedef typename tuple_element<k,NodeStorage>::type Storage;

          //! The const storage type of the child.
          typedef shared_ptr<const typename tuple_element<k,ChildTypes>::type> ConstStorage;
        };


        //! @name Child Access
        //! @{

        //! Returns the i-th child (const version).
        /**
         * \returns a const reference to the i-th child.
         */
        template<std::size_t k>
        const typename Child<k>::Type& child() const
        {
          return *get<k>(_children);
        }

        //! Returns the i-th child (const version).
        /**
         * \returns a const reference to the i-th child.
         */
        template<std::size_t k>
        const typename Child<k>::Type& getChild() const
        {
          return child<k>();
        }

        //! Returns the i-th child.
        /**
         * \returns a reference to the i-th child.
         */
        template<std::size_t k>
        typename Child<k>::Type& child()
        {
          return *get<k>(_children);
        }

        //! Returns the i-th child.
        /**
         * \returns a reference to the i-th child.
         */
        template<std::size_t k>
        typename Child<k>::Type& getChild()
        {
          return child<k>();
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
          return get<k>(_children);
        }

        //! Returns the storage of the i-th child.
        /**
         * \returns a copy of the object storing the i-th child.
         */
        template<std::size_t k>
        typename Child<k>::Storage childStorage()
        {
          return get<k>(_children);
        }

        //! Sets the i-th child to the passed-in value.
        template<std::size_t k>
        void setChild(typename Child<k>::type& child)
        {
          get<k>(_children) = stackobject_to_shared_ptr(child);
        }

        //! Sets the stored value representing the i-th child to the passed-in value.
        template<std::size_t k>
        void setChild(typename Child<k>::storage_type child)
        {
          get<k>(_children) = child;
        }

        //! @}

      private:

        //! Helper function to correctly handle empty nodes in the constructor.
        /**
         * The default implementation assumes the passed-in object to be located on the stack
         * and wraps it in a shared_ptr with a no-op deleter.
         */
        template<typename T>
        static shared_ptr<T> guarded_wrap_object(T& t)
        {
          return stackobject_to_shared_ptr(t);
        }

        //! Helper function specialization for empty nodes.
        /**
         * For efficiency reasons, this returns a static shared_ptr, allowing the EmptyNode object
         * and its reference counter block to be shared by all empty children.
         */
        static shared_ptr<EmptyNode> guarded_wrap_object(EmptyNode& en)
        {
          return emptyNodePtr;
        }

      protected:

        //! Default constructor.
        /**
         * The default constructor is protected, as CompositeNode is a utility
         * class that needs to be filled with meaning by subclassing it
         * and adding useful functionality to the subclass.
         *
         * \warning When using the default constructor, make sure to set ALL children
         * by means of the setChild() methods!
         */
        CompositeNode()
        {}

        //! Initializes the CompositeNode with the passed-in child objects.
        CompositeNode(C0& c0,
                      typename OptionalChild<C1>::type c1 = typename OptionalChild<C1>::type(),
                      typename OptionalChild<C2>::type c2 = typename OptionalChild<C2>::type(),
                      typename OptionalChild<C3>::type c3 = typename OptionalChild<C3>::type(),
                      typename OptionalChild<C4>::type c4 = typename OptionalChild<C4>::type(),
                      typename OptionalChild<C5>::type c5 = typename OptionalChild<C5>::type(),
                      typename OptionalChild<C6>::type c6 = typename OptionalChild<C6>::type(),
                      typename OptionalChild<C7>::type c7 = typename OptionalChild<C7>::type(),
                      typename OptionalChild<C8>::type c8 = typename OptionalChild<C8>::type(),
                      typename OptionalChild<C9>::type c9 = typename OptionalChild<C9>::type())
          : _children(stackobject_to_shared_ptr(c0),
                      guarded_wrap_object(c1),
                      guarded_wrap_object(c2),
                      guarded_wrap_object(c3),
                      guarded_wrap_object(c4),
                      guarded_wrap_object(c5),
                      guarded_wrap_object(c6),
                      guarded_wrap_object(c7),
                      guarded_wrap_object(c8),
                      guarded_wrap_object(c9))
        {}

        //! Initializes the CompositeNode with copies of the passed-in storage objects.
        CompositeNode(shared_ptr<C0> c0,
                      shared_ptr<C1> c1,
                      shared_ptr<C2> c2,
                      shared_ptr<C3> c3,
                      shared_ptr<C4> c4,
                      shared_ptr<C5> c5,
                      shared_ptr<C6> c6,
                      shared_ptr<C7> c7,
                      shared_ptr<C8> c8,
                      shared_ptr<C9> c9)
          : _children(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
        {}

        //! Initializes the CompositeNode from the passed-in NodeStorage object.
        CompositeNode(const NodeStorage& children)
          : _children(children)
        {}

      private:
        NodeStorage _children;
      };

      //! \} group Nodes

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_COMPOSITENODE_HH
