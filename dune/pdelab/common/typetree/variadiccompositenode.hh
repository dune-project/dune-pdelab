// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_VARIADICCOMPOSITENODE_HH
#define DUNE_PDELAB_COMMON_TYPETREE_VARIADICCOMPOSITENODE_HH

#if !(HAVE_VARIADIC_TEMPLATES || DOXYGEN)
#error The class VariadicCompositeNode requires compiler support for variadic templates, which your compiler lacks.
#endif

#include <dune/pdelab/common/typetree/nodetags.hh>
#include <dune/common/tuples.hh>
#include <dune/common/deprecated.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Nodes
       *  \ingroup TypeTree
       *  \{
       */

      //! Base class for composite nodes based on variadic templates.
      template<typename... Children>
      class VariadicCompositeNode
      {

      public:

        //! The type tag that describes a VariadicCompositeNode.
        typedef VariadicCompositeNodeTag NodeTag;

        //! The type used for storing the children.
        typedef tuple<shared_ptr<Children>... > NodeStorage;

        //! A tuple storing the types of all children.
        typedef tuple<Children...> ChildTypes;

        //! Mark this class as non leaf in the \ref TypeTree.
        static const bool isLeaf = false;

        //! Mark this class as a non power in the \ref TypeTree.
        static const bool isPower = false;

        //! Mark this class as a composite in the \ref TypeTree.
        static const bool isComposite = true;

        //! The number of children.
        static const std::size_t CHILDREN = sizeof...(Children);

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

        //! Returns the i-th child.
        /**
         * \returns a reference to the i-th child.
         */
        template<std::size_t k>
        typename Child<k>::Type& child()
        {
          return *get<k>(_children);
        }

        //! Returns the i-th child (const version).
        /**
         * \returns a const reference to the i-th child.
         */
        template<std::size_t k>
        const typename Child<k>::Type& child() const
        {
          return *get<k>(_children);
        }

        //! Returns the i-th child.
        /**
         * \returns a reference to the i-th child.
         */
        template<std::size_t k>
        typename Child<k>::Type& DUNE_DEPRECATED getChild()
        {
          return child<k>();
        }

        //! Returns the i-th child (const version).
        /**
         * \returns a const reference to the i-th child.
         */
        template<std::size_t k>
        const typename Child<k>::Type& DUNE_DEPRECATED getChild() const
        {
          return child<k>();
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

        //! Sets the i-th child to the passed-in value.
        template<std::size_t k>
        void setChild(typename Child<k>::type& child)
        {
          get<k>(_children) = stackobject_to_shared_ptr(child);
        }

        //! Sets the storage of the i-th child to the passed-in value.
        template<std::size_t k>
        void setChild(typename Child<k>::storage_type child)
        {
          get<k>(_children) = child;
        }

        const NodeStorage& nodeStorage() const
        {
          return _children;
        }

        //! @}

      protected:

        //! @name Constructors
        //! @{

        //! Default constructor.
        /**
         * This constructor requires the storage type to be default
         * constructible.
         * \warning If the storage type is a pointer, the resulting object
         * will not be usable before its children are set using any of the
         * setChild(...) methods!
         */
        VariadicCompositeNode()
        {}

#if HAVE_RVALUE_REFERENCES && HAVE_VARIADIC_CONSTRUCTOR_SFINAE
        //! Initialize all children with the passed-in objects.
        template<typename... Args, typename = typename enable_if<(sizeof...(Args) == CHILDREN)>::type>
        VariadicCompositeNode(Args&&... args)
          : _children(convert_arg(std::forward<Args>(args))...)
        {}

#else

        VariadicCompositeNode(Children&... children)
          : _children(convert_arg(children)...)
        {}

#endif

        //! Initialize the VariadicCompositeNode with copies of the passed in Storage objects.
        VariadicCompositeNode(shared_ptr<Children>... children)
          : _children(children...)
        {}

        //! Initialize the VariadicCompositeNode with a copy of the passed-in storage type.
        VariadicCompositeNode(const NodeStorage& children)
          : _children(children)
        {}

        //! @}

      private:
        NodeStorage _children;
      };

      //! \} group Nodes

    } // namespace TypeTree

  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_VARIADICCOMPOSITENODE_HH
