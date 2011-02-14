// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_POWERNODE_HH
#define DUNE_PDELAB_COMMON_TYPETREE_POWERNODE_HH

#include <dune/common/deprecated.hh>
#include <dune/common/array.hh>
#include <dune/common/typetraits.hh>

#include <dune/pdelab/common/typetree/nodetags.hh>
#include <dune/pdelab/common/typetree/utility.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Nodes
       *  \ingroup TypeTree
       *  \{
       */

#if HAVE_VARIADIC_TEMPLATES && HAVE_RVALUE_REFERENCES

#ifndef DOXYGEN

      namespace {

        // prototype and end of recursion
        template<typename T, typename It, typename... Args>
        void assign_reference_pack_to_shared_ptr_array_unpack(It it, Args&&... args) {}

        template<typename T, typename It, typename Arg, typename... Args>
        void assign_reference_pack_to_shared_ptr_array_unpack(It it, Arg&& arg, Args&&... args)
        {
          static_assert(is_same<T,typename remove_const<typename remove_reference<Arg>::type>::type>::value,"type mismatch during array conversion");
          *it = convert_arg(std::forward<Arg>(arg));
          assign_reference_pack_to_shared_ptr_array_unpack<T>(++it,std::forward<Args>(args)...);
        }

        template<typename T, std::size_t n, typename... Args>
        void assign_reference_pack_to_shared_ptr_array(array<shared_ptr<T>,n>& res, Args&&... args)
        {
          static_assert(sizeof...(Args) == n, "invalid number of arguments");
          return assign_reference_pack_to_shared_ptr_array_unpack<T>(res.begin(),std::forward<Args>(args)...);
        }

      } // anonymous namespace

#endif

#endif

      template<typename PowerNode, typename T, std::size_t k>
      struct PowerNodeChildCountCheck
        : public enable_if<is_same<
                             typename PowerNode::ChildType,
                             T>::value &&
                           PowerNode::CHILDREN == k,
                           T>
      {};

      /** \brief Collect k instances of type T within a \ref TypeTree.
       *
       *  \tparam T The base type
       *  \tparam k The number of instances this node should collect
       */
      template<typename T, std::size_t k>
      class PowerNode
      {

      public:

        //! Mark this class as non leaf in the \ref TypeTree.
        static const bool isLeaf = false;

        //! Mark this class as a power in the \ref TypeTree.
        static const bool isPower = true;

        //! Mark this class as a non composite in the \ref TypeTree.
        static const bool isComposite = false;

        //! The number of children.
        static const std::size_t CHILDREN = k;

        //! The type tag that describes a PowerNode.
        typedef PowerNodeTag NodeTag;

        //! The type of each child.
        typedef T ChildType;

        //! The storage type of each child.
        typedef shared_ptr<T> ChildStorageType;

        //! The const version of the storage type of each child.
        typedef shared_ptr<const T> ChildConstStorageType;

        //! The type used for storing the children.
        typedef array<ChildStorageType,k> NodeStorage;


        //! Access to the type and storage type of the i-th child.
        template<std::size_t i>
        struct Child
        {
          //! The type of the child.
          typedef T Type;

          //! The storage type of the child.
          typedef ChildStorageType Storage;

          //! The const storage type of the child.
          typedef ChildConstStorageType ConstStorage;
        };

        //! @name Child Access (templated methods)
        //! @{

        //! Returns the i-th child.
        /**
         * \returns a reference to the i-th child.
         */
        template<std::size_t i>
        T& child ()
        {
          return *_children[i];
        }

        //! Returns the i-th child (const version).
        /**
         * \returns a const reference to the i-th child.
         */
        template<std::size_t i>
        const T& child () const
        {
          return *_children[i];
        }

        //! Returns the i-th child.
        /**
         * \returns a reference to the i-th child.
         */
        template<int i>
        T& DUNE_DEPRECATED getChild ()
        {
          return *_children[i];
        }

        //! Returns the i-th child (const version).
        /**
         * \returns a const reference to the i-th child.
         */
        template<int i>
        const T& DUNE_DEPRECATED getChild () const
        {
          return *_children[i];
        }

        //! Returns the storage of the i-th child.
        /**
         * \returns a copy of the object storing the i-th child.
         */
        template<int i>
        ChildStorageType childStorage()
        {
          return _children[i];
        }

        //! Returns the storage of the i-th child (const version).
        /**
         * This method is only important if the child is stored as
         * some kind of pointer, as this allows the pointee type to
         * become const.
         * \returns a copy of the object storing the i-th child.
         */
        template<int i>
        ChildConstStorageType childStorage() const
        {
          return _children[i];
        }

        //! Sets the i-th child to the passed-in value.
        template<int i>
        void setChild (T& t)
        {
          _children[i] = stackobject_to_shared_ptr(t);
        }

        //! Sets the stored value representing the i-th child to the passed-in value.
        template<int i>
        void setChild (ChildStorageType st)
        {
          _children[i] = st;
        }

        //! @}


        //! @name Child Access (Dynamic methods)
        //! @{

        //! Returns the i-th child.
        /**
         * \returns a reference to the i-th child.
         */
        T& child (int i)
        {
          return *_children[i];
        }

        //! Returns the i-th child (const version).
        /**
         * \returns a const reference to the i-th child.
         */
        const T& child (int i) const
        {
          return *_children[i];
        }

        //! Returns the i-th child.
        /**
         * \returns a reference to the i-th child.
         */
        T& getChild (int i) DUNE_DEPRECATED
        {
          return *_children[i];
        }

        //! Returns the i-th child (const version).
        /**
         * \returns a const reference to the i-th child.
         */
        const T& getChild (int i) const DUNE_DEPRECATED
        {
          return *_children[i];
        }

        //! Returns the storage of the i-th child.
        /**
         * \returns a copy of the object storing the i-th child.
         */
        ChildStorageType childStorage(int i)
        {
          return _children[i];
        }

        //! Returns the storage of the i-th child (const version).
        /**
         * This method is only important if the child is stored as
         * some kind of pointer, as this allows the pointee type to
         * become const.
         * \returns a copy of the object storing the i-th child.
         */
        ChildConstStorageType childStorage (int i) const
        {
          return (_children[i]);
        }

        //! Sets the i-th child to the passed-in value.
        void setChild (int i, T& t)
        {
          _children[i] = stackobject_to_shared_ptr(t);
        }

        //! Sets the stored value representing the i-th child to the passed-in value.
        void setChild (int i, ChildStorageType st)
        {
          _children[i] = st;
        }

        const NodeStorage& nodeStorage() const
        {
          return _children;
        }

        //! @}

        //! @name Constructors
        //! @{

      protected:

        //! Default constructor.
        /**
         * The default constructor is protected, as PowerNode is a utility
         * class that needs to be filled with meaning by subclassing it
         * and adding useful functionality to the subclass.
         *
         * \warning When using the default constructor, make sure to set ALL children
         * by means of the setChild() methods!
         */
        PowerNode()
        {}

        //! Initialize the PowerNode with a copy of the passed-in storage type.
        PowerNode(const NodeStorage& children)
          : _children(children)
        {}

        //! Initialize all children with copies of a storage object constructed from the parameter \c t.
        PowerNode (T& t, bool distinct_objects = true)
        {
          if (distinct_objects)
            {
              for (typename NodeStorage::iterator it = _children.begin(); it != _children.end(); ++it)
                *it = make_shared<T>(t);
            }
          else
            {
              shared_ptr<T> sp = stackobject_to_shared_ptr(t);
              std::fill(_children.begin(),_children.end(),sp);
            }
        }

#ifdef DOXYGEN

        //! Initialize all children with the passed-in objects.
        /**
         * The availability of this constructor depends on the number of children and
         * compiler support for C++0x: For 1 <= k <= 10, it is always present, but for
         * k > 10, it requires C++0x support in the compiler. If your compiler doesn't,
         * use PowerNode(const Storage& children) instead.
         *
         * Moreover, the C++0x-based version also supports passing in temporary objects
         * and will move those objects into the node. Attempting to do so with the legacy
         * version will result in a compile error.
         */
        PowerNode(T& t1, T& t2, ...)
        {}

#else

#if HAVE_VARIADIC_TEMPLATES && HAVE_RVALUE_REFERENCES

        template<typename... Children>
        PowerNode (Children&&... children)
        {
          assign_reference_pack_to_shared_ptr_array(_children,std::forward<Children>(children)...);
        }

#else

        template<typename U>
        PowerNode (typename PowerNodeChildCountCheck<PowerNode,U,2>::type& c0,
                   U& c1)
        {
          _children[0] = stackobject_to_shared_ptr(c0);
          _children[1] = stackobject_to_shared_ptr(c1);
        }

        template<typename U>
        PowerNode (typename PowerNodeChildCountCheck<PowerNode,U,3>::type& c0,
                   U& c1,
                   U& c2)
        {
          _children[0] = stackobject_to_shared_ptr(c0);
          _children[1] = stackobject_to_shared_ptr(c1);
          _children[2] = stackobject_to_shared_ptr(c2);
        }

        template<typename U>
        PowerNode (typename PowerNodeChildCountCheck<PowerNode,U,4>::type& c0,
                   U& c1,
                   U& c2,
                   U& c3)
        {
          _children[0] = stackobject_to_shared_ptr(c0);
          _children[1] = stackobject_to_shared_ptr(c1);
          _children[2] = stackobject_to_shared_ptr(c2);
          _children[3] = stackobject_to_shared_ptr(c3);
        }

        template<typename U>
        PowerNode (typename PowerNodeChildCountCheck<PowerNode,U,5>::type& c0,
                   U& c1,
                   U& c2,
                   U& c3,
                   U& c4)
        {
          _children[0] = stackobject_to_shared_ptr(c0);
          _children[1] = stackobject_to_shared_ptr(c1);
          _children[2] = stackobject_to_shared_ptr(c2);
          _children[3] = stackobject_to_shared_ptr(c3);
          _children[4] = stackobject_to_shared_ptr(c4);
        }

        template<typename U>
        PowerNode (typename PowerNodeChildCountCheck<PowerNode,U,6>::type& c0,
                   U& c1,
                   U& c2,
                   U& c3,
                   U& c4,
                   U& c5)
        {
          _children[0] = stackobject_to_shared_ptr(c0);
          _children[1] = stackobject_to_shared_ptr(c1);
          _children[2] = stackobject_to_shared_ptr(c2);
          _children[3] = stackobject_to_shared_ptr(c3);
          _children[4] = stackobject_to_shared_ptr(c4);
          _children[5] = stackobject_to_shared_ptr(c5);
        }

        template<typename U>
        PowerNode (typename PowerNodeChildCountCheck<PowerNode,U,7>::type& c0,
                   U& c1,
                   U& c2,
                   U& c3,
                   U& c4,
                   U& c5,
                   U& c6)
        {
          _children[0] = stackobject_to_shared_ptr(c0);
          _children[1] = stackobject_to_shared_ptr(c1);
          _children[2] = stackobject_to_shared_ptr(c2);
          _children[3] = stackobject_to_shared_ptr(c3);
          _children[4] = stackobject_to_shared_ptr(c4);
          _children[5] = stackobject_to_shared_ptr(c5);
          _children[6] = stackobject_to_shared_ptr(c6);
        }

        template<typename U>
        PowerNode (typename PowerNodeChildCountCheck<PowerNode,U,8>::type& c0,
                   U& c1,
                   U& c2,
                   U& c3,
                   U& c4,
                   U& c5,
                   U& c6,
                   U& c7)
        {
          _children[0] = stackobject_to_shared_ptr(c0);
          _children[1] = stackobject_to_shared_ptr(c1);
          _children[2] = stackobject_to_shared_ptr(c2);
          _children[3] = stackobject_to_shared_ptr(c3);
          _children[4] = stackobject_to_shared_ptr(c4);
          _children[5] = stackobject_to_shared_ptr(c5);
          _children[6] = stackobject_to_shared_ptr(c6);
          _children[7] = stackobject_to_shared_ptr(c7);
        }

        template<typename U>
        PowerNode (typename PowerNodeChildCountCheck<PowerNode,U,9>::type& c0,
                   U& c1,
                   U& c2,
                   U& c3,
                   U& c4,
                   U& c5,
                   U& c6,
                   U& c7,
                   U& c8)
        {
          _children[0] = stackobject_to_shared_ptr(c0);
          _children[1] = stackobject_to_shared_ptr(c1);
          _children[2] = stackobject_to_shared_ptr(c2);
          _children[3] = stackobject_to_shared_ptr(c3);
          _children[4] = stackobject_to_shared_ptr(c4);
          _children[5] = stackobject_to_shared_ptr(c5);
          _children[6] = stackobject_to_shared_ptr(c6);
          _children[7] = stackobject_to_shared_ptr(c7);
          _children[8] = stackobject_to_shared_ptr(c8);
        }

        template<typename U>
        PowerNode (typename PowerNodeChildCountCheck<PowerNode,U,10>::type& c0,
                   U& c1,
                   U& c2,
                   U& c3,
                   U& c4,
                   U& c5,
                   U& c6,
                   U& c7,
                   U& c8,
                   U& c9)
        {
          _children[0] = stackobject_to_shared_ptr(c0);
          _children[1] = stackobject_to_shared_ptr(c1);
          _children[2] = stackobject_to_shared_ptr(c2);
          _children[3] = stackobject_to_shared_ptr(c3);
          _children[4] = stackobject_to_shared_ptr(c4);
          _children[5] = stackobject_to_shared_ptr(c5);
          _children[6] = stackobject_to_shared_ptr(c6);
          _children[7] = stackobject_to_shared_ptr(c7);
          _children[8] = stackobject_to_shared_ptr(c8);
          _children[9] = stackobject_to_shared_ptr(c9);
        }

#endif // C++0x

#endif // DOXYGEN

        //! @}

      private:
        NodeStorage _children;
      };

      //! \} group Nodes

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_POWERNODE_HH
