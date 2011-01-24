// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_COMPOSITENODE_HH
#define DUNE_PDELAB_COMMON_TYPETREE_COMPOSITENODE_HH

#include <dune/pdelab/common/nodetags.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup TypeTree
       *  \ingroup PDELab
       *  \{
       */

      namespace {

        template<typename Children, std::size_t i, std::size_t n, bool atEnd = false>
        struct count_children
        {

          static const bool emptyNode = is_same<typename tuple_element<i,Children>::type,EmptyNode>::value;

          dune_static_assert(atEnd ? emptyNode : true,"invalid child structure (EmptyNode followed by real node)");

          static const std::size_t value = count_children<Children,i+1,n,emptyNode>::value + (emptyNode ? 0 : 1);

        };

        template<typename Children, std::size_t n, bool atEnd>
        struct count_children<Children,n,n,atEnd>
        {

          static const std::size_t value = 0;

        };

        shared_ptr<EmptyNode> emptyNodePtr(make_shared<EmptyNode>());

      } // anonymous namespace

      template<typename T>
      struct OptionalChild
      {
        typedef T& type;
      };

#ifndef DOXYGEN

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

      template<typename C0, typename C1 = EmptyNode, typename C2 = EmptyNode, typename C3 = EmptyNode, typename C4 = EmptyNode,
               typename C5 = EmptyNode, typename C6 = EmptyNode, typename C7 = EmptyNode, typename C8 = EmptyNode, typename C9 = EmptyNode>
      class CompositeNode
      {

      public:

        // type used for storing the children
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

        typedef tuple<C0,C1,C2,C3,C4,C5,C6,C7,C8,C9> ChildTypes;

        static const bool isLeaf = false;
        static const bool isComposite = true;
        static const bool isPower = false;
        typedef CompositeNodeTag NodeTag;

        static const std::size_t CHILDREN = count_children<ChildTypes,0,tuple_size<ChildTypes>::value>::value;

        // standard methods for child access
        template<std::size_t k>
        struct Child {
          typedef typename tuple_element<k,ChildTypes>::type Type;
          typedef typename tuple_element<k,NodeStorage>::type Storage;
          typedef shared_ptr<const typename tuple_element<k,ChildTypes>::type> ConstStorage;
        };

        template<std::size_t k>
        const typename Child<k>::Type& child() const
        {
          return *get<k>(_children);
        }

        template<std::size_t k>
        const typename Child<k>::Type& getChild() const
        {
          return child<k>();
        }

        template<std::size_t k>
        typename Child<k>::Type& child()
        {
          return *get<k>(_children);
        }

        template<std::size_t k>
        typename Child<k>::Type& getChild()
        {
          return child<k>();
        }

        template<std::size_t k>
        typename Child<k>::ConstStorage childStorage() const
        {
          return get<k>(_children);
        }

        template<std::size_t k>
        typename Child<k>::Storage childStorage()
        {
          return get<k>(_children);
        }

        template<std::size_t k>
        void setChild(typename Child<k>::type& child)
        {
          get<k>(_children) = stackobject_to_shared_ptr(child);
        }

        template<std::size_t k>
        void setChild(typename Child<k>::storage_type child)
        {
          get<k>(_children) = child;
        }

      private:

        template<typename T>
        static shared_ptr<T> guarded_wrap_object(T& t)
        {
          return stackobject_to_shared_ptr(t);
        }

        static shared_ptr<EmptyNode> guarded_wrap_object(EmptyNode& en)
        {
          return emptyNodePtr;
        }

      protected:

        CompositeNode()
        {}

        CompositeNode(C0& c0,
                      typename optional_child<C1>::type c1 = typename optional_child<C1>::type(),
                      typename optional_child<C2>::type c2 = typename optional_child<C2>::type(),
                      typename optional_child<C3>::type c3 = typename optional_child<C3>::type(),
                      typename optional_child<C4>::type c4 = typename optional_child<C4>::type(),
                      typename optional_child<C5>::type c5 = typename optional_child<C5>::type(),
                      typename optional_child<C6>::type c6 = typename optional_child<C6>::type(),
                      typename optional_child<C7>::type c7 = typename optional_child<C7>::type(),
                      typename optional_child<C8>::type c8 = typename optional_child<C8>::type(),
                      typename optional_child<C9>::type c9 = typename optional_child<C9>::type())
          : _children(wrap_stack_object(c0),
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

        CompositeNode(const NodeStorage& children)
          : _children(children)
        {}

      private:
        NodeStorage _children;
      };

      //! \} group TypeTree

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif DUNE_PDELAB_COMMON_TYPETREE_COMPOSITENODE_HH
