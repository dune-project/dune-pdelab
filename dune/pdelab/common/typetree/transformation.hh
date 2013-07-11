// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_TRANSFORMATION_HH
#define DUNE_PDELAB_COMMON_TYPETREE_TRANSFORMATION_HH

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/array.hh>
#include <dune/common/tuples.hh>
#include <dune/pdelab/common/typetraits.hh>
#include <dune/pdelab/common/typetree/nodetags.hh>
#include <dune/pdelab/common/typetree/utility.hh>


namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Transformation
       *  \ingroup TypeTree
       *  \{
       */

#ifdef DOXYGEN

      //! Look up transformation descriptor to transform SourceNode with Transformation.
      /**
       * The tree transformation engine expects this function to return a struct describing
       * how to perform the Transformation for The type SourceNode, which has ImplementationTag Tag.
       * This function has to be specialized for every combination of Transformation and Tag that
       * the transformation engine should support.
       *
       * \note The specialization does not have to placed in the namespace Dune::PDELab::TypeTree,
       *       it can simply reside in the same namespace as either the SourceNode or the Tag.
       *
       * \note This function will never be really called, the engine only extracts the return type.
       *       It is thus not necessary to actually implement the function, it is sufficient to
       *       declare it.
       *
       * \tparam SourceNode     The type of the node in the source tree that should be transformed.
       * \tparam Transformation The type of the transformation to apply to the source node.
       * \tparam Tag            The implementation tag of the source node.
       */
      template<typename SourceNode, typename Transformation, typename Tag>
      void lookupNodeTransformation(SourceNode* s, Transformation* t, Tag tag);

#else // DOXYGEN

      /**
       * \tparam S   C++ type of source node.
       * \tparam T   Tag identifying the transformation.
       * \tparam Tag Tag identifying the source type.
       *
       * Tag may be identical for different implementation of the same concept
       * (i.e. all leaf GridFunctionSpace), but this is not required.  This
       * allows you to handle different leaf GridFunctionSpace implementation
       * differently.  Tag should be extracted from S::ImplementationTag.
       */
      template<typename S, typename T, typename Tag>
      struct LookupNodeTransformation
      {
        // TODO: add configure test and replace __typeof__ with a macro
        typedef __typeof__(lookupNodeTransformation(declptr<S>(),declptr<T>(),Tag())) type;
        dune_static_assert((!is_same<type,void>::value), "Unable to find valid transformation descriptor");
      };

      struct EmptyNodeTransformation;

      // Specialization for EmptyNode. This is mainly here to save the user from possible
      // ambiguities when looking up lookupNodeTransformation().
      template<typename S, typename T>
      struct LookupNodeTransformation<S,T,EmptyNodeTag>
      {
        typedef EmptyNodeTransformation type;
      };

#endif // DOXYGEN


      //! Transform a TypeTree.
      /**
       * This struct can be used to apply a transformation to a given TypeTree. It exports the type of
       * the resulting (transformed) tree and contains methods to actually transform tree instances.
       *
       * \tparam SourceTree     = The TypeTree that should be transformed.
       * \tparam Transformation = The Transformation to apply to the TypeTree.
       * \tparam Tag            = This parameter is an implementation detail and must always be set to its default value.
       * \tparam recursive      = This parameter is an implementation detail and must always be set to its default value.
       */
      template<typename SourceTree, typename Transformation, typename Tag = StartTag, bool recursive = true>
      struct TransformTree
      {

#ifndef DOXYGEN

        typedef typename LookupNodeTransformation<SourceTree,Transformation,typename SourceTree::ImplementationTag>::type NodeTransformation;

        // the type of the new tree that will result from this transformation
        typedef typename TransformTree<SourceTree,Transformation,typename SourceTree::NodeTag,NodeTransformation::recursive>::transformed_type transformed_type;

        // the storage type of the new tree that will result from this transformation
        typedef typename TransformTree<SourceTree,Transformation,typename SourceTree::NodeTag,NodeTransformation::recursive>::transformed_storage_type transformed_storage_type;

#endif // DOXYGEN

        //! The type of the transformed tree.
        typedef transformed_type Type;


        //! Apply transformation to an existing tree s.
        static transformed_type transform(const SourceTree& s, const Transformation& t = Transformation())
        {
          return TransformTree<SourceTree,Transformation,typename SourceTree::NodeTag,NodeTransformation::recursive>::transform(s,t);
        }

        //! Apply transformation to an existing tree s.
        static transformed_type transform(const SourceTree& s, Transformation& t)
        {
          return TransformTree<SourceTree,Transformation,typename SourceTree::NodeTag,NodeTransformation::recursive>::transform(s,t);
        }

        //! Apply transformation to storage type of an existing tree, returning a heap-allocated storage type
        //! instance of the transformed tree.
        static transformed_storage_type transform_storage(shared_ptr<const SourceTree> sp, Transformation& t)
        {
          return TransformTree<SourceTree,Transformation,typename SourceTree::NodeTag,NodeTransformation::recursive>::transform_storage(sp,t);
        }

        //! Apply transformation to storage type of an existing tree, returning a heap-allocated storage type
        //! instance of the transformed tree.
        static transformed_storage_type transform_storage(shared_ptr<const SourceTree> sp, const Transformation& t)
        {
          return TransformTree<SourceTree,Transformation,typename SourceTree::NodeTag,NodeTransformation::recursive>::transform_storage(sp,t);
        }


      };

#ifndef DOXYGEN // internal per-node implementations of the transformation algorithm

      // handle a leaf node - this is easy
      template<typename S, typename T, bool recursive>
      struct TransformTree<S,T,LeafNodeTag,recursive>
      {
        // get transformed type from specification
        typedef typename LookupNodeTransformation<S,T,typename S::ImplementationTag>::type NodeTransformation;

        typedef typename NodeTransformation::transformed_type transformed_type;
        typedef typename NodeTransformation::transformed_storage_type transformed_storage_type;

        // delegate instance transformation to per-node specification
        static transformed_type transform(const S& s, T& t)
        {
          return NodeTransformation::transform(s,t);
        }

        // delegate instance transformation to per-node specification
        static transformed_type transform(const S& s, const T& t)
        {
          return NodeTransformation::transform(s,t);
        }

        static transformed_storage_type transform_storage(shared_ptr<const S> sp, T& t)
        {
          return NodeTransformation::transform_storage(sp,t);
        }

        static transformed_storage_type transform_storage(shared_ptr<const S> sp, const T& t)
        {
          return NodeTransformation::transform_storage(sp,t);
        }

      };


      // common implementation for non-recursive transformation of non-leaf nodes
      template<typename S, typename T>
      struct TransformTreeNonRecursive
      {
        // get transformed type from specification
        typedef typename LookupNodeTransformation<S,T,typename S::ImplementationTag>::type NodeTransformation;

        typedef typename NodeTransformation::transformed_type transformed_type;
        typedef typename NodeTransformation::transformed_storage_type transformed_storage_type;

        // delegate instance transformation to per-node specification
        static transformed_type transform(const S& s, T& t)
        {
          return NodeTransformation::transform(s,t);
        }

        // delegate instance transformation to per-node specification
        static transformed_type transform(const S& s, const T& t)
        {
          return NodeTransformation::transform(s,t);
        }

        static transformed_storage_type transform_storage(shared_ptr<const S> sp, T& t)
        {
          return NodeTransformation::transform_storage(sp,t);
        }

        static transformed_storage_type transform_storage(shared_ptr<const S> sp, const T& t)
        {
          return NodeTransformation::transform_storage(sp,t);
        }

      };


      // handle power tag - a little more tricky
      template<typename S, typename T>
      struct TransformTree<S,T,PowerNodeTag,true>
      {
        // get transformed type from specification
        // Handling this transformation in a way that makes the per-node specification easy to write
        // is a little involved:
        // The problem is that the transformed power node must be parameterized on the transformed child
        // type. So we need to transform the child type and pass the transformed child type to an inner
        // template of the node transformation struct called result (see example of such a specification
        // further down).
        typedef typename LookupNodeTransformation<S,T,typename S::ImplementationTag>::type NodeTransformation;
        typedef typename LookupNodeTransformation<typename S::ChildType,T,typename S::ChildType::ImplementationTag>::type ChildNodeTransformation;

        typedef typename NodeTransformation::template result<typename TransformTree<typename S::ChildType,
                                                                                    T,
                                                                                    typename S::ChildType::NodeTag,
                                                                                    ChildNodeTransformation::recursive>::transformed_type
                                                             >::type transformed_type;

        typedef typename NodeTransformation::template result<typename TransformTree<typename S::ChildType,
                                                                                    T,
                                                                                    typename S::ChildType::NodeTag,
                                                                                    ChildNodeTransformation::recursive>::transformed_type
                                                             >::storage_type transformed_storage_type;

        // Transform an instance of S.
        static transformed_type transform(const S& s, T& t)
        {
          // transform children
          typedef TransformTree<typename S::ChildType,T,typename S::ChildType::NodeTag,ChildNodeTransformation::recursive> ChildTreeTransformation;
          typedef typename ChildTreeTransformation::transformed_type transformed_child;
          const std::size_t child_count = S::CHILDREN;
          array<shared_ptr<transformed_child>,child_count> children;
          for (std::size_t k = 0; k < child_count; ++k) {
            children[k] = ChildTreeTransformation::transform_storage(s.childStorage(k),t);
          }
          // transform node
          return NodeTransformation::transform(s,t,children);
        }

        static transformed_type transform(const S& s, const T& t)
        {
          // transform children
          typedef TransformTree<typename S::ChildType,T,typename S::ChildType::NodeTag,ChildNodeTransformation::recursive> ChildTreeTransformation;
          typedef typename ChildTreeTransformation::transformed_type transformed_child;
          const std::size_t child_count = S::CHILDREN;
          array<shared_ptr<transformed_child>,child_count> children;
          for (std::size_t k = 0; k < child_count; ++k) {
            children[k] = ChildTreeTransformation::transform_storage(s.childStorage(k),t);
          }
          // transform node
          return NodeTransformation::transform(s,t,children);
        }

        static transformed_storage_type transform_storage(shared_ptr<const S> sp, T& t)
        {
          // transform children
          typedef TransformTree<typename S::ChildType,T,typename S::ChildType::NodeTag,ChildNodeTransformation::recursive> ChildTreeTransformation;
          typedef typename ChildTreeTransformation::transformed_storage_type transformed_child_storage;
          const std::size_t child_count = S::CHILDREN;
          array<transformed_child_storage,child_count> children;
          for (std::size_t k = 0; k < child_count; ++k) {
            children[k] = ChildTreeTransformation::transform_storage(sp->childStorage(k),t);
          }
          return NodeTransformation::transform_storage(sp,t,children);
        }

        static transformed_storage_type transform_storage(shared_ptr<const S> sp, const T& t)
        {
          // transform children
          typedef TransformTree<typename S::ChildType,T,typename S::ChildType::NodeTag,ChildNodeTransformation::recursive> ChildTreeTransformation;
          typedef typename ChildTreeTransformation::transformed_storage_type transformed_child_storage;
          const std::size_t child_count = S::CHILDREN;
          array<transformed_child_storage,child_count> children;
          for (std::size_t k = 0; k < child_count; ++k) {
            children[k] = ChildTreeTransformation::transform_storage(sp->childStorage(k),t);
          }
          return NodeTransformation::transform_storage(sp,t,children);
        }

      };

      // non-recursive version of the PowerNode transformation.
      template<typename S, typename T>
      struct TransformTree<S,T,PowerNodeTag,false>
        : public TransformTreeNonRecursive<S,T>
      {};



      // non-variadic version

      // helper struct that does the actual transformation for a composite node. We need this additional struct
      // to extract the template argument list with the types of all children from the node, which we cannot do
      // directly in the transformation<> template, as the type passed to transformation<> will usually be a
      // derived type and will normally have more template arguments than just the children. This declaration
      // just introduces the type of the helper struct, we always instantiate the specialization defined below;
      template<typename S, typename Children, typename T>
      struct transform_composite_node;


      // specialized version of the helper struct which extracts the template argument list with the children from
      // its second template parameter, which has to be CompositeNode::ChildTypes. Apart from that, the struct is
      // similar to the one for a PowerNode, but it obviously delegates transformation of the children to the TMP.
      template<typename S, typename T, typename C0, typename C1, typename C2, typename C3, typename C4, typename C5, typename C6, typename C7, typename C8, typename C9>
      struct transform_composite_node<S,tuple<C0,C1,C2,C3,C4,C5,C6,C7,C8,C9>,T>
      {

        // transformed type, using the same nested struct trick as the PowerNode
        typedef typename LookupNodeTransformation<S,T,typename S::ImplementationTag>::type NodeTransformation;

        typedef typename NodeTransformation::template result<
          typename TransformTree<C0,T,typename C0::NodeTag,LookupNodeTransformation<C0,T,typename C0::ImplementationTag>::type::recursive>::transformed_type,
          typename TransformTree<C1,T,typename C1::NodeTag,LookupNodeTransformation<C1,T,typename C1::ImplementationTag>::type::recursive>::transformed_type,
          typename TransformTree<C2,T,typename C2::NodeTag,LookupNodeTransformation<C2,T,typename C2::ImplementationTag>::type::recursive>::transformed_type,
          typename TransformTree<C3,T,typename C3::NodeTag,LookupNodeTransformation<C3,T,typename C3::ImplementationTag>::type::recursive>::transformed_type,
          typename TransformTree<C4,T,typename C4::NodeTag,LookupNodeTransformation<C4,T,typename C4::ImplementationTag>::type::recursive>::transformed_type,
          typename TransformTree<C5,T,typename C5::NodeTag,LookupNodeTransformation<C5,T,typename C5::ImplementationTag>::type::recursive>::transformed_type,
          typename TransformTree<C6,T,typename C6::NodeTag,LookupNodeTransformation<C6,T,typename C6::ImplementationTag>::type::recursive>::transformed_type,
          typename TransformTree<C7,T,typename C7::NodeTag,LookupNodeTransformation<C7,T,typename C7::ImplementationTag>::type::recursive>::transformed_type,
          typename TransformTree<C8,T,typename C8::NodeTag,LookupNodeTransformation<C8,T,typename C8::ImplementationTag>::type::recursive>::transformed_type,
          typename TransformTree<C9,T,typename C9::NodeTag,LookupNodeTransformation<C9,T,typename C9::ImplementationTag>::type::recursive>::transformed_type
          > resulttypes;

        typedef typename resulttypes::type transformed_type;
        typedef typename resulttypes::storage_type transformed_storage_type;
        typedef typename S::ImplementationTag Tag;

        static transformed_type transform(const S& s, T& t)
        {
          return NodeTransformation::transform(s,
                                               t,
                                               TransformTree<C0,T,typename C0::NodeTag,LookupNodeTransformation<C0,T,typename C0::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<0>(),t),
                                               TransformTree<C1,T,typename C1::NodeTag,LookupNodeTransformation<C1,T,typename C1::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<1>(),t),
                                               TransformTree<C2,T,typename C2::NodeTag,LookupNodeTransformation<C2,T,typename C2::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<2>(),t),
                                               TransformTree<C3,T,typename C3::NodeTag,LookupNodeTransformation<C3,T,typename C3::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<3>(),t),
                                               TransformTree<C4,T,typename C4::NodeTag,LookupNodeTransformation<C4,T,typename C4::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<4>(),t),
                                               TransformTree<C5,T,typename C5::NodeTag,LookupNodeTransformation<C5,T,typename C5::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<5>(),t),
                                               TransformTree<C6,T,typename C6::NodeTag,LookupNodeTransformation<C6,T,typename C6::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<6>(),t),
                                               TransformTree<C7,T,typename C7::NodeTag,LookupNodeTransformation<C7,T,typename C7::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<7>(),t),
                                               TransformTree<C8,T,typename C8::NodeTag,LookupNodeTransformation<C8,T,typename C8::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<8>(),t),
                                               TransformTree<C9,T,typename C9::NodeTag,LookupNodeTransformation<C9,T,typename C9::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<9>(),t));
        }

        static transformed_type transform(const S& s, const T& t)
        {
          return NodeTransformation::transform(s,
                                               t,
                                               TransformTree<C0,T,typename C0::NodeTag,LookupNodeTransformation<C0,T,typename C0::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<0>(),t),
                                               TransformTree<C1,T,typename C1::NodeTag,LookupNodeTransformation<C1,T,typename C1::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<1>(),t),
                                               TransformTree<C2,T,typename C2::NodeTag,LookupNodeTransformation<C2,T,typename C2::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<2>(),t),
                                               TransformTree<C3,T,typename C3::NodeTag,LookupNodeTransformation<C3,T,typename C3::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<3>(),t),
                                               TransformTree<C4,T,typename C4::NodeTag,LookupNodeTransformation<C4,T,typename C4::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<4>(),t),
                                               TransformTree<C5,T,typename C5::NodeTag,LookupNodeTransformation<C5,T,typename C5::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<5>(),t),
                                               TransformTree<C6,T,typename C6::NodeTag,LookupNodeTransformation<C6,T,typename C6::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<6>(),t),
                                               TransformTree<C7,T,typename C7::NodeTag,LookupNodeTransformation<C7,T,typename C7::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<7>(),t),
                                               TransformTree<C8,T,typename C8::NodeTag,LookupNodeTransformation<C8,T,typename C8::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<8>(),t),
                                               TransformTree<C9,T,typename C9::NodeTag,LookupNodeTransformation<C9,T,typename C9::ImplementationTag>::type::recursive>::transform_storage(s.template childStorage<9>(),t));
        }

        static transformed_storage_type transform_storage(shared_ptr<const S> sp, T& t)
        {
          return NodeTransformation::transform_storage(sp,
                                                       t,
                                                       TransformTree<C0,T,typename C0::NodeTag,LookupNodeTransformation<C0,T,typename C0::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<0>(),t),
                                                       TransformTree<C1,T,typename C1::NodeTag,LookupNodeTransformation<C1,T,typename C1::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<1>(),t),
                                                       TransformTree<C2,T,typename C2::NodeTag,LookupNodeTransformation<C2,T,typename C2::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<2>(),t),
                                                       TransformTree<C3,T,typename C3::NodeTag,LookupNodeTransformation<C3,T,typename C3::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<3>(),t),
                                                       TransformTree<C4,T,typename C4::NodeTag,LookupNodeTransformation<C4,T,typename C4::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<4>(),t),
                                                       TransformTree<C5,T,typename C5::NodeTag,LookupNodeTransformation<C5,T,typename C5::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<5>(),t),
                                                       TransformTree<C6,T,typename C6::NodeTag,LookupNodeTransformation<C6,T,typename C6::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<6>(),t),
                                                       TransformTree<C7,T,typename C7::NodeTag,LookupNodeTransformation<C7,T,typename C7::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<7>(),t),
                                                       TransformTree<C8,T,typename C8::NodeTag,LookupNodeTransformation<C8,T,typename C8::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<8>(),t),
                                                       TransformTree<C9,T,typename C9::NodeTag,LookupNodeTransformation<C9,T,typename C9::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<9>(),t));
        }

        static transformed_storage_type transform_storage(shared_ptr<const S> sp, const T& t)
        {
          return NodeTransformation::transform_storage(sp,
                                                       t,
                                                       TransformTree<C0,T,typename C0::NodeTag,LookupNodeTransformation<C0,T,typename C0::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<0>(),t),
                                                       TransformTree<C1,T,typename C1::NodeTag,LookupNodeTransformation<C1,T,typename C1::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<1>(),t),
                                                       TransformTree<C2,T,typename C2::NodeTag,LookupNodeTransformation<C2,T,typename C2::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<2>(),t),
                                                       TransformTree<C3,T,typename C3::NodeTag,LookupNodeTransformation<C3,T,typename C3::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<3>(),t),
                                                       TransformTree<C4,T,typename C4::NodeTag,LookupNodeTransformation<C4,T,typename C4::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<4>(),t),
                                                       TransformTree<C5,T,typename C5::NodeTag,LookupNodeTransformation<C5,T,typename C5::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<5>(),t),
                                                       TransformTree<C6,T,typename C6::NodeTag,LookupNodeTransformation<C6,T,typename C6::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<6>(),t),
                                                       TransformTree<C7,T,typename C7::NodeTag,LookupNodeTransformation<C7,T,typename C7::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<7>(),t),
                                                       TransformTree<C8,T,typename C8::NodeTag,LookupNodeTransformation<C8,T,typename C8::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<8>(),t),
                                                       TransformTree<C9,T,typename C9::NodeTag,LookupNodeTransformation<C9,T,typename C9::ImplementationTag>::type::recursive>::transform_storage(sp->template childStorage<9>(),t));
        }
      };


      // the specialization of transformation<> for the CompositeNode. This just extracts the
      // CompositeNode::ChildType member and forwards to the helper struct
      template<typename S, typename T>
      struct TransformTree<S,T,CompositeNodeTag,true>
      {
        typedef typename transform_composite_node<S,typename S::ChildTypes,T>::transformed_type transformed_type;
        typedef typename transform_composite_node<S,typename S::ChildTypes,T>::transformed_storage_type transformed_storage_type;

        static transformed_type transform(const S& s, T& t)
        {
          return transform_composite_node<S,typename S::ChildTypes,T>::transform(s,t);
        }

        static transformed_type transform(const S& s, const T& t)
        {
          return transform_composite_node<S,typename S::ChildTypes,T>::transform(s,t);
        }

        static transformed_storage_type transform_storage(shared_ptr<const S> sp, T& t)
        {
          return transform_composite_node<S,typename S::ChildTypes,T>::transform_storage(sp,t);
        }

        static transformed_storage_type transform_storage(shared_ptr<const S> sp, const T& t)
        {
          return transform_composite_node<S,typename S::ChildTypes,T>::transform_storage(sp,t);
        }

      };

      // non-recursive version of the CompositeNode transformation.
      template<typename S, typename T>
      struct TransformTree<S,T,CompositeNodeTag,false>
        : public TransformTreeNonRecursive<S,T>
      {};

#if HAVE_VARIADIC_TEMPLATES

      // helper struct that does the actual transformation for a composite node. We need this additional struct
      // to extract the template argument list with the types of all children from the node, which we cannot do
      // directly in the transformation<> template, as the type passed to transformation<> will usually be a
      // derived type and will normally have more template arguments than just the children. This declaration
      // just introduces the type of the helper struct, we always instantiate the specialization defined below;
      template<typename S, typename Children, typename T>
      struct transform_variadic_composite_node;

      // specialized version of the helper struct which extracts the template argument list with the children from
      // its second template parameter, which has to be CompositeNode::ChildTypes. Apart from that, the struct is
      // similar to the one for a PowerNode, but it obviously delegates transformation of the children to the TMP.
      template<typename S, typename T, typename... C>
      struct transform_variadic_composite_node<S,tuple<C...>,T>
      {

        // transformed type, using the same nested struct trick as the PowerNode
        typedef typename S::ImplementationTag Tag;
        typedef typename LookupNodeTransformation<S,T,Tag>::type NodeTransformation;
        typedef typename NodeTransformation::template result<typename TransformTree<C,
                                                                                    T,
                                                                                    typename C::NodeTag,
                                                                                    LookupNodeTransformation<C,T,typename C::ImplementationTag>::type::recursive
                                                                                    >::transformed_type...
                                                             >::type transformed_type;

        typedef typename NodeTransformation::template result<typename TransformTree<C,
                                                                                    T,
                                                                                    typename C::NodeTag,
                                                                                    LookupNodeTransformation<C,T,typename C::ImplementationTag>::type::recursive
                                                                                    >::transformed_type...
                                                             >::storage_type transformed_storage_type;

        // Retrieve the transformation descriptor for the child with index i.
        // This little helper improves really improves the readability of the
        // transformation functions.
        template<std::size_t i>
        struct ChildTransformation
          : public TransformTree<typename S::template Child<i>::Type,
                                 T,
                                 typename S::template Child<i>::Type::NodeTag,
                                 LookupNodeTransformation<
                                   typename S::template Child<i>::Type,
                                   T,
                                   typename S::template Child<i>::Type::ImplementationTag
                                   >::type::recursive
                                 >
        {};


        template<std::size_t... i>
        static transformed_type transform(const S& s, T& t, index_pack<i...> indices)
        {
          return NodeTransformation::transform(s,t,ChildTransformation<i>::transform_storage(s.template childStorage<i>(),t)...);
        }

        template<std::size_t... i>
        static transformed_type transform(const S& s, const T& t, index_pack<i...> indices)
        {
          return NodeTransformation::transform(s,t,ChildTransformation<i>::transform_storage(s.template childStorage<i>(),t)...);
        }

        template<std::size_t... i>
        static transformed_storage_type transform_storage(shared_ptr<const S> sp, T& t, index_pack<i...> indices)
        {
          return NodeTransformation::transform_storage(sp,t,ChildTransformation<i>::transform_storage(sp->template childStorage<i>(),t)...);
        }

        template<std::size_t... i>
        static transformed_storage_type transform_storage(shared_ptr<const S> sp, const T& t, index_pack<i...> indices)
        {
          return NodeTransformation::transform_storage(sp,t,ChildTransformation<i>::transform_storage(sp->template childStorage<i>(),t)...);
        }

      };


      // the specialization of transformation<> for the CompositeNode. This just extracts the
      // CompositeNode::ChildTypes member and forwards to the helper struct
      template<typename S, typename T>
      struct TransformTree<S,T,VariadicCompositeNodeTag,true>
      {

      private:

        typedef typename S::ChildTypes ChildTypes;

        static typename tuple_index_pack_builder<ChildTypes>::type child_indices()
        {
          return typename tuple_index_pack_builder<ChildTypes>::type();
        }

      public:

        typedef typename transform_variadic_composite_node<S,ChildTypes,T>::transformed_type transformed_type;
        typedef typename transform_variadic_composite_node<S,ChildTypes,T>::transformed_storage_type transformed_storage_type;

        static transformed_type transform(const S& s, T& t)
        {
          return transform_variadic_composite_node<S,ChildTypes,T>::transform(s,t,child_indices());
        }

        static transformed_type transform(const S& s, const T& t)
        {
          return transform_variadic_composite_node<S,ChildTypes,T>::transform(s,t,child_indices());
        }

        static transformed_storage_type transform_storage(shared_ptr<const S> sp, T& t)
        {
          return transform_variadic_composite_node<S,ChildTypes,T>::transform_storage(sp,t,child_indices());
        }

        static transformed_storage_type transform_storage(shared_ptr<const S> sp, const T& t)
        {
          return transform_variadic_composite_node<S,ChildTypes,T>::transform_storage(sp,t,child_indices());
        }

      };

      // non-recursive version of the VariadicCompositeNode transformation.
      template<typename S, typename T>
      struct TransformTree<S,T,VariadicCompositeNodeTag,false>
        : public TransformTreeNonRecursive<S,T>
      {};


#endif // HAVE_VARIADIC_TEMPLATES

      // generic transformation descriptor for empty nodes
      struct EmptyNodeTransformation
      {
        // there is nothing to recurse into here
        static const bool recursive = false;
      };

      // handle empty nodes
      template<typename T, bool recursive>
      struct TransformTree<EmptyNode,T,EmptyNodeTag,recursive>
      {
        // get transformed type from specification
        typedef EmptyNode transformed_type;
        typedef shared_ptr<EmptyNode> transformed_storage_type;

        // delegate instance transformation to per-node specification
        static transformed_type transform(const EmptyNode& s, const T& t)
        {
          DUNE_THROW(NotImplemented,"this should never get called!");
        }

        static transformed_storage_type transform_storage(shared_ptr<const EmptyNode> en, const T& t)
        {
          //return const_pointer_cast<transformed_type>(en); // Dune built-in shared_ptr does not support this!
          return emptyNodePtr();
        }
      };

#endif // DOXYGEN

      //! \} group Traversal

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_TRANSFORMATION_HH
