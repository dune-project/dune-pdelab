// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_TRANSFORMATIONUTILITIES_HH
#define DUNE_PDELAB_COMMON_TYPETREE_TRANSFORMATIONUTILITIES_HH

#include <dune/pdelab/common/typetree/nodetags.hh>
#include <dune/common/exceptions.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup TypeTree
       *  \ingroup PDELab
       *  \{
       */

      template<typename SourceNode, typename Transformation, typename TransformedNode>
      struct WrappingLeafNodeTransformation
      {

        typedef TransformedNode transformed_type;
        typedef shared_ptr<transformed_type> transformed_storage_type;

        static transformed_type transform(const SourceNode& s, const Transformation& t)
        {
          return transformed_type(s);
        }

        static transformed_storage_type transform_storage(shared_ptr<const SourceNode> s, const Transformation& t)
        {
          return make_shared<transformed_type>(s);
        }

      };


      template<typename SourceNode, typename Transformation, template<typename Child> class TransformedNodeTemplate>
      struct TemplatizedWrappingPowerNodeTransformation
      {

        template<typename TC>
        struct result
        {
          typedef typename TransformedNodeTemplate<TC>::type type;
          typedef shared_ptr<type> storage_type;
        };

        template<typename TC>
        static typename result<TC>::type transform(const SourceNode& s, const Transformation& t, const array<shared_ptr<TC>,result<TC>::type::CHILDREN>& children)
        {
          return typename result<TC>::type(s,children);
        }

        template<typename TC>
        static typename result<TC>::storage_type transform_storage(shared_ptr<const SourceNode> s, const Transformation& t, const array<shared_ptr<TC>,result<TC>::type::CHILDREN>& children)
        {
          return make_shared<typename result<TC>::type>(s,children);
        }

      };


      template<typename SourceNode, template<typename,typename,std::size_t> class TransformedNode>
      struct WrappingPowerNodeTransformationTemplate
      {
        template<typename TC>
        struct result
        {
          typedef TransformedNode<SourceNode,TC,SourceNode::CHILDREN> type;
        };
      };

      template<typename SourceNode, typename Transformation, template<typename,typename,std::size_t> class TransformedNode>
      struct WrappingPowerNodeTransformation
        : public TemplatizedWrappingPowerNodeTransformation<SourceNode,
                                                            Transformation,
                                                            WrappingPowerNodeTransformationTemplate<SourceNode,TransformedNode>::template result
                                                            >
      {};


      template<typename SourceNode, typename Transformation, template<typename...> class TransformedNodeTemplate>
      struct TemplatizedWrappingVariadicCompositeNodeTransformation
      {

        template<typename... TC>
        struct result
        {
          typedef typename TransformedNodeTemplate<TC...>::type type;
          typedef shared_ptr<type> storage_type;
        };

        template<typename... TC>
        static typename result<TC...>::type transform(const SourceNode& s, const Transformation& t, shared_ptr<TC>... children)
        {
          return typename result<TC...>::type(s,children...);
        }

        template<typename... TC>
        static typename result<TC...>::storage_type transform_storage(shared_ptr<const SourceNode> s, const Transformation& t, shared_ptr<TC>... children)
        {
          return make_shared<typename result<TC...>::type>(s,children...);
        }

      };

      template<typename SourceNode, template<typename...> class TransformedNode>
      struct WrappingVariadicCompositeNodeTransformationTemplate
      {
        template<typename... TC>
        struct result
        {
          typedef TransformedNode<SourceNode,TC...> type;
        };
      };

      template<typename SourceNode, typename Transformation, template<typename...> class TransformedNode>
      struct WrappingVariadicCompositeNodeTransformation
        : public TemplatizedWrappingVariadicCompositeNodeTransformation<SourceNode,
                                                                        Transformation,
                                                                        WrappingVariadicCompositeNodeTransformationTemplate<SourceNode,TransformedNode>::template result
                                                                        >
      {};

      //! \} group TypeTree

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_TRANSFORMATIONUTILITIES_HH
