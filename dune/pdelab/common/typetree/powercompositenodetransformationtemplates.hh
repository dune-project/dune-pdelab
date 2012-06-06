// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_POWERCOMPOSITENODETRANSFORMATIONTEMPLATES_HH
#define DUNE_PDELAB_COMMON_TYPETREE_POWERCOMPOSITENODETRANSFORMATIONTEMPLATES_HH

#include <cstddef>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Transformation
       *  \ingroup TypeTree
       *  \{
       */

      template<typename SourceNode, typename Transformation, template<typename,typename,std::size_t> class TransformedNode>
      struct GenericPowerNodeTransformationTemplate
      {
        template<typename TC>
        struct result
        {
          typedef TransformedNode<SourceNode,TC,SourceNode::CHILDREN> type;
        };
      };


#if HAVE_VARIADIC_TEMPLATES

      template<typename SourceNode, typename Transformation, template<typename,typename...> class TransformedNode>
      struct GenericVariadicCompositeNodeTransformationTemplate
      {
        template<typename... TC>
        struct result
        {
          typedef TransformedNode<SourceNode,TC...> type;
        };
      };

#endif // HAVE_VARIADIC_TEMPLATES


      template<typename SourceNode,
               typename Transformation,
               template<typename SourceNode_,
                        typename C0,
                        typename C1,
                        typename C2,
                        typename C3,
                        typename C4,
                        typename C5,
                        typename C6,
                        typename C7,
                        typename C8,
                        typename C9> class TransformedNode>
      struct GenericCompositeNodeTransformationTemplate
      {
        template<typename TC0,
                 typename TC1,
                 typename TC2,
                 typename TC3,
                 typename TC4,
                 typename TC5,
                 typename TC6,
                 typename TC7,
                 typename TC8,
                 typename TC9>
        struct result
        {
          typedef TransformedNode<SourceNode,TC0,TC1,TC2,TC3,TC4,TC5,TC6,TC7,TC8,TC9> type;
        };
      };


      //! \} group Transformation

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_POWERCOMPOSITENODETRANSFORMATIONTEMPLATES_HH
