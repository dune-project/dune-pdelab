// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_NODEINTERFACE_HH
#define DUNE_PDELAB_COMMON_TYPETREE_NODEINTERFACE_HH

#include <cstddef>

#include <dune/common/documentation.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Nodes
       *  \ingroup TypeTree
       *  \{
       */

      /** \brief Interface for nodes in a \ref TypeTree.
       *
       * This class cannot be used itself, it is for documentation purposes
       * only.
       *
       * \note Constructor signatures are explicitly not specified by this
       *       interface.
       * \note In addition, every node in a tree must be derived from one of
       *       the node base classes LeafNode, PowerNode, CompositeNode, or
       *       VariadicCompositeNode, or from a base class for a
       *       yet-to-be-defined new node type.
       */
      struct NodeInterface
      {
        //! Whether this is a leaf node in a \ref TypeTree.
        static const bool isLeaf = implementationDefined;

        //! Whether this is a power node in the \ref TypeTree.
        static const bool isPower = implementationDefined;

        //! Whether this is a composite node in the \ref TypeTree.
        static const bool isComposite = implementationDefined;

        //! Number of children of this node in the \ref TypeTree
        static const std::size_t CHILDREN = implementationDefined;

        //! The type tag that describes what kind of node this is
        /**
         * One of LeafNodeTag, PowerNodeTag, CompositeNodeTag, or
         * VariadicCompositeNodeTag.  Other tags are also possible when new
         * kinds of nodes are defined.
         */
        typedef ImplementationDefined NodeTag;

        //! container type to pass around a collection of children
        /**
         * \note This typedef is not present for leaf nodes.
         */
        typedef ImplementationDefined NodeStorage;
      };

      //! \} group Nodes

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif //  DUNE_PDELAB_COMMON_TYPETREE_NODEINTERFACE_HH
