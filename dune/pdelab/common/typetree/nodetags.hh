// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_NODETAGS_HH
#define DUNE_PDELAB_COMMON_TYPETREE_NODETAGS_HH

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup Nodes
       *  \ingroup TypeTree
       *  \{
       */

      //! Tag designating a leaf node.
      struct LeafNodeTag {};

      //! Tag designating a power node.
      struct PowerNodeTag {};

      //! Tag designating a composite node.
      struct CompositeNodeTag {};

      //! Tag designating a composite node that is based on variadic templates.
      struct VariadicCompositeNodeTag {};

      //! Tag designating an empty child node in a non-variadic composite node.
      struct EmptyNodeTag {};

      //! This is a marker to fill empty spaces in composite nodes if support for variadic templates is missing.
      struct EmptyNode
      {
        typedef EmptyNodeTag NodeTag;
      };

#ifndef DOXYGEN

      //! Special tag used as start value in algorithms.
      struct StartTag {};



#endif // DOXYGEN

      //! \} group Nodes

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_NODETAGS_HH
