// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_TAGS_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_TAGS_HH

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    struct PowerGridFunctionSpaceTag {};

    struct VectorGridFunctionSpaceTag : public PowerGridFunctionSpaceTag {};

    struct CompositeGridFunctionSpaceTag {};

    struct LeafGridFunctionSpaceTag {};

    //! Tag for the intermediate base class of the CompositeGridFunctionSpace.
    struct CompositeGridFunctionSpaceBaseTag {};

    struct EntityBlockedOrderingTag {};

    //! \brief Indicate lexicographic ordering of the unknowns of non-leaf
    //!        grid function spaces.
    /**
     * This class instructs the non-leaf GridFunctionSpaces to order the dofs
     * of the child-GridFunctionSpaces in a lexicographic manner in the
     * combined dof-vector, i.e. first all dofs of child 0, then all dofs of
     * child 1 and so on.
     */
    struct LexicographicOrderingTag { };


    //! Tag denoting a PowerLocalFunctionSpace
    struct PowerLocalFunctionSpaceTag {};

    //! Tag denoting a CompositeLocalFunctionSpace
    struct CompositeLocalFunctionSpaceTag {};

    //! Tag denoting a LeafLocalFunctionSpace
    struct LeafLocalFunctionSpaceTag {};


    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_TAGS_HH
