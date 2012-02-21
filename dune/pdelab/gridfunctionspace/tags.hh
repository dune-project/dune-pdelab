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

    /** \brief Tag indicating an arbitrary number of unkowns per entity.
     *
     * class used to pass compile-time parameter to the GridFunctionSpace.
     */
    struct GridFunctionGeneralMapper {};

    //! Tag indicating a function space with a single unknown attached to every
    //! entity of a exactly one single codimension.
    struct SingleCodimMapper {};

    //! Tag denoting a PowerLocalFunctionSpace
    struct PowerLocalFunctionSpaceTag {};

    //! Tag denoting a CompositeLocalFunctionSpace
    struct CompositeLocalFunctionSpaceTag {};

    //! Tag denoting a LeafLocalFunctionSpace
    struct LeafLocalFunctionSpaceTag {};

    //! Tag for denoting possibly nested containers, requiring a recursive
    //! allocation algorithm.
    struct HierarchicContainerAllocationTag {};

    //! Tag for denoting that a backend / ordering will always spawn flat
    //! containers.
    struct FlatContainerAllocationTag {};

    //! Tag denoting that an ordering will work with the default implementation
    //! of the LFSIndexCache.
    struct DefaultLFSCacheTag {};

    //! Tag denoting that an ordering will work with the simplified version of
    //! the LFSIndexCache.
    struct SimpleLFSCacheTag {};

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_TAGS_HH
