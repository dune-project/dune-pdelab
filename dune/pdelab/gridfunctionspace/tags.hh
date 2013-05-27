// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_TAGS_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_TAGS_HH

#include <dune/grid/common/gridenums.hh>

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

    //! \brief Indicate blocking of the unknowns by grid entity.
    /**
     * This class instructs the non-leaf GridFunctionSpaces to block the dofs
     * of the child-GridFunctionSpaces by grid entity, i.e. first all dofs of
     * all children that belong to vertex 0, then all dofs associated with vertex
     * 1 etc.
     *
     * The EntityBlockedOrdering correctly handles different block sizes for different
     * GeometryTypes as well as GridFunctionSpaces with variable sizes, e.g. for
     * $p$-adaptivity and in a MultiDomain context.
     */
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

    //! Mixin indicating whether a leaf GridFunctionSpace should never assume a const ordering size.
    template<bool v>
    struct NoConstOrderingSize
    {
      static const bool no_const_ordering_size = v;
    };

    namespace {

      // Use this compile-time bridge to shup up GCC warnings

      template<int i>
      struct shift_if_nonnegative
      {
        static const unsigned int value = 1 << i;
      };

      template<>
      struct shift_if_nonnegative<-1>
      {
        static const unsigned int value = 0;
      };

    }

    //! Helper for building the bitmask describing the grid partitions contained in a GFS.
    /**
     * This struct should be used to construct the bitmask for use by the
     * PartitionInfoProvider.
     */
    template<int p0 = -1, int p1 = -1, int p2 = -1, int p3 = -1, int p4 = -1>
    struct PartitionSelector
    {

      static const unsigned int partition_mask =
        shift_if_nonnegative<p0>::value |
        shift_if_nonnegative<p1>::value |
        shift_if_nonnegative<p2>::value |
        shift_if_nonnegative<p3>::value |
        shift_if_nonnegative<p4>::value;

    };

    //! Tag indicating a standard ordering for a leaf GridfunctionSpace.
    /**
     * Any additional policies regarding the ordering should be passed via
     * the template parameter Params. By itself, this tag indicates that the
     * user wants to use the standard, MultiIndex-based ordering infrastructure
     * for this GridFunctionSpace.
     *
     * \tparam Params  Parameter struct for passing additional static information
     *                 to the ordering. This parameter will become the base class
     *                 of the tag.
     */
    template<typename Params>
    struct LeafOrderingTag
      : public Params
    {};

#ifndef DOXYGEN

    typedef PartitionSelector<
      InteriorEntity,
      BorderEntity,
      OverlapEntity,
      FrontEntity,
      GhostEntity
      > AllPartitionSelector;

    typedef PartitionSelector<
      InteriorEntity,
      BorderEntity
      > NonOverlappingPartitionSelector;

#endif // DOXYGEN

    //! Leaf ordering parameters for standard function spaces.
    struct DefaultLeafOrderingParams
      : public NoConstOrderingSize<false>
      , public AllPartitionSelector
    {};

    //! Leaf ordering parameters for non-overlapping function spaces.
    struct NonOverlappingLeafOrderingParams
      : public NoConstOrderingSize<true>
      , public NonOverlappingPartitionSelector
    {};

    //! Default ordering tag for a MultiIndex-based ordering with standard behavior.
    typedef LeafOrderingTag<
      DefaultLeafOrderingParams
      > DefaultLeafOrderingTag;

    //! GridFunctionGeneralMapper is deprecated, use DefaultLeafOrderingTag instead.
    /**
     * \deprecated  Use DefaultLeafOrdering instead.
     */
    typedef DefaultLeafOrderingTag GridFunctionGeneralMapper;

    //! Ordering tag for a MultiIndex-based ordering on nonoverlapping grids with standard behavior.
    typedef LeafOrderingTag<
      NonOverlappingLeafOrderingParams
      > NonOverlappingLeafOrderingTag;

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
