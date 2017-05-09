// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_TAGS_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_TAGS_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/typetree/utility.hh>

#include <dune/pdelab/common/dofindex.hh>
#include <dune/pdelab/common/simpledofindex.hh>

#include <numeric>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    struct FunctionSpaceTag {};

    struct GridFunctionSpaceTag : public FunctionSpaceTag {};

    struct PowerGridFunctionSpaceTag : public GridFunctionSpaceTag {};

    struct VectorGridFunctionSpaceTag : public PowerGridFunctionSpaceTag {};

    struct CompositeGridFunctionSpaceTag : public GridFunctionSpaceTag {};

    struct LeafGridFunctionSpaceTag : public GridFunctionSpaceTag {};

    template<typename ProxiedGFSTag>
    struct GridFunctionSubSpaceTag
      : public ProxiedGFSTag
    {};

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

    //! \brief Indicate interleaved ordering of the unknowns of non-leaf
    //!        grid function spaces according to a given blocking pattern.
    /**
     * This class instructs the non-leaf GridFunctionSpaces to order the dofs
     * of the child-GridFunctionSpaces in an interleaved manner in the
     * combined dof-vector. The sizes of the individual blocks have to be passed
     * to the constructor of the tag.
     *
     * \note In the vast majority of scenarios, you will want to use the
     *       EntityBlockedOrderingTag instead of this one, as it is much less error-prone
     *       and works in a wider variety of settings. Only use the InterleavedOrderingTag
     *       if you know that the EntityBlockedOrderingTag will not work for you!
     */
    struct InterleavedOrderingTag
    {

      //! Constructs an InterleavedOrderingTag with a block structure given by the initializer list sizes.
      /**
       * If you have a sufficiently recent compiler, this constructor enables a much more readable
       * syntax when creating a GridFunctionSpace. Assuming that GFS is a PowerGridFunctionSpace and VBE its
       * associated vector backend, you can shorten the verbose
       *
       * \code
       * std::vector<std::size_t> sizes(3);
       * sizes[0] = 2;
       * sizes[1] = 5;
       * sizes[2] = 3;
       * GFS gfs(child_gfs,VBE(),InterleavedOrderingTag(sizes));
       * \endcode
       *
       * to the much shorter and more readable
       *
       * \code
       * GFS gfs(child_gfs,VBE(),{2,5,3});
       * \endcode
       */
      InterleavedOrderingTag(std::initializer_list<std::size_t> sizes)
        : _offsets(sizes.size() + 1,0)
      {
        std::partial_sum(sizes.begin(),sizes.end(),_offsets.begin() + 1);
      }

      //! Constructs an InterleavedOrderingTag with a block structure given by the std::vector sizes.
      InterleavedOrderingTag(std::vector<std::size_t> sizes)
        : _offsets(sizes.size() + 1,0)
      {
        std::partial_sum(sizes.begin(),sizes.end(),_offsets.begin() + 1);
      }

      //! Returns a list of offsets for the child blocks.
      const std::vector<std::size_t>& offsets() const
      {
        return _offsets;
      }

    private:

      std::vector<std::size_t> _offsets;
    };

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

    struct EmptyParams
    {};

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

    using DefaultLeafOrderingTag = LeafOrderingTag<EmptyParams>;

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

    namespace Experimental {

      struct DuneFunctionsCacheTag {};

    }

    template<typename GFS, typename Tag>
    struct _build_dof_index_type
    {
      typedef Dune::PDELab::DOFIndex<std::size_t,TypeTree::TreeInfo<GFS>::depth,2> type;
    };

    template<typename GFS>
    struct _build_dof_index_type<GFS,SingleCodimMapper>
    {
      typedef SimpleDOFIndex<typename GFS::Traits::SizeType> type;
    };


    template<typename GFS>
    struct build_dof_index_type
    {
      typedef typename _build_dof_index_type<GFS,typename GFS::OrderingTag>::type type;
    };


#ifndef DOXYGEN

    //! GridFunctionSpace to LocalFunctionSpace transformation.
    /**
     * gfs_to_lfs describes the transformation of a GridFunctionSpace tree to its corresponding
     * LocalFunctionSpace tree and holds any information that may be required for performing
     * the transformation.
     *
     * \warning The exact meaning of the template parameter is an implementation detail
     *          and may change at any time, as the only class that is supposed to instantiate
     *          the transformation is LocalFunctionSpace. Implementors of custom transformation
     *          descriptors should only use information exported by gfs_to_lfs. In particular,
     *          the registration declaration should not make any assumptions on GFS and just
     *          treat it as some kind of opaque parameter type.
     *
     * \tparam GFS  the root GridFunctionSpace that the resulting LocalFunctionSpace tree
     *              will be based on.
     */
    template<typename GFS>
    struct gfs_to_lfs {

      //! The MultiIndex type that will be used in the resulting LocalFunctionSpace tree.
      //typedef Dune::PDELab::MultiIndex<std::size_t,TypeTree::TreeInfo<GFS>::depth> MultiIndex;
      typedef typename build_dof_index_type<GFS>::type DOFIndex;

    };

#endif // DOXYGEN

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_TAGS_HH
