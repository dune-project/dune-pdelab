#ifndef DUNE_PDELAB_COMMON_PARTITION_SIMPLE_HH
#define DUNE_PDELAB_COMMON_PARTITION_SIMPLE_HH

#include <dune/pdelab/common/partition/halo.hh>
#include <dune/pdelab/common/partition/coloring.hh>
// #include <dune/pdelab/common/trace.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/concepts/gridview.hh>

#include <dune/common/iteratorrange.hh>

#include <vector>
#include <memory>
#include <array>

namespace Dune::PDELab::inline Experimental::EntitySetPartition {

namespace Impl {

/**
 * @brief Partition set with entities evenly split on several patches
 * @details This is a naive partition on an grid view where iterators of the grid
 * view are evenly split and assigned to different patches of one label set.
 *
 * @tparam GV     Grid view to partition
 * @tparam codim  Codimension to partition
 * @tparam pit    Parallel partition to operate on
 */
template<Dune::Concept::GridView GV, int codim = 0, PartitionIteratorType pit = PartitionIteratorType::All_Partition>
class SimpleMixin {
public:
  //! Uderlying grid view
  using GridView = GV;
  //! Entity being partitioned
  using Entity = typename GridView::Codim<codim>::Entity;
  //! Range of entities grouped by a patch
  using PatchSet = IteratorRange<typename GridView::Codim<codim>::Iterator>;
  //! Range of patches grouped by a label
  using LabelSet = std::vector<PatchSet>;
  //! Range of labels
  using PartitionSet = std::array<LabelSet,1>;

  /**
   * @brief Construct a Simple partition
   *
   * @param grid_view  Grid view to operate on
   * @param patches Number of patches to split the entity set
   */
  explicit SimpleMixin(const GridView& grid_view, std::size_t patches)
   : _partition_set{std::make_shared<PartitionSet>()}
   , _grid_view{grid_view}
  {
    update(patches);
  }

  //! Uderlying grid view
  [[nodiscard]] auto entitySet() const noexcept { return _grid_view; }

    //! Range of the partition set
  [[nodiscard]] const PartitionSet& range() const noexcept { return *_partition_set; }

protected:
  /**
   * @brief Sets a new set of patches
   *
   * @param patches  Number of patches in the label set
   */
  void update(std::size_t patches) {
    TRACE_EVENT("dune", "EntitySet::updatePartition");
    if (patches == 0)
      DUNE_THROW(InvalidStateException, "There cannot be zero partitions");
    (*_partition_set)[0].clear();
    (*_partition_set)[0].reserve(patches);

    auto begin_it = _grid_view.template begin<codim,pit>();
    auto end_it = _grid_view.template end<codim,pit>();
    auto dist = std::distance(begin_it, end_it);
    auto chunk = dist / patches;
    auto remainder = dist % patches;

    for (size_t i = 0; i < patches-1; ++i) {
      auto next_end = std::next(begin_it, chunk + (remainder ? 1 : 0));
      (*_partition_set)[0].emplace_back(begin_it, next_end);

      begin_it = next_end;
      if (remainder) remainder -= 1;
    }

    // last chunk
    (*_partition_set)[0].emplace_back(begin_it, end_it);
  }

private:
  std::shared_ptr<PartitionSet> _partition_set;
  GridView _grid_view;
};

} // namespace Impl

/**
 * @brief Partition set with entities evenly split on several patches
 * @details This is a naive partition on an grid view where iterators of the grid
 * view are evenly split and assigned to different patches of one label set.
 * The patches may be neigboring each other, thus, inducing a shared memory region on
 * several entities of the set.
 * This class models the EntitySetPartition concept.
 *
 * @tparam GV     Grid view to partition
 * @tparam codim  Codimension to partition
 * @tparam pit    Parallel partition to operate on
 */
template<Dune::Concept::GridView GV, int codim = 0, PartitionIteratorType pit = PartitionIteratorType::All_Partition>
struct Simple
  : public Impl::SimpleMixin<GV, codim, pit>
  , public Impl::UncoloredOverlapMixin<GV>
{

  /**
   * @brief Construct a Simple uncolored entiy set partition
   *
   * @param grid_view   The grid view to operate on
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  explicit Simple(const GV& grid_view, std::size_t patches, std::size_t halo)
   : Impl::SimpleMixin<GV, codim, pit>{grid_view, patches}
   , Impl::UncoloredOverlapMixin<GV>::UncoloredOverlapMixin{grid_view}
  {
    Impl::UncoloredOverlapMixin<GV>::updatePartitionOverlap(*this, halo);
  }

  /**
   * @brief Update partition with new patches
   *
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  void update(std::size_t partitions, std::size_t halo) {
    Impl::SimpleMixin<GV, codim, pit>::update(partitions);
    Impl::UncoloredOverlapMixin<GV>::updatePartitionOverlap(*this, halo);
  }
};

/**
 * @brief Colored partition set with entities evenly split on several patches
 * @details This is a naive partition on an grid view where iterators of the grid
 * view are evenly split and assigned to different patches.
 * Labels are colored so the memory region of every entity in the set is private.
 * This class models the EntitySetPartition concept.
 *
 * @tparam GV     Grid view to partition
 * @tparam codim  Codimension to partition
 * @tparam pit    Parallel partition to operate on
 */
template<Dune::Concept::GridView GridView, int codim = 0, PartitionIteratorType pit = PartitionIteratorType::All_Partition>
struct SimpleColored
  : public Impl::ColoringAdaptor<Impl::SimpleMixin<GridView, codim, pit>>
{
  /**
   * @brief Construct a Simple colored entiy set partition
   *
   * @param grid_view   The grid view to operate on
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  explicit SimpleColored(const GridView& grid_view, std::size_t partitions, std::size_t halo)
   : Impl::ColoringAdaptor<Impl::SimpleMixin<GridView, codim, pit>>{Impl::SimpleMixin<GridView, codim, pit>{grid_view, partitions}, halo}
  {}

  /**
   * @brief Update partition with new patches and labels
   *
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  void update(std::size_t partitions, std::size_t halo) {
    Impl::ColoringAdaptor<Impl::SimpleMixin<GridView, codim, pit>>::update(Impl::SimpleMixin<GridView, codim, pit>{this->entitySet(), partitions}, halo);
  }
};

} // namespace Dune::PDELab::EntitySetPartition

#endif // DUNE_PDELAB_COMMON_PARTITION_SIMPLE_HH
