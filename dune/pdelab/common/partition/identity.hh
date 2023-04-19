#ifndef DUNE_PDELAB_COMMON_PARTITION_IDENTITY_HH
#define DUNE_PDELAB_COMMON_PARTITION_IDENTITY_HH

#include <dune/assembler/common/partition/overlap.hh>

#include <dune/grid/concepts/gridview.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/partitionset.hh>

#include <dune/common/iteratorrange.hh>

#include <type_traits>
#include <array>
#include <memory>

namespace Dune::Assembler::EntitySetPartition {

/**
 * @brief Entity set partition with all entities assigned to one patch
 * @details This is the simples partition on an grid view
 * This class models the EntitySetPartition concept
 *
 * @tparam GV     Grid view to partition
 * @tparam codim  Codimension to partition
 * @tparam pit    Parallel partition to operate on
 */
template<Dune::Concept::GridView GV, int codim = 0, PartitionIteratorType pit = PartitionIteratorType::All_Partition>
class Identity : public Impl::ColoredOverlapMixin {
public:
  //! Uderlying grid view
  using GridView = GV;
  //! Entity being partitioned
  using Entity = typename GridView::Codim<codim>::Entity;
  //! Range of entities grouped by a patch
  using PatchSet = IteratorRange<typename GridView::Codim<codim>::Iterator>;
  //! Range of patches grouped by a label
  using LabelSet = std::array<PatchSet,1>;
  //! Range of labels
  using PartitionSet = std::array<LabelSet,1>;

  //! Construct Identity partion from a grid view
  explicit Identity(const GridView& grid_view)
    : _grid_view{grid_view}
    , _partition_set{LabelSet{PatchSet{_grid_view.template begin<codim,pit>(), _grid_view.template end<codim,pit>()}}}
  {}

  //! Uderlying grid view
  [[nodiscard]] GridView gridView() const noexcept { return _grid_view; }

  //! Range of the partition set
  [[nodiscard]] const PartitionSet& range() const noexcept { return _partition_set; }

private:
  GridView _grid_view;
  PartitionSet _partition_set;
};

} // namespace Dune::Assembler::EntitySetPartition

#endif // DUNE_PDELAB_COMMON_PARTITION_IDENTITY_HH
