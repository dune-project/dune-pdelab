#ifndef DUNE_PDELAB_COMMON_PARTITION_IDENTITY_HH
#define DUNE_PDELAB_COMMON_PARTITION_IDENTITY_HH

#include <dune/pdelab/common/partition/halo.hh>

#include <dune/grid/concepts/gridview.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/partitionset.hh>

#include <dune/common/iteratorrange.hh>

#include <type_traits>
#include <array>
#include <memory>

namespace Dune::PDELab::inline Experimental::EntitySetPartitioner {

/**
 * @brief Entity set partition with all entities assigned to one patch
 * @details This is the simples partition on an grid view
 * This class models the EntitySetPartition concept
 *
 * @tparam ES     Grid view to partition
 * @tparam codim  Codimension to partition
 * @tparam pit    Parallel partition to operate on
 */
template<Dune::Concept::GridView ES, int codim = 0, PartitionIteratorType pit = PartitionIteratorType::All_Partition>
class Identity : public Impl::ColoredHaloMixin {
public:
  //! Uderlying grid view
  using EntitySet = ES;
  //! Entity being partitioned
  using Entity = typename EntitySet::template Codim<codim>::Entity;
  //! Range of entities grouped by a patch
  using PatchSet = IteratorRange<typename EntitySet::template Codim<codim>::Iterator>;
  //! Range of patches grouped by a label
  using LabelSet = std::array<PatchSet,1>;
  //! Range of labels
  using PartitionSet = std::array<LabelSet,1>;

  constexpr static auto isRegionAlwaysUnique() { return std::true_type(); }

  [[nodiscard]] static auto maxLabels() {
    return std::integral_constant<size_t, 1>{};
  }

  [[nodiscard]] static auto maxPatches() {
    return std::integral_constant<size_t, 1>{};
  }

  //! Construct Identity partion from a grid view
  explicit Identity(const EntitySet& entity_set)
    : _entity_set{entity_set}
    , _partition_set{LabelSet{PatchSet{_entity_set.template begin<codim,pit>(), _entity_set.template end<codim,pit>()}}}
  {}

  //! Uderlying grid view
  [[nodiscard]] EntitySet entitySet() const noexcept { return _entity_set; }

  //! Range of the partition set
  [[nodiscard]] const PartitionSet& range() const noexcept { return _partition_set; }

private:
  EntitySet _entity_set;
  PartitionSet _partition_set;
};

} // namespace Dune::PDELab::inline Experimental::EntitySetPartitioner

#endif // DUNE_PDELAB_COMMON_PARTITION_IDENTITY_HH
