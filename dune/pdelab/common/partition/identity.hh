#ifndef DUNE_PDELAB_COMMON_PARTITION_IDENTITY_HH
#define DUNE_PDELAB_COMMON_PARTITION_IDENTITY_HH

#include <dune/pdelab/common/partition/halo/region.hh>

#include <dune/grid/concepts/gridview.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/partitionset.hh>

#include <dune/common/iteratorrange.hh>

#include <type_traits>
#include <array>
#include <memory>

namespace Dune::PDELab::EntitySetPartition {

/**
 * @brief Entity set partition with all entities assigned to one patch (i.e. no partitioning)
 * @details This is the simplest partition on an entity set, but with the interface of a EntitySetPartition.
 *
 * @tparam ES     Entity set being partitioned
 *
 */
template<Dune::Concept::GridView ES>
class Identity {
public:
  //! Uderlying entity set
  using EntitySet = ES;
  //! Entity being partitioned
  using Element = typename EntitySet::template Codim<0>::Entity;
  //! Range of entities grouped by a patch
  using PatchSet = IteratorRange<typename EntitySet::template Codim<0>::Iterator>;
  //! Range of patches grouped by a label
  using LabelSet = std::array<PatchSet,1>;
  //! Range of labels
  using PartitionSet = std::array<LabelSet,1>;

  //! Construct Identity partition from a entity set
  explicit Identity(const EntitySet& entity_set)
    : _entity_set{entity_set}
  {
    update(entity_set);
  }

  //! Halo region of the partition for an entity
  [[nodiscard]] constexpr static auto haloRegion(const Dune::Concept::Entity auto& entity) {
    return interior_halo_region;
  }

  //! Underlying entity set
  [[nodiscard]] EntitySet entitySet() const noexcept { return _entity_set; }

  //! begin of the partition set range
  [[nodiscard]] auto begin() const { return _partition_set.begin(); }

  //! end of the partition set range
  [[nodiscard]] auto end() const { return _partition_set.end(); }

  //! Update the partition with a new entity set
  void update(EntitySet entity_set) {
    _entity_set = entity_set;
    _partition_set = PartitionSet{LabelSet{PatchSet{entity_set.template begin<0>(), entity_set.template end<0>()}}};
  }

private:
  EntitySet _entity_set;
  PartitionSet _partition_set;
};

} // namespace Dune::PDELab::EntitySetPartition

#endif // DUNE_PDELAB_COMMON_PARTITION_IDENTITY_HH
