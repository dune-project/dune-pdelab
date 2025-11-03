#ifndef DUNE_PDELAB_COMMON_PARTITION_ITERATORSPLIT_HH
#define DUNE_PDELAB_COMMON_PARTITION_ITERATORSPLIT_HH

#include <dune/pdelab/common/partition/halo/region.hh>
#include <dune/pdelab/common/partition/halo/colored.hh>
#include <dune/pdelab/common/partition/halo/uncolored.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/concepts/gridview.hh>

#include <dune/common/iteratorrange.hh>

#include <vector>
#include <memory>
#include <array>
#include <algorithm>

namespace Dune::PDELab::EntitySetPartition {

namespace Impl {

/**
 * @brief Partition set with entities evenly split on several patches
 * @details This is a naive partition on an entity set where iterators of the entity
 * set are evenly split and assigned to different patches of one label set.
 *
 * @tparam ES     Entity set to partition
 */
template<Dune::Concept::GridView ES>
class IteratorSplitMixin {
public:
  //! Uderlying entity set
  using EntitySet = ES;
  //! Entity being partitioned
  using Element = typename EntitySet::template Codim<0>::Entity;
  //! Range of entities grouped by a patch
  using PatchSet = IteratorRange<typename EntitySet::template Codim<0>::Iterator>;
  //! Range of patches grouped by a label
  using LabelSet = std::vector<PatchSet>;
  //! Range of labels
  using PartitionSet = std::array<LabelSet,1>;

  /**
   * @brief Construct a IteratorSplit partition
   *
   * @param entity_set  Grid view to operate on
   * @param patches Number of patches to split the entity set (0 means automatic choice)
   */
  explicit IteratorSplitMixin(const EntitySet& entity_set, std::size_t patches)
    : _entity_set{entity_set}
    , _patches{patches}
  {
    update(entity_set);
  }

  //! Uderlying entity set
  [[nodiscard]] auto entitySet() const noexcept { return _entity_set; }

  //! Begin of the partition set
  [[nodiscard]] auto begin() const {
    return _partition_set->begin();
  }

  //! End of the partition set
  [[nodiscard]] auto end() const {
    return _partition_set->end();
  }

protected:
  /**
   * @brief Update partition with a new entity set
   *
   * @param entity_set  Entity set to operate on
   */
  void update(EntitySet entity_set) {
    _entity_set = entity_set;

    if (_patches == 0) {
      std::size_t entities_per_patch = 100;
      switch(EntitySet::dimension) {
        case 1: entities_per_patch = 125; break;
        case 2: entities_per_patch = 250; break;
        default: entities_per_patch = 500; break;
      }
      _patches = std::max<std::size_t>(1, entity_set.size(0) / entities_per_patch);
    }

    PartitionSet partition_set;
    partition_set[0].clear();
    partition_set[0].reserve(_patches);

    auto begin_it = _entity_set.template begin<0>();
    auto end_it = _entity_set.template end<0>();
    auto dist = std::distance(begin_it, end_it);
    auto chunk = dist / _patches;
    auto remainder = dist % _patches;

    for (size_t i = 0; i < _patches-1; ++i) {
      auto next_end = std::next(begin_it, chunk + (remainder ? 1 : 0));
      partition_set[0].emplace_back(begin_it, next_end);

      begin_it = next_end;
      if (remainder) remainder -= 1;
    }

    // last chunk
    partition_set[0].emplace_back(begin_it, end_it);
    _partition_set = std::make_shared<PartitionSet>(std::move(partition_set));
  }

private:
  std::shared_ptr<PartitionSet> _partition_set;
  EntitySet _entity_set;
  std::size_t _patches;
};

} // namespace Impl

/**
 * @brief Partition set with entities evenly split on several patches
 * @details This is a naive partition on an entity set where iterators of the grid
 * view are evenly split and assigned to different patches of one label set.
 * The patches are neighboring each other, thus, inducing a halo overlap in some entities of the partition.
 *
 * @tparam EntitySet     Grid view to partition
 */
template<Dune::Concept::GridView EntitySet>
struct IteratorSplit
  : public Impl::IteratorSplitMixin<EntitySet>
  , public Impl::UncoloredHaloMixin<EntitySet>
{

  /**
   * @brief Construct a IteratorSplit uncolored entity set partition
   *
   * @param entity_set      The entity set to operate on
   * @param patches         Number of patches to have in total (0 means automatic choice)
   * @param halo_distance   Neighbor distance where another entity in the same label set is considered connected
   */
  explicit IteratorSplit(const EntitySet& entity_set, std::size_t patches = 0, std::size_t halo_distance = all_interior_halo_region)
   : Impl::IteratorSplitMixin<EntitySet>{entity_set, patches}
   , Impl::UncoloredHaloMixin<EntitySet>::UncoloredHaloMixin{halo_distance}
  {
    Impl::UncoloredHaloMixin<EntitySet>::updateHalo(*this);
  }

  //! Update the partition with a new entity set
  void update(const EntitySet& entity_set) {
    Impl::IteratorSplitMixin<EntitySet>::update(entity_set);
    Impl::UncoloredHaloMixin<EntitySet>::updateHalo(*this);
  }
};

/**
 * @brief Colored partition set with entities evenly split on several patches
 * @details This is a naive partition on an entity set where iterators of the grid
 * view are evenly split and assigned to different patches.
 * Labels are colored so the halo region of every entity in the set is interior.
 *
 * @tparam EntitySet     Entity set to partition
 */
template<Dune::Concept::GridView EntitySet>
struct IteratorSplitColored
  : public Impl::ColoredHaloAdaptor<Impl::IteratorSplitMixin<EntitySet>>
{
  /**
   * @brief Construct a IteratorSplit colored entity set partition
   *
   * @param entity_set  The entity set to operate on
   * @param patches     Number of patches to have in total
   * @param halo        Distance another entity in the same label set is considered connected
   */
  explicit IteratorSplitColored(const EntitySet& entity_set, std::size_t patches = 0, std::size_t halo_distance = all_interior_halo_region)
    : Impl::ColoredHaloAdaptor<Impl::IteratorSplitMixin<EntitySet>>{Impl::IteratorSplitMixin<EntitySet>{entity_set, patches}, halo_distance}
    , _patches{patches}
  {}

  //! Update the partition with a new entity set
  void update(const EntitySet& entity_set) {
    Impl::IteratorSplitMixin<EntitySet> base(entity_set, _patches);
    Impl::ColoredHaloAdaptor<Impl::IteratorSplitMixin<EntitySet>>::update(std::move(base));
  }

private:
  std::size_t _patches;
};

} // namespace Dune::PDELab::EntitySetPartition

#endif // DUNE_PDELAB_COMMON_PARTITION_ITERATORSPLIT_HH
