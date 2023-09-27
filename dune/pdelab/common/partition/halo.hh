#ifndef DUNE_PDELAB_COMMON_PARTITION_HALO_HH
#define DUNE_PDELAB_COMMON_PARTITION_HALO_HH


#include <dune/pdelab/concepts/entityset_partition.hh>
#include <dune/pdelab/common/trace.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/concepts/gridview.hh>

#include <dune/common/hybridutilities.hh>

#include <vector>
#include <memory>
#include <atomic>
#include <mutex>
#include <execution>

// splits a grid view into several entity sets grouped by color
namespace Dune::PDELab::inline Experimental::EntitySetPartitioner::Impl {

/**
 * @brief Colored halo mixin
 * @details Colored partitions have no entities in the halo
 * @see EntitySetPartition
 */
struct ColoredHaloMixin {

  /**
   * @brief Whether an entity is in the halo
   * @param entity   The entity to test
   * @return std::false_type
   */
  constexpr static Region region(auto&& entity) { return unique_region; }

  constexpr static auto isRegionAlwaysUnique() { return std::true_type(); }

  [[nodiscard]] static auto maxLabels() {
    return std::numeric_limits<std::size_t>::max();
  }

  [[nodiscard]] static auto maxPatches() {
    return std::numeric_limits<std::size_t>::max();
  }

  [[nodiscard]] static auto maxEntities() {
    return std::numeric_limits<std::size_t>::max();
  }
};

/**
 * @brief Uncolored halo mixin
 * @details Calculates and stores the overlap of the patches wrt halo
 * @see EntitySetPartitioner
 */
template<Dune::Concept::GridView GV>
class UncoloredHaloMixin {
public:

  explicit UncoloredHaloMixin(const GV& entity_set)
    : _mapper{std::make_shared<Dune::MultipleCodimMultipleGeomTypeMapper<GV>>(entity_set, [](GeometryType gt, int dimgrid) { return true; })}
    , _entity_halo{std::make_shared<std::vector<bool>>()}
  {}

  /**
   * @brief Whether an entity is in the halo
   * @param entity   The entity to test
   * @return true if entity is in the halo, false otherwise.
   */
  template<Dune::Concept::Entity Entity>
  [[nodiscard]] Region region(const Entity& entity) const {
    assert(_entity_halo);
    assert((Entity::codimension == 0 or Entity::dimension == 0) && "NotImplemented");
    return (*_entity_halo)[_mapper->index(entity)] ? EntitySetPartitioner::Region::Shared : EntitySetPartitioner::Region::Unique;
  }

  constexpr static auto isRegionAlwaysUnique() { return std::false_type(); }

  [[nodiscard]] static auto maxLabels() {
    return std::numeric_limits<typename GV::IndexSet::IndexType>::max();
  }

  [[nodiscard]] static auto maxPatches() {
    return std::numeric_limits<typename GV::IndexSet::IndexType>::max();
  }

  [[nodiscard]] static auto maxEntities() {
    return std::numeric_limits<typename GV::IndexSet::IndexType>::max();
  }

protected:

  //! @todo replace with deduce this from C++23...
  template<Concept::EntitySetPartition Super>
  requires std::is_base_of_v<UncoloredHaloMixin<GV>, Super>
  void updatePartitionHalo(const Super& super, std::size_t halo) {
    static_assert(std::same_as<GV, typename Super::EntitySet>);
    TRACE_EVENT("dune", "EntitySet::updatePartitionHalo");
    static_assert(Super::Entity::dimension == Super::EntitySet::dimension,
      "Not implemented: Partition overlap information is only available on codim == 0 entities");

    using Entity = typename Super::Entity;
    const std::size_t no_id = std::numeric_limits<std::size_t>::max();
    using mo = std::memory_order;

    // if halo is too big, mark all entities to belong to the halo
    bool max_overlap = (halo == std::numeric_limits<std::size_t>::max());
    _entity_halo->assign(_mapper->size(), max_overlap);
    if (max_overlap) return;

    std::vector<std::size_t> entity_owner(_mapper->size(), no_id);
    std::vector<char> entity_in_overlap(_mapper->size(), 0);

    auto mark_owner = [&](const Entity& entity, std::size_t id, auto) {
      Hybrid::forEach(std::make_index_sequence<Entity::dimension+1>{}, [&](auto codim){
        for (const auto& sub_entity : subEntities(entity, Dune::Codim<codim>{})) {
          auto entity_index = _mapper->index(sub_entity);
          std::atomic_ref(entity_owner[entity_index]).store(id, mo::relaxed);
        }
      });
    };

    auto mark_entity_overlap = [&](const Entity& entity, std::size_t id) {
      Hybrid::forEach(std::make_index_sequence<Entity::dimension+1>{}, [&](auto codim){
        for (const auto& sub_entity : subEntities(entity, Dune::Codim<codim>{})) {
          auto entity_index = _mapper->index(sub_entity);
          auto owner = entity_owner[entity_index];
          if (owner != id)
            std::atomic_ref(entity_in_overlap[entity_index]).store(true, mo::relaxed);
        }
      });
    };

    std::function<void(const Entity&,std::size_t,std::size_t)> mark_overlap;
    mark_overlap = [&](const Entity& entity, std::size_t id, std::size_t current_halo) {
      mark_entity_overlap(entity, id);
      if ((current_halo--) != 0)
        for (const auto& intersection : intersections(super.entitySet(), entity))
          if (intersection.neighbor())
            mark_overlap(intersection.outside(), id, current_halo);
    };

    // mark all sub entities belonging to this element as shared
    std::mutex marker_mutex;
    auto propagate_overlap = [&](const Entity& entity, auto, auto) {
      auto entity_index = _mapper->index(entity);
      Hybrid::forEach(std::make_index_sequence<Entity::dimension+1>{}, [&](auto codim){
        for (const auto& sub_entity : subEntities(entity, Dune::Codim<codim>{})) {
          auto sub_entity_index = _mapper->index(sub_entity);
          bool in_overlap = entity_in_overlap[sub_entity_index];
          if (in_overlap) {
            std::unique_lock lg{marker_mutex};
            // TODO: propagate value to other codimensions
            (*_entity_halo)[entity_index] = true;
            (*_entity_halo)[sub_entity_index] = true;
          }
        }
      });
    };

    auto for_each_element = [halo](const auto& entity_sets, auto apply){
      auto patches = std::distance(entity_sets.begin(), entity_sets.end());
      auto partitions = std::ranges::iota_view<std::size_t, std::size_t>(0, patches);
      std::for_each(std::execution::par, partitions.begin(), partitions.end(), [&](auto patch_id){
        for (const auto& entity : entity_sets[patch_id])
          apply(entity, patch_id, halo);
      });
    };

    for (const auto& concurrent_entity_sets : super.range()) {
      // phase 1: assign an owner to each used sub-entity
      for_each_element(concurrent_entity_sets, mark_owner);

      // phase 2: mark the overlap based in incompatible owners
      for_each_element(concurrent_entity_sets, mark_overlap);

      // phase 3: entities with an overlap are marked as shared region
      for_each_element(concurrent_entity_sets, propagate_overlap);
    }
  }

private:
  std::shared_ptr<Dune::MultipleCodimMultipleGeomTypeMapper<GV>> _mapper;
  std::shared_ptr<std::vector<bool>> _entity_halo;
};

} // namespace Dune::Assembler::EntitySetPartitioner::Impl

#endif // DUNE_PDELAB_COMMON_PARTITION_HALO_HH
