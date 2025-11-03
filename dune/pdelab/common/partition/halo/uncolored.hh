#ifndef DUNE_PDELAB_COMMON_PARTITION_HALO_UNCOLORED_HH
#define DUNE_PDELAB_COMMON_PARTITION_HALO_UNCOLORED_HH

#include <dune/pdelab/common/partition/halo/region.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/concepts/gridview.hh>

#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>

#include <vector>
#include <memory>
#include <atomic>
#include <version>

#if __cpp_lib_execution >= 201603L
#include <execution>
#endif

namespace Dune::PDELab::EntitySetPartition::Impl {

/**
 * @brief Uncolored halo mixin
 * @details Calculates and stores the overlap of the patches wrt halo
 * @see EntitySetPartition
 */
template<Dune::Concept::GridView GV>
class UncoloredHaloMixin {
public:

  explicit UncoloredHaloMixin(std::size_t halo_distance)
    : _halo_distance{halo_distance}
  {}

  /**
   * @brief Whether an entity is in the interior or overlap halo region of the partition
   * @param entity   The entity to test
   */
  template<Dune::Concept::Entity Entity>
  requires (Entity::codimension == 0)
  [[nodiscard]] HaloRegion haloRegion(const Entity& entity) const {
    if (_halo_distance == all_overlap_halo_region)
      return EntitySetPartition::HaloRegion::Overlap;
    else if (_halo_distance == all_interior_halo_region)
      return EntitySetPartition::HaloRegion::Interior;
    else
      return (*_entity_in_halo)[_element_mapper->index(entity)]
                ? EntitySetPartition::HaloRegion::Overlap
                : EntitySetPartition::HaloRegion::Interior;
  }

protected:

  //! @todo replace with deduce this from C++23...
  template<class Super>
  requires std::is_base_of_v<UncoloredHaloMixin<GV>, Super>
  void updateHalo(const Super& super) {
    static_assert(std::same_as<GV, typename Super::EntitySet>);
    static_assert(Super::Element::dimension == Super::EntitySet::dimension,
      "Not implemented: Partition overlap information is only available on codim == 0 entities");

    using Element = typename Super::Element;

    _element_mapper = std::make_shared<Dune::MultipleCodimMultipleGeomTypeMapper<GV>>(super.entitySet(), mcmgElementLayout());

    // if halo is too big or too small, no need to store region
    if (_halo_distance == all_overlap_halo_region or _halo_distance == all_interior_halo_region) {
      _entity_in_halo = nullptr;
      _element_mapper = nullptr;
      return;
    }

    auto all_entities_layout = [](GeometryType gt, int dimgrid) { return true; };
    Dune::MultipleCodimMultipleGeomTypeMapper<GV> all_mapper(super.entitySet(), all_entities_layout);

    // temporary storage where every entry can be updated concurrently
    // note that proxy objects in std::vector<bool> are not thread-safe
    std::vector<char> entity_in_halo(_element_mapper->size(), false);
    std::vector<std::size_t> entity_owner;
    std::vector<char> entity_in_overlap;

    // clang (<19) still does not support std::atomic_ref, so use __atomic builtins if not available
#if __cpp_lib_atomic_ref >= 201806L
    static_assert(std::atomic_ref<char>::is_always_lock_free);
    static_assert(std::atomic_ref<std::size_t>::is_always_lock_free);
    constexpr auto atomic_load_relaxed = []<class T>(T& obj) {
      return std::atomic_ref<T>{obj}.load(std::memory_order::relaxed);
    };
    constexpr auto atomic_store_relaxed = []<class T>(T& obj, T val) {
      std::atomic_ref{obj}.store(val, std::memory_order::relaxed);
    };
#else
    static_assert(__atomic_always_lock_free(sizeof(char), 0));
    static_assert(__atomic_always_lock_free(sizeof(std::size_t), 0));
    constexpr auto atomic_load_relaxed = []<class T>(T& obj) {
      T tmp{};
      __atomic_load(&obj, &tmp, __ATOMIC_RELAXED);
      return tmp;
    };
    constexpr auto atomic_store_relaxed = []<class T>(T& obj, T val) {
      __atomic_store(&obj, &val, __ATOMIC_RELAXED);
    };
#endif

    // mark every entity with a unique owner
    auto mark_owner = [&](const Element& entity, std::size_t id, auto) {
      Hybrid::forEach(std::make_index_sequence<Element::dimension+1>{}, [&](auto codim){
        for (const auto& sub_entity : subEntities(entity, Dune::Codim<codim>{})) {
          auto entity_index = all_mapper.index(sub_entity);
          atomic_store_relaxed(entity_owner[entity_index], id);
        }
      });
    };

    // if an entity is owned by different patches, it is in the overlap
    auto mark_entity_overlap = [&](const Element& entity, std::size_t id) {
      Hybrid::forEach(std::make_index_sequence<Element::dimension+1>{}, [&](auto codim){
        for (const auto& sub_entity : subEntities(entity, Dune::Codim<codim>{})) {
          auto entity_index = all_mapper.index(sub_entity);
          auto owner = atomic_load_relaxed(entity_owner[entity_index]);
          if (owner != id)
            atomic_store_relaxed(entity_in_overlap[entity_index], char{true});
        }
      });
    };

    // propagate the overlap on neighboring entities recursively until the halo distance is reached
    std::function<void(const Element&,std::size_t,std::size_t)> mark_overlap;
    mark_overlap = [&](const Element& entity, std::size_t id, std::size_t halo_distance) {
      mark_entity_overlap(entity, id);
      if (halo_distance != 0) {
        for (const auto& intersection : intersections(super.entitySet(), entity))
          if (intersection.neighbor())
            mark_overlap(intersection.outside(), id, halo_distance-1);
      }
    };

    // if any of the sub-entities is in the overlap, the entity is in the halo
    auto propagate_overlap = [&](const Element& entity, auto, auto) {
      auto entity_index = _element_mapper->index(entity);
      Hybrid::forEach(std::make_index_sequence<Element::dimension+1>{}, [&](auto codim){
        for (const auto& sub_entity : subEntities(entity, Dune::Codim<codim>{})) {
          auto sub_entity_index = all_mapper.index(sub_entity);
          bool in_overlap = atomic_load_relaxed(entity_in_overlap[sub_entity_index]);
          if (in_overlap) {
            atomic_store_relaxed(entity_in_halo[entity_index], char{true});
          }
        }
      });
    };

    // helper to run a function on all entities of all patches in parallel
    auto for_each_element = [halo_distance = _halo_distance](const auto& entity_sets, auto apply){
      auto patches = std::distance(entity_sets.begin(), entity_sets.end());
      auto partitions = Dune::range(std::size_t{0}, static_cast<std::size_t>(patches));
      std::for_each(
#if __cpp_lib_execution >= 201603L
        std::execution::par,
#endif
        partitions.begin(), partitions.end(), [&](auto patch_id){
        for (const auto& entity : entity_sets[patch_id])
          apply(entity, patch_id, halo_distance);
      });
    };

    for (const auto& concurrent_entity_sets : super) {
      // clean up overlap info
      entity_owner.assign(all_mapper.size(), std::numeric_limits<std::size_t>::max());
      entity_in_overlap.assign(all_mapper.size(), false);

      // phase 1: assign an owner to each used sub-entity
      for_each_element(concurrent_entity_sets, mark_owner);

      // phase 2: mark the overlap based in incompatible owners
      for_each_element(concurrent_entity_sets, mark_overlap);

      // phase 3: entities with an overlap are marked as shared region
      for_each_element(concurrent_entity_sets, propagate_overlap);
    }

    _entity_in_halo = std::make_shared<std::vector<bool>>(entity_in_halo.size());
    std::copy(entity_in_halo.begin(), entity_in_halo.end(), _entity_in_halo->begin());
  }

private:
  std::shared_ptr<Dune::MultipleCodimMultipleGeomTypeMapper<GV>> _element_mapper;
  std::shared_ptr<std::vector<bool>> _entity_in_halo;
  std::size_t _halo_distance;
};

} // namespace Dune::PDELab::EntitySetPartition::Impl

#endif // DUNE_PDELAB_COMMON_PARTITION_HALO_UNCOLORED_HH
