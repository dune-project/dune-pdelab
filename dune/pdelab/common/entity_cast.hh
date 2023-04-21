#ifndef DUNE_PDELAB_COMMON_ENTITY_CAST_HH
#define DUNE_PDELAB_COMMON_ENTITY_CAST_HH

#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/gridview.hh>

namespace Dune::PDELab::inline Experimental {

namespace Impl {

inline namespace Default {

/**
 * @brief Default implementation of entity cast
 * In case of multidomain grids, this implementation automatically casts
 * multidomain entities into a subdomain entities.
 * @param target_entity_set  The entity set into which cast the entity
 * @param entity             Entity to cast
 * @return decltype(auto)    The new entity
 */
template<Dune::Concept::GridView TargetEntitySet, class Entity>
requires Dune::Concept::Entity<std::remove_cvref_t<Entity>>
decltype(auto) entityCast(const TargetEntitySet& target_entity_set, Entity&& source_entity) {
  // multidomain to subdomain conversion
  if constexpr (requires { target_entity_set.grid().subDomainEntity(std::forward<Entity>(source_entity)); } )
    return target_entity_set.grid().subDomainEntity(std::forward<Entity>(source_entity));
  else
    return source_entity;
}

} // namespace Default

//! Customization point object on entity casts
struct EntityCast {
  /**
   * @brief Invokes an entity cast
   * The source entity must be casted into an entity belonging to the
   * target entity set.
   *
   * @param target_entity_set   Target entity set to cast entity to
   * @param source_entity       Source entity
   * @return decltype(auto)     Casted entity
   */
  template<Dune::Concept::GridView TargetEntitySet, class Entity>
  requires Dune::Concept::Entity<std::remove_cvref_t<Entity>>
  decltype(auto) operator()(const TargetEntitySet& target_entity_set, Entity&& source_entity) const {
    return entityCast(target_entity_set, std::forward<Entity>(source_entity));
  }
};

} // namespace Impl

namespace Default = Impl::Default;

//! Customization point object for entity casts
inline const Impl::EntityCast             entityCast {};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_ENTITY_CAST_HH
