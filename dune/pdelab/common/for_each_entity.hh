#ifndef DUNE_PDELAB_COMMON_FOR_EACH_ENTITY_HH
#define DUNE_PDELAB_COMMON_FOR_EACH_ENTITY_HH

#include <dune/pdelab/common/partition/identity.hh>
#include <dune/pdelab/common/for_each.hh>

#include <dune/pdelab/common/trace.hh>

#include <dune/pdelab/concepts/entityset_partition.hh>

#include <execution>
#include <concepts>

namespace Dune::PDELab::inline Experimental {

template<typename ExecutionPolicy, Concept::EntitySetPartition EntitySetPartition>
requires std::is_execution_policy_v<ExecutionPolicy>
void forEachEntity(const ExecutionPolicy& policy, const EntitySetPartition& entity_set_partition, std::invocable<const typename EntitySetPartition::Entity&> auto fapply) {
  for (const auto& label : entity_set_partition.range()) {
    TRACE_EVENT("dune", "forEachEntity::Label");
    std::for_each(policy, label.begin(), label.end(), [&](const auto& patch) mutable {
      auto _fapply = fapply;
      // TRACE_EVENT("dune", "forEachEntity::Patch");
      for (auto&& entity : patch)
        _fapply(entity);
    });
  }
}

template<Concept::EntitySetPartition EntitySetPartition>
void forEachEntity(const EntitySetPartition& entity_set_partition, std::invocable<const typename EntitySetPartition::Entity&> auto fapply) {
  forEachEntity(std::execution::seq, entity_set_partition, fapply);
}

template<Dune::Concept::GridView GridView>
void forEachEntity(const GridView& entity_set, std::invocable<const typename GridView::template Codim<0>::Entity&> auto fapply) {
  forEachEntity(EntitySetPartitioner::Identity{entity_set}, fapply);
}

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_COMMON_FOR_EACH_ENTITY_HH
