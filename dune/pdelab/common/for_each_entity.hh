#ifndef DUNE_PDELAB_COMMON_FOR_EACH_ENTITY_HH
#define DUNE_PDELAB_COMMON_FOR_EACH_ENTITY_HH

#include <dune/pdelab/common/partition/identity.hh>
#include <dune/pdelab/common/for_each.hh>
#include <dune/pdelab/common/execution.hh>
#include <dune/pdelab/common/trace.hh>

#include <dune/pdelab/concepts/entityset_partition.hh>

#include <concepts>

namespace Dune::PDELab::inline Experimental {

template<typename ExecutionPolicy, Concept::EntitySetPartition EntitySetPartition>
requires Execution::is_execution_policy_v<ExecutionPolicy>
void forEachEntity(const ExecutionPolicy& policy, const EntitySetPartition& entity_set_partition, std::invocable<const typename EntitySetPartition::Entity&> auto fapply) {
  for (const auto& label : entity_set_partition.range()) {
    TRACE_EVENT("dune", "forEachEntity::Label");
    forEach(policy, label, [&](const auto& patch) mutable {
      TRACE_EVENT_BEGIN("dune", "forEachEntity::Patch");
      auto _fapply = fapply;
      for (auto&& entity : patch)
        _fapply(entity);
      TRACE_EVENT_END("dune");
    });
  }
}

template<Concept::EntitySetPartition EntitySetPartition>
void forEachEntity(const EntitySetPartition& entity_set_partition, std::invocable<const typename EntitySetPartition::Entity&> auto fapply) {
  forEachEntity(Execution::seq, entity_set_partition, fapply);
}

template<Dune::Concept::GridView GridView>
void forEachEntity(const GridView& entity_set, std::invocable<const typename GridView::template Codim<0>::Entity&> auto fapply) {
  forEachEntity(EntitySetPartitioner::Identity{entity_set}, fapply);
}

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_COMMON_FOR_EACH_ENTITY_HH
