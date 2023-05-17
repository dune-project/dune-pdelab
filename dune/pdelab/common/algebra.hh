#ifndef DUNE_PDELAB_COMMON_ALGEBRA_HH
#define DUNE_PDELAB_COMMON_ALGEBRA_HH

#include <dune/pdelab/common/container_entry.hh>
#include <dune/pdelab/common/container_traversal.hh>

namespace Dune::PDELab::inline Experimental {

template<class ExecutionPolicy, class Y, class Field, class X>
requires std::is_execution_policy_v<std::decay_t<ExecutionPolicy>>
void axpy(ExecutionPolicy&& policy, Y& y, Field alpha, const X& x){
  forEachContainerEntry(policy, x, [&](auto x_i, auto i){
    containerEntry(y, i) += alpha*x_i;
  });
}


template<class Y, class Field, class X>
void axpy(Y& y, Field alpha, const X& x){
  forEachContainerEntry(std::execution::seq, x, [&](auto x_i, auto i){
    containerEntry(y, i) += alpha*x_i;
  });
}


namespace Impl {

template<class ExecutionPolicy, class Domain, class Map, class Range, class Callable>
void linearTransformation(
  ExecutionPolicy,
  Domain&&    domain,
  Map&&       map,
  Range&&     range,
  Callable&&  callable,
  auto map_mi, auto domain_mi, auto range_mi)
{
  // This function takes a map (i.e. tensor product between domain and range, think
  // about it as a tensor or matrix) and invokes the callable on each leaf entry
  // of the map/domain/range (e.g. matrix values). A recursive loop is needed on all
  // nested loops within them map.
  // Each map loop should match with either a loop in domain or range.
  // Depending on the constness of domain/range, we forward the execution policy or not.
  // If one of range/domain is const, loops to these ranges need be serialized.
  // (imagine a matrix-vector multiplication range[i] += map[i][j]*domain[j]).
  // If none of both domain/range are const, all loops need to be serialized.

  // TODO check that the following assertion holds (without calling split)
  // auto [old_domain_mi, old_range_mi] = containerIndexSplit(map, map_mi);
  // static_assert(std::same_as<decltype(old_domain_mi), decltype(domain_mi)>);
  // static_assert(std::same_as<decltype(old_range_mi),  decltype(range_mi)>);

  if constexpr (Concept::Range<Map>) {
    using namespace Dune::Indices;
    using std::as_const;
    auto [t_domain_mi, t_range_mi] = containerIndexSplit(as_const(map), push_back(map_mi, _0));
    const auto domain_loop = std::same_as<decltype(t_range_mi),  decltype(range_mi)>;
    const auto range_loop  = std::same_as<decltype(t_domain_mi), decltype(domain_mi)>;
    static_assert(domain_loop xor range_loop);

    auto loop_policy = []{
      constexpr auto const_map    = std::is_const_v<std::remove_reference_t<Map>>;
      constexpr auto const_domain = std::is_const_v<std::remove_reference_t<Domain>>;
      constexpr auto const_range  = std::is_const_v<std::remove_reference_t<Range>>;
      if constexpr (!const_map)
        return std::execution::seq;
      if constexpr (const_domain and const_range)
        return ExecutionPolicy{};
      if constexpr (!const_domain and !const_range)
        return std::execution::seq;
      else if constexpr ((range_loop and !const_range) or (domain_loop and !const_domain))
        return ExecutionPolicy{};
      else
        return std::execution::seq;
    }();

    PDELab::forEach(loop_policy, std::forward<Map>(map), [&]<class MapEntry>(MapEntry&& map_entry, auto i){
      // containerEntry only works if index split does not reorder the indices
      auto new_map_mi = push_back(map_mi, i);
      auto [new_domain_mi, new_range_mi] = containerIndexSplit(as_const(map), new_map_mi);
      if constexpr (domain_loop) {
        // assert(new_domain_mi == push_back(domain_mi, i));
        Impl::linearTransformation(
          ExecutionPolicy{},
          PDELab::containerEntry(std::forward<Domain>(domain), TypeTree::treePath(i)),
          std::forward<MapEntry>(map_entry),
          std::forward<Range>(range),
          std::forward<Callable>(callable),
          new_map_mi, new_domain_mi, new_range_mi);
      } else {
        // assert(new_range_mi == push_back(range_mi, i));
        Impl::linearTransformation(
          ExecutionPolicy{},
          std::forward<Domain>(domain),
          std::forward<MapEntry>(map_entry),
          PDELab::containerEntry(std::forward<Range>(range), TypeTree::treePath(i)),
          std::forward<Callable>(callable),
          new_map_mi, new_domain_mi, new_range_mi);
      }
    });
  } else {
    std::invoke(
      std::forward<Callable>(callable),
      std::forward<Domain>(domain),
      std::forward<Map>(map),
      std::forward<Range>(range)
    );
  }
}

} // namespace Impl

// TODO document, category value is very important here!!
template<class ExecutionPolicy, class Domain, class Map, class Range, class Callable>
requires std::is_execution_policy_v<std::decay_t<ExecutionPolicy>>
void linearTransformation(ExecutionPolicy&& policy, Domain&& domain, Map&& map, Range&& range, Callable&& callable){
  Impl::linearTransformation(
    policy,
    std::forward<Domain>(domain),
    std::forward<Map>(map),
    std::forward<Range>(range),
    std::forward<Callable>(callable),
    TypeTree::treePath(), TypeTree::treePath(), TypeTree::treePath()
  );
}

template<class Domain, class Map, class Range, class Callable>
requires (not std::is_execution_policy_v<std::decay_t<Domain>>)
void linearTransformation(Domain&& domain, Map&& map, Range&& range, Callable&& callable){
  linearTransformation(
    std::execution::seq,
    std::forward<Domain>(domain),
    std::forward<Map>(map),
    std::forward<Range>(range),
    std::forward<Callable>(callable)
  );
}

template<class ExecutionPolicy, class Domain, class Map, class Range>
requires std::is_execution_policy_v<std::decay_t<ExecutionPolicy>>
void linearTransformation(ExecutionPolicy&& policy, const Domain& domain, const Map& map, Range& range){
  linearTransformation(policy, domain, map, range, [](const auto& x, const auto& a, auto& y){y += a*x;});
}

template<class Domain, class Map, class Range>
void linearTransformation(const Domain& domain, const Map& map, Range& range){
  linearTransformation(std::execution::seq, domain, map, range, [](const auto& x, const auto& a, auto& y){y += a*x;});
}

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_ALGEBRA_HH
