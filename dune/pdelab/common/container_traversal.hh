#ifndef DUNE_PDELAB_COMMON_CONTAINER_TRAVERSAL_HH
#define DUNE_PDELAB_COMMON_CONTAINER_TRAVERSAL_HH

#include <dune/pdelab/concepts/multiindex.hh>

#include <dune/pdelab/common/for_each.hh>

#include <dune/typetree/treepath.hh>
#include <dune/typetree/traversal.hh>

#include <utility>
#include <execution>

namespace Dune::PDELab::inline Experimental {

/**
 * @brief Traverse each entry of a nested container
 * @details Functors may accept one or two values. The first one being the entry
 * to be evaluated, and the second one the multi-index of such entry
 *
 * @tparam ExcecutionPolicy   The execution policy. @see std::execution
 * @tparam Container          Recursive container to be traversed
 * @tparam AtValue            Functor to be applied at each value
 * @tparam PreValue           Functor to be applied before each value
 * @tparam PostValue          Functor to be applied after each value
 * @tparam MultiIndex         Multi-index representing the current position of the container
 */
template<
  class ExcecutionPolicy,
  class Container,
  class AtValue,
  class PreValue,
  class PostValue,
  Concept::MultiIndex MultiIndex>
requires (std::is_execution_policy_v<std::decay_t<ExcecutionPolicy>>)
constexpr void
forEachContainerEntry(ExcecutionPolicy&& policy,
            Container&& container,
            AtValue&& at_value,
            PreValue&& pre_value,
            PostValue&& post_value,
            MultiIndex multiindex)
{
  auto invoke = [&multiindex]<class Callable, class Entry>(Callable&& callable, Entry&& entry){
    static_assert(std::invocable<Callable&&, Entry&&, const MultiIndex&> or std::invocable<Callable&&, Entry&&>);
    if constexpr (std::invocable<Callable&&, Entry&&, const MultiIndex&>)
      callable(std::forward<Entry>(entry), std::as_const(multiindex));
    else
      callable(std::forward<Entry>(entry));
  };

  invoke(std::forward<PreValue>(pre_value), std::forward<Container>(container));

  if constexpr (Concept::Range<std::remove_cvref_t<Container>> )
    Dune::PDELab::forEach(std::forward<ExcecutionPolicy>(policy), std::forward<Container>(container), [&]<class Entry>(Entry&& entry, auto i){
      forEachContainerEntry(
        std::forward<ExcecutionPolicy>(policy),
        std::forward<Entry>(entry),
        std::forward<AtValue>(at_value),
        std::forward<PreValue>(pre_value),
        std::forward<PostValue>(post_value),
        push_back(multiindex, i)
      );
    });
  else
    invoke(std::forward<AtValue>(at_value), std::forward<Container>(container));

  invoke(std::forward<PostValue>(post_value), std::forward<Container>(container));
}


/**
 * @brief Traverse each entry of a nested container
 * @details Functors may accept one or two values. The first one being the entry
 * to be evaluated, and the second one the multi-index of such entry
 *
 * @tparam Container          Recursive container to be traversed
 * @tparam AtValue            Functor to be applied at each value
 * @tparam PreValue           Functor to be applied before each value
 * @tparam PostValue          Functor to be applied after each value
 * @tparam MultiIndex         Multi-index representing the current position of the container
 */
template<class Container, class AtValue, class PreValue, class PostValue>
constexpr void
forEachContainerEntry(
            Container&& container,
            AtValue&& at_value,
            PreValue&& pre_value,
            PostValue&& post_value,
            Concept::MultiIndex auto multiindex)
{
  forEachContainerEntry(
    std::execution::seq,
    std::forward<Container>(container),
    std::forward<AtValue>(at_value),
    std::forward<PreValue>(pre_value),
    std::forward<PostValue>(post_value),
    multiindex
  );
}

/**
 * @brief Traverse each entry of a nested container
 * @details Functors may accept one or two values. The first one being the entry
 * to be evaluated, and the second one the multi-index of such entry
 *
 * @tparam ExcecutionPolicy   The execution policy. @see std::execution
 * @tparam Container          Recursive container to be traversed
 * @tparam AtValue            Functor to be applied at each value
 */
template<class ExcecutionPolicy, class Container, class AtValue>
requires (std::is_execution_policy_v<std::decay_t<ExcecutionPolicy>>)
void forEachContainerEntry(ExcecutionPolicy&& policy, Container&& container, AtValue&& at_value) {
  forEachContainerEntry(
    std::forward<ExcecutionPolicy>(policy),
    std::forward<Container>(container),
    std::forward<AtValue>(at_value),
    TypeTree::NoOp{},
    TypeTree::NoOp{},
    TypeTree::treePath()
  );
}


/**
 * @brief Traverse each entry of a nested container
 * @details Functors may accept one or two values. The first one being the entry
 * to be evaluated, and the second one the multi-index of such entry
 *
 * @tparam Container          Recursive container to be traversed
 * @tparam AtValue            Functor to be applied at each value
 */
template<class Container, class AtValue>
void forEachContainerEntry(Container&& container, AtValue&& at_value) {
  forEachContainerEntry(
    std::execution::seq,
    std::forward<Container>(container),
    std::forward<AtValue>(at_value),
    TypeTree::NoOp{},
    TypeTree::NoOp{},
    TypeTree::treePath()
  );
}

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_CONTAINER_TRAVERSAL_HH
