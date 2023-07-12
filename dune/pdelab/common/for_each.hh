#ifndef DUNE_PDELAB_COMMON_FOR_EACH_HH
#define DUNE_PDELAB_COMMON_FOR_EACH_HH

#include <dune/pdelab/concepts/range.hh>
#include <dune/pdelab/concepts/treenode.hh>

#include <dune/pdelab/common/multiindex.hh>

#include <dune/common/indices.hh>

#if HAVE_TBB
#include <tbb/parallel_for.h>
#include <tbb/task_group.h>
#endif

#include <utility>
#include <execution>

namespace Dune {

// forward declaration
template<class, int> class FieldVector;

// forward declaration
template<class, int, int> class FieldMatrix;

} // namespace Dune

namespace Dune::PDELab::inline Experimental {

// helper class to exclude types from using a parallel forEach
template<class> struct ExcludeParallelForEach : std::false_type {};
template<class T, std::size_t n> struct ExcludeParallelForEach<Dune::FieldVector<T,n>> : std::true_type {};
template<class T, std::size_t rows, std::size_t cols> struct ExcludeParallelForEach<Dune::FieldMatrix<T,rows, cols>> : std::true_type {};

namespace Impl {

inline namespace Default {

//! Default for each implementation
template<class Container, class Callable>
requires Concept::Range<std::remove_cvref_t<Container>>
constexpr void forEach(std::execution::sequenced_policy, Container&& container, Callable&& at_value)
{
  using UContainer = std::remove_cvref_t<Container>;

  auto invoke = [&at_value]<class Value, class Index>(Value&& value, Index index){
    static_assert(std::invocable<Callable&&, Value&&, Index> || std::invocable<Callable&&, Value&&>);
    if constexpr (std::invocable<Callable&&, Value&&, Index>)
      std::invoke(std::forward<Callable>(at_value), std::forward<Value>(value), index);
    else
      std::invoke(std::forward<Callable>(at_value), std::forward<Value>(value));
  };

  if constexpr (Concept::SparseDynamicRange<UContainer> ) {
    using std::end;
    auto end_it = end(container);
    using std::begin;
    for (auto it = begin(container); it != end_it; ++it)
      invoke(*it, it.index());
  } else if constexpr (Concept::DenseDynamicRange<UContainer>) {
    std::size_t i = 0;
    for (auto&& entry : container)
      invoke(std::forward<decltype(entry)>(entry), i++);
  } else if constexpr (Concept::DenseStaticRange<UContainer>) {
    [&]<class I, std::size_t... i>(std::integer_sequence<I,i...>) {
      (invoke(container[std::integral_constant<I,i>{}], std::integral_constant<I,i>{}), ...);
    }(std::make_index_sequence<UContainer::size()>{});
  } else if constexpr (Concept::SparseStaticRange<UContainer>) {
    [&]<class I, std::size_t... i>(std::integer_sequence<I,i...>) {
      (invoke(container[std::integral_constant<I,i>{}], std::integral_constant<I,i>{}), ...);
    }(UContainer::enumerate());
  } else {
    static_assert(AlwaysFalse<Container>{});
  }
}


#if HAVE_TBB

//! Default for each implementation of parallel policy
template<class Container, class Callable>
requires Concept::Range<std::remove_cvref_t<Container>>
constexpr void forEach(std::execution::parallel_policy, Container&& container, Callable&& at_value)
{
  using UContainer = std::remove_cvref_t<Container>;

  auto invoke = [&at_value]<class Value, class Index>(Value&& value, Index index){
    static_assert(std::invocable<Callable&&, Value&&, Index> || std::invocable<Callable&&, Value&&>);
    if constexpr (std::invocable<Callable&&, Value&&, Index>)
      std::invoke(std::forward<Callable>(at_value), std::forward<Value>(value), index);
    else
      std::invoke(std::forward<Callable>(at_value), std::forward<Value>(value));
  };

  // forward field vector/matrix to sequential loop
  if constexpr (ExcludeParallelForEach<UContainer>::value) {
    forEach(std::execution::seq, std::forward<Container>(container), std::forward<Callable>(at_value));
  } else {
    tbb::task_group task_group;

    if constexpr (Concept::SparseDynamicRange<UContainer> or Concept::DenseDynamicRange<UContainer>) {

      // simple forward range that, on split, forward advances until middle is found
      using T = decltype(std::begin(container));
      class forward_range {
        T _begin;
        T _end;
        std::size_t _size;
      public:
        T begin() const {return _begin;}
        T end() const {return _end;}
        bool empty() const {return _begin==_end;}
        bool is_divisible() const {return _size>1;}
        forward_range( T first, T last ) : _begin(first), _end(last), _size(std::distance(first, last)) {}
        forward_range(forward_range& r, tbb::split) {
          size_t h = r._size/2;
          _end = r._end;
          _begin = r._begin;
          using std::advance;
          advance(_begin, h);
          _size = r._size-h;
          r._end = _begin;
          r._size = h;
        }
      };

      auto tbb_range = [&container]{
        // if possible use native TBB blocked range otherwise use custom forward range
        if constexpr (requires { tbb::blocked_range{std::begin(container), std::end(container), 1000}; })
          return tbb::blocked_range{std::begin(container), std::end(container), 1000};
        else
          return forward_range{std::begin(container), std::end(container)};
      }();
      if constexpr (Concept::SparseDynamicRange<UContainer>) {
        tbb::parallel_for(tbb_range ,[&](auto range){
          for (auto it = std::begin(range); it != std::end(range); ++it)
            invoke(*it, it.index());
        });
      } else if constexpr (Concept::DenseDynamicRange<UContainer>) {
        tbb::parallel_for(tbb_range, [&](auto range){
          std::size_t i = std::distance(std::begin(container), range.begin());
          for (auto&& entry : range)
            invoke(std::forward<decltype(entry)>(entry), i++);
        });
      }
    } else if constexpr (Concept::DenseStaticRange<UContainer>) {
      [&]<class I, std::size_t... i>(std::integer_sequence<I,i...>) {
        (task_group.run([&]{invoke(container[std::integral_constant<I,i>{}], std::integral_constant<I,i>{});}), ...);
      }(std::make_index_sequence<UContainer::size()>{});
    } else if constexpr (Concept::SparseStaticRange<UContainer>) {
      [&]<class I, std::size_t... i>(std::integer_sequence<I,i...>) {
        (task_group.run([&]{invoke(container[std::integral_constant<I,i>{}], std::integral_constant<I,i>{});}), ...);
      }(UContainer::enumerate());
    } else {
      static_assert(AlwaysFalse<Container>{});
    }
  }
}

//! Default for each implementation of parallel unsequenced policy
template<class Container, class Callable>
requires Concept::Range<std::remove_cvref_t<Container>>
constexpr void forEach(std::execution::parallel_unsequenced_policy, Container&& container, Callable& at_value)
{
  forEach(std::execution::par, std::forward<Container>(container), std::forward<Callable>(at_value));
}

#endif // HAVE_TBB

/**
 * @brief Traverse each entry of a tree container
 * @details Functors may accept one or two values. The first one being the child
 * node to be evaluated, and the second one the index of such entry
 *
 * @tparam Container          Tree container to be traversed
 * @tparam Callable           Functor to be applied at each child
 */
template<class Container, class Callable>
requires Concept::ParentTreeNode<std::remove_cvref_t<Container>>
constexpr void forEach(std::execution::sequenced_policy, Container&& container, Callable&& at_value)
{
  auto invoke = [&at_value]<class Value, class Index>(Value&& value, Index index){
    static_assert(std::invocable<Callable&&, Value&&, Index> || std::invocable<Callable&&, Value&&>);
    if constexpr (std::invocable<Callable&&, Value&&, Index>)
      std::invoke(std::forward<Callable>(at_value), std::forward<Value>(value), index);
    else
      std::invoke(std::forward<Callable>(at_value), std::forward<Value>(value));
  };

  if constexpr (Concept::TupleTreeNode<Container>)
    Dune::unpackIntegerSequence(
      [&](auto... i) { (invoke(container.child(i), i), ...); },
      std::make_index_sequence<std::remove_cvref_t<Container>::degree()>{});
  else
    for (std::size_t i = 0; i != container.degree(); ++i)
      invoke(container.child(i), i);
}

//! Default for each implementation
template<class ExecutionPolicy, class Container, class Callable>
requires (std::is_execution_policy_v<std::decay_t<ExecutionPolicy>> && !std::same_as<std::decay_t<ExecutionPolicy>, std::execution::sequenced_policy>)
constexpr void forEach(ExecutionPolicy&& policy, Container&& container, Callable&& at_value)
{
  // if policy is not tag dispatched, use sequential policy
  forEach(std::execution::seq, std::forward<Container>(container), std::forward<Callable>(at_value));
}

//! Default for each implementation
template<class Container, class Callable>
constexpr void forEach(Container&& container, Callable&& at_value)
{
  // use sequential policy if none is asked
  forEach(std::execution::seq, std::forward<Container>(container), std::forward<Callable>(at_value));
}

} // namespace Default

//! Customization point object on for each traversal
struct ForEach {

  /**
   * @brief Applies the given function object at_value to every entry in the container.
   *
   * @param container    Container to iterate
   * @param at_value     Function object that receives every entity, and optionally an iteration index
   */
  template<class Container, class Callable>
  constexpr void operator()(Container&& container, Callable&& at_value) const {
    forEach(std::forward<Container>(container), std::forward<Callable>(at_value));
  }


  /**
   * @brief Applies the given function object at_value to every entry in the container.
   *
   * @param policy       Execution policy to use. @see std::execition.
   * @param container    Container to iterate
   * @param at_value     Function object that receives every entity, and optionally an iteration index
   */
  template<class ExecutionPolicy, class Container, class Callable>
  requires std::is_execution_policy_v<std::decay_t<ExecutionPolicy>>
  constexpr void operator()(ExecutionPolicy&& policy, Container&& container, Callable&& at_value) const {
    forEach(policy, std::forward<Container>(container), std::forward<Callable>(at_value));
  }
};

} // namespace Impl


namespace Default = Impl::Default;

//! Customization point object on for each traversal
inline const Impl::ForEach         forEach {};


} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_FOR_EACH_HH
