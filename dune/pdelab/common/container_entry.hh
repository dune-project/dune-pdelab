#ifndef DUNE_PDELAB_COMMON_CONTAINER_ENTRY_HH
#define DUNE_PDELAB_COMMON_CONTAINER_ENTRY_HH

#include <dune/pdelab/concepts/indexable.hh>
#include <dune/pdelab/concepts/multiindex.hh>

#include <dune/pdelab/common/multiindex.hh>

#include <dune/typetree/treepath.hh>

#include <utility>

namespace Dune::PDELab::inline Experimental {

namespace Impl {

inline namespace Default {

/**
 * @brief Access an entry in a container
 * @details When multi-index is empty returns container
 *
 * @tparam Container        Container type
 * @tparam MultiIndex       Multi-index type
 * @param container         Container with the desired entry
 * @param multiindex        Multi-index to reach the target entry
 * @return decltype(auto)   Entry at multi-index
 */
template<class Container, Concept::FixedSizeMultiIndex MultiIndex>
requires (Concept::FixedSizeMultiIndex<MultiIndex> && MultiIndex::size() == 0)
constexpr decltype(auto) containerEntry(Container&& container, MultiIndex) {
  return std::forward<Container>(container);
}

/**
 * @brief Access an entry in a container
 * @details Use front-most index to return a bracket-indexed entry from the container
 *
 * @tparam Container        Container type
 * @tparam MultiIndex       Multi-index type
 * @param container         Container with the desired entry
 * @param multiindex        Multi-index to reach the target entry
 * @return decltype(auto)   Entry at multi-index
 */
template<class Container, Concept::FixedSizeMultiIndex MultiIndex>
constexpr decltype(auto) containerEntry(Container&& container, MultiIndex multiindex)
requires (Concept::FixedSizeMultiIndex<MultiIndex> && MultiIndex::size() != 0 && Concept::BracketIndexable<decltype(container),decltype(front(multiindex))>)
{
  auto index = front(multiindex);
  return containerEntry(std::forward<Container>(container)[index], pop_front(multiindex));
}


template<class Container, Concept::MultiIndex MultiIndex>
constexpr decltype(auto) containerEntry(Container&& container, MultiIndex multiindex)
requires (not Concept::FixedSizeMultiIndex<MultiIndex>)
{
  if constexpr (Concept::BracketIndexable<decltype(container),decltype(front(multiindex))>) {
    auto index = front(multiindex);
    return containerEntry(std::forward<Container>(container)[index], pop_front(multiindex));
  } else {
    return std::forward<Container>(container);
  }
}

} // namespace Default

//! Customization point object on container entry
struct ContainerEntry {

  /**
   * @brief Access an entry in a container
   *
   * @tparam Container        Container type
   * @tparam MultiIndex       Multi-index type
   * @param container         Container with the desired entry
   * @param multiindex        Multi-index to reach the target entry
   * @return decltype(auto)   Entry at multi-index
   */
  template<class Container>
  constexpr decltype(auto) operator()(Container&& container, Concept::MultiIndex auto multiindex) const {
    return containerEntry(std::forward<Container>(container), multiindex);
  }
};

} // namespace Impl

namespace Default = Impl::Default;

//! Customization point object on container entry
inline const Impl::ContainerEntry         containerEntry {};


namespace Impl {

// default implementation TODO CustomizationPoint
template<Concept::MultiIndex Path>
constexpr auto containerIndexSplit(const auto& map, Path path) {
  // split interleaved indices
  auto domain = unpackIntegerSequence([&](auto... i) {
    return TypeTree::treePath(path[std::integral_constant<std::size_t,i*2+1>{}]...);
  }, std::make_index_sequence<(Path::size()/2)>{});

  auto range = unpackIntegerSequence([&](auto... i) {
    return TypeTree::treePath(path[std::integral_constant<std::size_t,i*2>{}]...);
  }, std::make_index_sequence<(Path::size()/2+Path::size()%2)>{});

  return std::make_pair(domain, range);
}

// default implementation TODO CustomizationPoint
template<Concept::MultiIndex DomainPath, Concept::MultiIndex RangePath>
constexpr auto containerIndexMerge(const auto& map, DomainPath domain, RangePath range) {
  static_assert(DomainPath::size() == RangePath::size());
  if constexpr (DomainPath::size() == 0)
    return TypeTree::treePath();
  else {
    // merge indices interleaved one to the other
    return unpackIntegerSequence([&](auto... i) {
      return join(TypeTree::treePath(range[i], domain[i])...);
    }, std::make_index_sequence<RangePath::size()>{});
  }
}

} // namespace Impl

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_CONTAINER_ENTRY_HH
