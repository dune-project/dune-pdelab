#ifndef DUNE_PDELAB_COMMON_LOCAL_CONTAINER_ACCESS_HH
#define DUNE_PDELAB_COMMON_LOCAL_CONTAINER_ACCESS_HH

#include <dune/pdelab/concepts/local_basis.hh>

#include <dune/pdelab/common/container_entry.hh>

#include <utility>
#include <type_traits>

namespace Dune::PDELab::inline Experimental {

namespace Impl {

//! Customization point object on local container entry
struct LocalContainerEntry {

  /**
   * @brief Access a local degree of freedom entry in a container
   * @details Equivalent to `containerEntry(container, node.index(dof))`. Prefer
   * this version as local space nodes may know faster versions of this function (e.g.
   * due different internal data representation of the multi-index).
   *
   * @param container         Container from where to get the entry
   * @param leaf_space_node   Leaf local space node with knwoledge of degree of freedom position
   * @param dof               Index of the degree of freedom to access
   * @return constexpr decltype(auto)  The container entry
   */
  template<class Container, Concept::LocalIndexSetLeaf LocalIndexSetLeaf>
  constexpr decltype(auto) operator()(Container&& container, const LocalIndexSetLeaf& leaf_space_node, typename LocalIndexSetLeaf::size_type dof) const {
    return localContainerEntry(std::forward<Container>(container), leaf_space_node, dof);
  }

  template<class Container, Concept::LocalIndexSetLeaf TestLocalIndexSetLeaf, Concept::LocalIndexSetLeaf TrailLocalIndexSetLeaf>
  constexpr decltype(auto) operator()(
    Container&& container,
    const TestLocalIndexSetLeaf& test_node,
    typename TestLocalIndexSetLeaf::size_type test_dof,
    const TrailLocalIndexSetLeaf& trial_node,
    typename TrailLocalIndexSetLeaf::size_type trial_dof) const
  {
    return localContainerEntry(std::forward<Container>(container), test_node, test_dof, trial_node, trial_dof);
  }
};


} // namespace Impl

namespace Default = Impl::Default;

//! Customization point object on local container entry
inline const Impl::LocalContainerEntry    localContainerEntry {};


} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_LOCAL_CONTAINER_ACCESS_HH
