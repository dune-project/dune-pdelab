#ifndef DUNE_PDELAB_COMMON_TREE_ENTRY_HH
#define DUNE_PDELAB_COMMON_TREE_ENTRY_HH

#include <dune/pdelab/concepts/treenode.hh>
#include <dune/pdelab/concepts/multiindex.hh>

#include <dune/pdelab/common/container_entry.hh>

#include <utility>
#include <type_traits>

namespace Dune::PDELab::inline Experimental {

namespace Impl {

inline namespace Default {

/**
 * @brief Access an entry in a tree container
 * @details Use front-most index to return a child from the tree
 *
 * @tparam Container        Container type
 * @tparam MultiIndex       Multi-index type
 * @param container         Tree container with the desired entry
 * @param multiindex        Multi-index to reach the target entry
 * @return decltype(auto)   Entry at multi-index
 */
template<class TreeNode, Concept::FixedSizeMultiIndex MultiIndex>
requires (MultiIndex::max_size() != 0 && Concept::ParentTreeNode<std::remove_cvref_t<TreeNode>>)
constexpr decltype(auto) containerEntry(TreeNode&& tree_node, MultiIndex multiindex) {
  auto index = front(multiindex);
  return containerEntry(tree_node.child(index), pop_front(multiindex));
}

} // namespace Default

} // namespace Impl

namespace Default = Impl::Default;

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_TREE_ENTRY_HH
