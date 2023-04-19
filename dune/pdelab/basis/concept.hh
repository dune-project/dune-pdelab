#ifndef DUNE_ASSEMBLER_SPACE_CONCEPT_HH
#define DUNE_ASSEMBLER_SPACE_CONCEPT_HH

#include <dune/assembler/concepts/treenode.hh>

#include <dune/assembler/common/tree_traversal.hh>

#include <concepts>

namespace Dune::Assembler::Concept {

namespace Impl {

template<class Node>
concept SpaceNode = requires(Node node)
{
  { node.mergingStrategy() } -> std::convertible_to<typename std::remove_cvref_t<Node>::Traits::MergingStrategy>;
  { node.name() } -> std::convertible_to<std::string>;
  node.name(std::string{});
  requires std::copy_constructible<Node>;
  requires std::move_constructible<Node>;
};

void
requireSpaceNode(SpaceNode auto&& node);

template<class Leaf>
concept SpaceLeaf = requires(Leaf leaf)
{
  requires SpaceNode<Leaf>;
  { leaf.finiteElementMap() } -> std::convertible_to<typename std::remove_cvref_t<Leaf>::Traits::FiniteElementMap>;
  { leaf.constraintsOperator() } -> std::convertible_to<typename std::remove_cvref_t<Leaf>::Traits::ConstraintsOperator>;
  requires std::equality_comparable<Leaf>;
};

void
requireSpaceLeaf(SpaceLeaf auto&& leaf);

template<class DFS>
concept SpaceTree = requires(DFS discrete_function_space)
{
  requires TreeNode<DFS>;
  Dune::Assembler::forEachNode(
    discrete_function_space, [](auto&& node, auto tp) {
      Impl::requireSpaceNode(node);
      if constexpr (Concept::LeafTreeNode<decltype(node)>)
        Impl::requireSpaceLeaf(node);
    });
};
}

} // namespace Dune::Assembler::Concept

#endif // DUNE_ASSEMBLER_SPACE_CONCEPT_HH
