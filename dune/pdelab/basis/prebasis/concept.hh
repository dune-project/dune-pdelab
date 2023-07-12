#ifndef DUNE_PDELAB_BASIS_PREBASIS_CONCEPT_HH
#define DUNE_PDELAB_BASIS_PREBASIS_CONCEPT_HH

#include <dune/pdelab/concepts/treenode.hh>
#include <dune/pdelab/common/tree_traversal.hh>

#include <concepts>

namespace Dune::PDELab::inline Experimental::Concept::Impl {

template<class Node>
concept PreBasisNode = requires(Node node)
{
  { node.mergingStrategy() } -> std::convertible_to<typename std::remove_cvref_t<Node>::Traits::MergingStrategy>;
  { node.name() } -> std::convertible_to<std::string_view>;
  node.name(std::string{});
  node.name(std::string_view{});
  requires std::copy_constructible<Node>;
  requires std::move_constructible<Node>;
};

void
requirePreBasisNode(PreBasisNode auto&& node);

template<class Leaf>
concept PreBasisLeaf = requires(Leaf leaf)
{
  requires PreBasisNode<Leaf>;
  { leaf.finiteElementMap() } -> std::convertible_to<typename std::remove_cvref_t<Leaf>::Traits::FiniteElementMap>;
  { leaf.constraintsOperator() } -> std::convertible_to<typename std::remove_cvref_t<Leaf>::Traits::ConstraintsOperator>;
  requires std::equality_comparable<Leaf>;
};

void
requirePreBasisLeaf(PreBasisLeaf auto&& leaf);

template<class PBT>
concept PreBasisTree = requires(PBT pre_basis_tree)
{
  requires TreeNode<PBT>;
  Dune::PDELab::forEachNode(
    pre_basis_tree, [](auto&& node, auto tp) {
      Impl::requirePreBasisNode(node);
      if constexpr (Concept::LeafTreeNode<decltype(node)>)
        Impl::requirePreBasisLeaf(node);
    });
};

} // namespace Dune::PDELab::inline Experimental::Concept::Impl

#endif // DUNE_PDELAB_BASIS_PREBASIS_CONCEPT_HH
