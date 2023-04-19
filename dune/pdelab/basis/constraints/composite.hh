#ifndef DUNE_ASSEMBLER_SPACE_CONSTRAINTS_COMPOSITE_HH
#define DUNE_ASSEMBLER_SPACE_CONSTRAINTS_COMPOSITE_HH

#include <dune/assembler/concepts/tree.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/dynamicpowernode.hh>
#include <dune/typetree/powernode.hh>


#include <array>
#include <memory>
#include <tuple>
#include <vector>

namespace Dune::Assembler {

template<Concept::Tree Node, std::size_t degree>
struct ArrayConstraints
  : public TypeTree::PowerNode<Node, degree>
{
  ArrayConstraints(const std::array<std::shared_ptr<Node>, degree>& nodes)
    : TypeTree::PowerNode<Node, degree>{ nodes }
  {}

  ArrayConstraints(const ArrayConstraints&) = default;
};

template<Concept::Tree Node, std::size_t degree>
auto makeCompositeConstraints(const std::array<Node, degree>& nodes)
{
  std::array<std::shared_ptr<Node>, degree> storage;
  for (std::size_t i = 0; i < degree; ++i)
    storage[i] = std::make_shared<Node>(nodes[i]);
  return ArrayConstraints<Node, degree>{ storage };
}

template<Concept::Tree Node>
struct VectorConstraints
  : public TypeTree::DynamicPowerNode<Node>
{
  VectorConstraints(const std::vector<std::shared_ptr<Node>>& nodes)
    : TypeTree::DynamicPowerNode<Node>{ nodes }
  {}

  VectorConstraints(const VectorConstraints&) = default;
};

template<Concept::Tree Node>
auto makeCompositeConstraints(const std::vector<Node>& nodes)
{
  std::vector<std::shared_ptr<Node>> storage(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i)
    storage[i] = std::make_shared<Node>(nodes[i]);
  return VectorConstraints<Node>{ storage };
}

template<Concept::Tree... Nodes>
class TupleConstraints
  : public TypeTree::CompositeNode<Nodes...>
{
  TupleConstraints(const std::tuple<std::shared_ptr<Nodes>...>& nodes)
    : TypeTree::CompositeNode<Nodes...>{ nodes }
  {}

  TupleConstraints(const TupleConstraints&) = default;
};

template<Concept::Tree... Nodes>
auto makeCompositeConstraints(const std::tuple<Nodes...>& nodes)
{
  using TypeTuple = std::tuple<Nodes...>;
  using Storage = std::tuple<std::shared_ptr<Nodes>...>;
  auto sequence = std::make_index_sequence<sizeof...(Nodes)>{};
  auto storage = unpackIntegerSequence(
    [&](auto... i) {
      return Storage{ std::make_shared<std::tuple_element_t<i, TypeTuple>>(
        std::get<i>(nodes))... };
    },
    sequence);
  return TupleConstraints<Nodes...>{ storage };
}

// namespace Impl {

// template<Concept::Tree Tree, Concept::MultiIndex Path, class Callable>
// Concept::Tree auto makeConstraintsTree(const Tree& tree, Path path, Callable callable)
// {
//   if constexpr (Concept::LeafTreeNode<Tree>) {
//     static_assert(std::invocable<Callable, Tree, Path>);
//     return callable(tree, path);
//   } else if constexpr (Concept::ArrayTreeNode<Tree>) {
//     using ChidlNode = std::decay_t<decltype(makeConstraintsTree(tree.child(0),push_back(path,0), callable))>;
//     using Node = ArrayConstraintsContainer<ChidlNode, Tree::degree()>;
//     std::array<ChidlNode, Tree::degree()> array_constraints;
//     for (std::size_t i = 0; i < tree.degree(); ++i)
//       array_constraints[i] = makeConstraintsTree(tree.child(i), push_back(path,i), callable);
//     return makeCompositeConstraints(array_constraints);
//   } else if constexpr (Concept::VectorTreeNode<Tree>) {
//     using ChidlNode = std::decay_t<decltype(makeConstraintsTree(tree.child(0), push_back(path,0), callable))>;
//     std::vector<ChidlNode> vector_constraints(tree.degree());
//     for (std::size_t i = 0; i < tree.degree(); ++i)
//       vector_constraints[i] = makeConstraintsTree(tree.child(i), push_back(path,i), callable);
//     return makeCompositeConstraints(vector_constraints);
//   } else {
//     static_assert(Concept::TupleTreeNode<Tree>);
//     return unpackIntegerSequence(
//       [&](auto... i) {
//         return makeCompositeConstraints(std::make_tuple(makeConstraintsTree(tree.child(i), push_back(path,i), callable)...));
//       },
//       std::make_index_sequence<Tree::degree()>{});
//   }
// }

// }

} // namespace Dune::Assembler

#endif // DUNE_ASSEMBLER_SPACE_CONSTRAINTS_COMPOSITE_HH
