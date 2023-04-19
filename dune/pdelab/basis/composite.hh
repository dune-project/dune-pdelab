#ifndef DUNE_ASSEMBLER_SPACE_COMPOSITE_HH
#define DUNE_ASSEMBLER_SPACE_COMPOSITE_HH

#include <dune/assembler/space/concept.hh>
#include <dune/assembler/space/base.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/dynamicpowernode.hh>
#include <dune/typetree/powernode.hh>

#include <dune/common/indices.hh>

#include <array>
#include <memory>
#include <tuple>
#include <vector>

namespace Dune::Assembler {

/**
 * @brief Array node of Space
 *
 * @tparam MergingStrategy  Merging strategy
 * @tparam Node             A Space node type
 * @tparam degree           Size of the array
 */
template<class MergingStrategy,
         Concept::Impl::SpaceNode Node,
         std::size_t degree>
class ArraySpace
  : public TypeTree::PowerNode<Node, degree>
  , public SpaceNode<SpaceNodeTraits<MergingStrategy>>
{
private:
  using TreeNode = TypeTree::PowerNode<Node, degree>;
  using BaseNode = SpaceNode<SpaceNodeTraits<MergingStrategy>>;

public:
  using typename BaseNode::Traits;

  ArraySpace(
    const MergingStrategy& merging_strategy,
    const std::array<std::shared_ptr<Node>, degree>& nodes)
    : TreeNode{ nodes }
    , BaseNode{ merging_strategy }
  {
  }

  ArraySpace(const ArraySpace&) = default;
};

/**
 * @brief Make an array node of Space nodes
 *
 * @param merging_strategy  Merging strategy
 * @param nodes             Array of Space
 * @return auto             ArraySpace
 */
template<class MergingStrategy, class Node, std::size_t degree>
auto
makeCompositeSpace(const MergingStrategy& merging_strategy, const std::array<Node, degree>& nodes)
{
  std::array<std::shared_ptr<Node>, degree> storage;
  for (std::size_t i = 0; i < degree; ++i)
    storage[i] = std::make_shared<Node>(nodes[i]);
  using DFS = ArraySpace<MergingStrategy, Node, degree>;
  return DFS{ merging_strategy, storage };
}

/**
 * @brief Vector node of Space
 *
 * @tparam MergingStrategy  Merging strategy
 * @tparam Node             A Space node type
 */
template<class MergingStrategy, Concept::Impl::SpaceNode Node>
class VectorSpace
  : public TypeTree::DynamicPowerNode<Node>
  , public SpaceNode<SpaceNodeTraits<MergingStrategy>>
{
private:
  using TreeNode = TypeTree::DynamicPowerNode<Node>;
  using BaseNode = SpaceNode<SpaceNodeTraits<MergingStrategy>>;

public:
  using typename BaseNode::Traits;

  VectorSpace(const MergingStrategy& merging_strategy, const std::vector<std::shared_ptr<Node>>& nodes)
    : TreeNode{ nodes }
    , BaseNode{ merging_strategy }
  {
  }

  VectorSpace(const VectorSpace&) = default;
};

/**
 * @brief Make an vector node of Space nodes
 *
 * @param merging_strategy  Merging strategy
 * @param nodes             Vector of Space
 * @return auto             VectorSpace
 */
template<class MergingStrategy, class Node>
auto
makeCompositeSpace(const MergingStrategy& merging_strategy, const std::vector<Node>& nodes)
{
  std::vector<std::shared_ptr<Node>> storage(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i)
    storage[i] = std::make_shared<Node>(nodes[i]);
  using DFS = VectorSpace<MergingStrategy, Node>;
  return DFS{ merging_strategy, storage };
}

/**
 * @brief Tuple node of Space
 *
 * @tparam MergingStrategy  Merging strategy
 * @tparam Nodes            Space node types
 */
template<class MergingStrategy, Concept::Impl::SpaceNode... Nodes>
class TupleSpace
  : public TypeTree::CompositeNode<Nodes...>
  , public SpaceNode<SpaceNodeTraits<MergingStrategy>>
{
private:
  using TreeNode = TypeTree::CompositeNode<Nodes...>;
  using BaseNode = SpaceNode<SpaceNodeTraits<MergingStrategy>>;

public:
  using typename BaseNode::Traits;

  TupleSpace(const MergingStrategy& merging_strategy, const std::tuple<std::shared_ptr<Nodes>...>& nodes)
    : TreeNode{ nodes }
    , BaseNode{ merging_strategy }
  {
  }

  TupleSpace(const TupleSpace&) = default;
};

/**
 * @brief Make a tuple node of Space nodes
 *
 * @param merging_strategy  Merging strategy
 * @param nodes             Tuple of Space
 * @return auto             TupleSpace
 */
template<class MergingStrategy, class... Nodes>
auto
makeCompositeSpace(const MergingStrategy& merging_strategy, const std::tuple<Nodes...>& nodes)
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
  using DFS = TupleSpace<MergingStrategy, Nodes...>;
  return DFS{ merging_strategy, storage };
}

} // namespace Dune::Assembler

#endif // DUNE_ASSEMBLER_SPACE_COMPOSITE_HH
