#ifndef DUNE_PDELAB_BASIS_PREBASIS_COMPOSITE_HH
#define DUNE_PDELAB_BASIS_PREBASIS_COMPOSITE_HH

#include <dune/pdelab/basis/prebasis/concept.hh>
#include <dune/pdelab/basis/prebasis/node.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/dynamicpowernode.hh>
#include <dune/typetree/powernode.hh>

#include <dune/common/indices.hh>

#include <array>
#include <memory>
#include <tuple>
#include <vector>

namespace Dune::PDELab::inline Experimental {

/**
 * @brief Array node of Space
 *
 * @tparam MergingStrategy  Merging strategy
 * @tparam Node             A Space node type
 * @tparam degree           Size of the array
 */
template<class MergingStrategy,
         Concept::Impl::PreBasisNode Node,
         std::size_t degree>
class PreBasisArray
  : public TypeTree::PowerNode<Node, degree>
  , public PreBasisNode<PreBasisNodeTraits<MergingStrategy>>
{
private:
  using TreeNode = TypeTree::PowerNode<Node, degree>;
  using BaseNode = PreBasisNode<PreBasisNodeTraits<MergingStrategy>>;

public:
  using typename BaseNode::Traits;

  template<std::same_as<void> = void>
  auto makeOrdering() const {
    return BaseNode::mergingStrategy().makeOrdering(*this);
  }

  template<std::same_as<void> = void>
  auto makeLocalOrdering() const {
    return BaseNode::mergingStrategy().makeLocalOrdering(*this);
  }

  PreBasisArray(
    const MergingStrategy& merging_strategy,
    const std::array<std::shared_ptr<Node>, degree>& nodes)
    : TreeNode{ nodes }
    , BaseNode{ merging_strategy }
  {
  }

  PreBasisArray(const PreBasisArray&) = default;
};

/**
 * @brief Make an array node of Space nodes
 *
 * @param merging_strategy  Merging strategy
 * @param nodes             Array of Space
 * @return auto             PreBasisArray
 */
template<class MergingStrategy, Concept::Impl::PreBasisNode Node, std::size_t degree>
auto
composite(const MergingStrategy& merging_strategy, const std::array<Node, degree>& nodes)
{
  std::array<std::shared_ptr<Node>, degree> storage;
  for (std::size_t i = 0; i < degree; ++i)
    storage[i] = std::make_shared<Node>(nodes[i]);
  using DFS = PreBasisArray<MergingStrategy, Node, degree>;
  return DFS{ merging_strategy, storage };
}

/**
 * @brief Vector node of Space
 *
 * @tparam MergingStrategy  Merging strategy
 * @tparam Node             A Space node type
 */
template<class MergingStrategy, Concept::Impl::PreBasisNode Node>
class PreBasisVector
  : public TypeTree::DynamicPowerNode<Node>
  , public PreBasisNode<PreBasisNodeTraits<MergingStrategy>>
{
private:
  using TreeNode = TypeTree::DynamicPowerNode<Node>;
  using BaseNode = PreBasisNode<PreBasisNodeTraits<MergingStrategy>>;

public:
  using typename BaseNode::Traits;

  template<std::same_as<void> = void>
  auto makeOrdering() const {
    return BaseNode::mergingStrategy().makeOrdering(*this);
  }

  template<std::same_as<void> = void>
  auto makeLocalOrdering() const {
    return BaseNode::mergingStrategy().makeLocalOrdering(*this);
  }

  PreBasisVector(const MergingStrategy& merging_strategy, const std::vector<std::shared_ptr<Node>>& nodes)
    : TreeNode{ nodes }
    , BaseNode{ merging_strategy }
  {
  }

  PreBasisVector(const PreBasisVector&) = default;
};

/**
 * @brief Make an vector node of Space nodes
 *
 * @param merging_strategy  Merging strategy
 * @param nodes             Vector of Space
 * @return auto             PreBasisVector
 */
template<class MergingStrategy, Concept::Impl::PreBasisNode Node>
auto
composite(const MergingStrategy& merging_strategy, const std::vector<Node>& nodes)
{
  std::vector<std::shared_ptr<Node>> storage(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i)
    storage[i] = std::make_shared<Node>(nodes[i]);
  using DFS = PreBasisVector<MergingStrategy, Node>;
  return DFS{ merging_strategy, storage };
}

/**
 * @brief Tuple node of Space
 *
 * @tparam MergingStrategy  Merging strategy
 * @tparam Nodes            Space node types
 */
template<class MergingStrategy, Concept::Impl::PreBasisNode... Nodes>
class PreBasisTuple
  : public TypeTree::CompositeNode<Nodes...>
  , public PreBasisNode<PreBasisNodeTraits<MergingStrategy>>
{
private:
  using TreeNode = TypeTree::CompositeNode<Nodes...>;
  using BaseNode = PreBasisNode<PreBasisNodeTraits<MergingStrategy>>;

public:
  using typename BaseNode::Traits;

  auto makeOrdering() const {
    return BaseNode::mergingStrategy().makeOrdering(*this);
  }

  PreBasisTuple(const MergingStrategy& merging_strategy, const std::tuple<std::shared_ptr<Nodes>...>& nodes)
    : TreeNode{ nodes }
    , BaseNode{ merging_strategy }
  {
  }

  PreBasisTuple(const PreBasisTuple&) = default;
};

/**
 * @brief Make a tuple node of Space nodes
 *
 * @param merging_strategy  Merging strategy
 * @param nodes             Tuple of Space
 * @return auto             PreBasisTuple
 */
template<class MergingStrategy, Concept::Impl::PreBasisNode... Nodes>
auto
composite(const MergingStrategy& merging_strategy, const std::tuple<Nodes...>& nodes)
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
  using DFS = PreBasisTuple<MergingStrategy, Nodes...>;
  return DFS{ merging_strategy, storage };
}

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_BASIS_PREBASIS_COMPOSITE_HH
