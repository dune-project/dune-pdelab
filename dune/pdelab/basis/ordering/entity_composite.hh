#ifndef DUNE_PDELAB_BASIS_ORDERING_ENTITY_COMPOSITE_HH
#define DUNE_PDELAB_BASIS_ORDERING_ENTITY_COMPOSITE_HH

#include <dune/pdelab/basis/ordering/entity_base.hh>

#include <dune/pdelab/common/multiindex.hh>
#include <dune/pdelab/common/tree_traversal.hh>

#include <dune/pdelab/concepts/multiindex.hh>
#include <dune/pdelab/concepts/treenode.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/dynamicpowernode.hh>
#include <dune/typetree/powernode.hh>
#include <dune/typetree/treepath.hh>

#include <vector>
#include <memory>

namespace Dune::PDELab::inline Experimental::Impl {

/**
 * @brief Array implementation of an EntityOrderingNode
 *
 * @tparam MergingStrategy   A merging strategy between the node indices
 * @tparam ChildOrdering     Child ordering type to store
 * @tparam degree            Number of child orderings to keep
 */
template<class MergingStrategy,
         Concept::TreeNode ChildOrdering,
         std::size_t degree>
class ArrayEntityOrdering
  : public TypeTree::PowerNode<ChildOrdering, degree>
  , public EntityOrderingNode<
      ArrayEntityOrdering<MergingStrategy, ChildOrdering, degree>,
      MergingStrategy>
{
  using TreeNode = TypeTree::PowerNode<ChildOrdering, degree>;
  using OrderingNode = EntityOrderingNode<
    ArrayEntityOrdering<MergingStrategy, ChildOrdering, degree>,
    MergingStrategy>;

  //! Array tree node of local orderings  @see Concepts::ArrayTreeNode
  template<class ChildLocalNode>
  struct LocalNode : public TypeTree::PowerNode<ChildLocalNode, degree>
  {
    LocalNode(typename TypeTree::PowerNode<ChildLocalNode,
                                         degree>::NodeStorage&& storage)
      : TypeTree::PowerNode<ChildLocalNode, degree>{ std::move(storage) }
    {
    }

    LocalNode(const LocalNode&) = delete;
    LocalNode(LocalNode&&) = default;

    LocalNode& operator=(const LocalNode&) = delete;
    LocalNode& operator=(LocalNode&&) = default;

    void bindElement(const auto& entity) noexcept {
      forEach(*this, [&](auto& child){ child.bindElement(entity); });
    }
  };

public:
  //! Constructs an array of entity orderings
  ArrayEntityOrdering(typename TreeNode::NodeStorage&& storage,
                      const MergingStrategy& merging_strategy)
    : TreeNode{ std::move(storage) }
    , OrderingNode{ merging_strategy }
  {
  }

  ArrayEntityOrdering(const ArrayEntityOrdering&) = delete;
  ArrayEntityOrdering(ArrayEntityOrdering&&) = default;

  ArrayEntityOrdering& operator=(const ArrayEntityOrdering&) = delete;
  ArrayEntityOrdering& operator=(ArrayEntityOrdering&&) = default;

  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubSpacePath>
  auto makeLocalView(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubSpacePath& sub_space_path) const
  {
    if constexpr (Prefix::size() < SubSpacePath::size()) {
      auto child = sub_space_path[index_constant<Prefix::size()>{}];
      return this->child(child).makeLocalView(ordering, push_back(prefix, child), sub_space_path);
    } else {
      using Child = std::decay_t<decltype(*this->child(0).makeLocalView(
        ordering, push_back(prefix, 0), sub_space_path))>;
      typename LocalNode<Child>::NodeStorage storage;
      for (std::size_t i = 0; i < degree; ++i)
        storage[i] =
          this->child(i).makeLocalView(ordering, push_back(prefix, i), sub_space_path);
      return std::make_unique<LocalNode<Child>>(std::move(storage));
    }
  }

  //! Constructs an array of the local entity orderings from the entity
  //! orderings
  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubSpacePath>
  auto makeLocalIndexSet(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubSpacePath& sub_space_path) const
  {
    if constexpr (Prefix::size() < SubSpacePath::size()) {
      auto child = sub_space_path[index_constant<Prefix::size()>{}];
      return this->child(child).makeLocalIndexSet(ordering, push_back(prefix, child), sub_space_path);
    } else {
      using Child = std::decay_t<decltype(*this->child(0).makeLocalIndexSet(
        ordering, push_back(prefix, 0), sub_space_path))>;
      typename LocalNode<Child>::NodeStorage storage;
      for (std::size_t i = 0; i < degree; ++i)
        storage[i] =
          this->child(i).makeLocalIndexSet(ordering, push_back(prefix, i), sub_space_path);
      return std::make_unique<LocalNode<Child>>(std::move(storage));
    }
  }

private:
  template<class Container>
  using StaticSize = decltype(Container::size());

public:
  template<class Traits>
  static constexpr auto makeVectorContainer()
  {
    using ChildContainer =
      decltype(ChildOrdering::template makeVectorContainer<Traits>());
    if constexpr (OrderingNode::containerBlocked()) {
      return Traits::template makeArray<ChildContainer, degree>();
    } else {
      using ChildBlock =
        typename Traits::template block_type<ChildContainer>::type;
      if constexpr (Std::is_detected<StaticSize, ChildContainer>{})
        return Traits::template makeArray<ChildBlock,
                                          ChildContainer::size() * degree>();
      else
        return Traits::template makeVector<ChildBlock>();
    }
  }
};

template<class MergingStrategy, Concept::TreeNode ChildOrdering>
class VectorEntityOrdering
  : public TypeTree::DynamicPowerNode<ChildOrdering>
  , public EntityOrderingNode<
      VectorEntityOrdering<MergingStrategy, ChildOrdering>,
      MergingStrategy>
{
  using TreeNode = TypeTree::DynamicPowerNode<ChildOrdering>;
  using OrderingNode =
    EntityOrderingNode<VectorEntityOrdering<MergingStrategy, ChildOrdering>,
                       MergingStrategy>;

  template<class ChildLocalNode>
  struct LocalNode : public TypeTree::DynamicPowerNode<ChildLocalNode>
  {
    LocalNode(typename TypeTree::DynamicPowerNode<
                  ChildLocalNode>::NodeStorage&& storage)
      : TypeTree::DynamicPowerNode<ChildLocalNode>{ std::move(storage) }
    {
    }

    LocalNode(const LocalNode&) = delete;
    LocalNode(LocalNode&&) = default;

    LocalNode& operator=(const LocalNode&) = delete;
    LocalNode& operator=(LocalNode&&) = default;

    void bindElement(const auto& entity) noexcept {
      forEach(*this, [&](auto& child){ child.bindElement(entity); });
    }
  };

public:
  VectorEntityOrdering(typename TreeNode::NodeStorage&& storage,
                       const MergingStrategy& merging_strategy)
    : TreeNode{ std::move(storage) }
    , OrderingNode{ merging_strategy }
  {
    for(std::size_t i = 0; i != this->degree(); ++i)
      assert(this->childStorage(i));
  }

  VectorEntityOrdering(const VectorEntityOrdering&) = delete;
  VectorEntityOrdering(VectorEntityOrdering&&) = default;

  VectorEntityOrdering& operator=(const VectorEntityOrdering&) = delete;
  VectorEntityOrdering& operator=(VectorEntityOrdering&&) = default;

  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubSpacePath>
  auto makeLocalIndexSet(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubSpacePath& sub_space_path) const
  {
    if constexpr (Prefix::size() < SubSpacePath::size()) {
      auto child = sub_space_path[index_constant<Prefix::size()>{}];
      assert(child < this->degree());
      return this->child(child).makeLocalIndexSet(ordering, push_back(prefix, child), sub_space_path);
    } else {
      using Child = std::decay_t<decltype(*this->child(0).makeLocalIndexSet(
        ordering, push_back(prefix, 0), sub_space_path))>;
      typename LocalNode<Child>::NodeStorage storage(this->degree());
      for (std::size_t i = 0; i < this->degree(); ++i)
        storage[i] =
          this->child(i).makeLocalIndexSet(ordering, push_back(prefix, i), sub_space_path);
      return std::make_unique<LocalNode<Child>>(std::move(storage));
    }
  }

  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubSpacePath>
  auto makeLocalView(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubSpacePath& sub_space_path) const
  {
    if constexpr (Prefix::size() < SubSpacePath::size()) {
      auto child = sub_space_path[index_constant<Prefix::size()>{}];
      assert(child < this->degree());
      return this->child(child).makeLocalView(ordering, push_back(prefix, child), sub_space_path);
    } else {
      using Child = std::decay_t<decltype(*this->child(0).makeLocalView(
        ordering, push_back(prefix, 0), sub_space_path))>;
      typename LocalNode<Child>::NodeStorage storage(this->degree());
      for (std::size_t i = 0; i < this->degree(); ++i)
        storage[i] =
          this->child(i).makeLocalView(ordering, push_back(prefix, i), sub_space_path);
      return std::make_unique<LocalNode<Child>>(std::move(storage));
    }
  }

  template<class Traits>
  static constexpr auto makeVectorContainer()
  {
    using ChildContainer =
      decltype(ChildOrdering::template makeVectorContainer<Traits>());
    if constexpr (OrderingNode::containerBlocked()) {
      return Traits::template makeVector<ChildContainer>();
    } else {
      using ChildBlock =
        typename Traits::template block_type<ChildContainer>::type;
      return Traits::template makeVector<ChildBlock>();
    }
  }
};

template<class MergingStrategy, Concept::TreeNode... ChildOrdering>
class TupleEntityOrdering
  : public TypeTree::CompositeNode<ChildOrdering...>
  , public EntityOrderingNode<
      TupleEntityOrdering<MergingStrategy, ChildOrdering...>,
      MergingStrategy>
{
  using TreeNode = TypeTree::CompositeNode<ChildOrdering...>;
  using OrderingNode =
    EntityOrderingNode<TupleEntityOrdering<MergingStrategy, ChildOrdering...>,
                       MergingStrategy>;

  template<class... ChildLocalNode>
  struct LocalNode : public TypeTree::CompositeNode<ChildLocalNode...>
  {
    LocalNode(typename TypeTree::CompositeNode<
                  ChildLocalNode...>::NodeStorage&& storage)
      : TypeTree::CompositeNode<ChildLocalNode...>{ std::move(storage) }
    {
    }

    LocalNode(const LocalNode&) = delete;
    LocalNode(LocalNode&&) = default;

    LocalNode& operator=(const LocalNode&) = delete;
    LocalNode& operator=(LocalNode&&) = default;

    void bindElement(const auto& entity) noexcept {
      forEach(*this, [&](auto& child){ child.bindElement(entity); });
    }
  };

public:
  TupleEntityOrdering(typename TreeNode::NodeStorage&& storage,
                      const MergingStrategy& merging_strategy)
    : TreeNode{ std::move(storage) }
    , OrderingNode{ merging_strategy }
  {
  }

  TupleEntityOrdering(const TupleEntityOrdering&) = delete;
  TupleEntityOrdering(TupleEntityOrdering&&) = default;

  TupleEntityOrdering& operator=(const TupleEntityOrdering&) = delete;
  TupleEntityOrdering& operator=(TupleEntityOrdering&&) = default;

  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubSpacePath>
  auto makeLocalIndexSet(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubSpacePath& sub_space_path) const
  {
    if constexpr (Prefix::size() < SubSpacePath::size()) {
      auto child = sub_space_path[index_constant<Prefix::size()>{}];
      return this->child(child).makeLocalIndexSet(ordering, push_back(prefix, child), sub_space_path);
    } else {
      auto unfold_children = [&](auto... i) {
        using LSpace =
          LocalNode<std::decay_t<decltype(*this->child(i).makeLocalIndexSet(
            ordering, push_back(prefix, i), sub_space_path))>...>;
        typename LSpace::NodeStorage storage{ this->child(i).makeLocalIndexSet(
          ordering, push_back(prefix, i), sub_space_path)... };
        return std::make_unique<LSpace>(std::move(storage));
      };
      auto indices = std::make_index_sequence<sizeof...(ChildOrdering)>{};
      return unpackIntegerSequence(unfold_children, indices);
    }
  }

  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubSpacePath>
  auto makeLocalView(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubSpacePath& sub_space_path) const
  {
    if constexpr (Prefix::size() < SubSpacePath::size()) {
      auto child = sub_space_path[index_constant<Prefix::size()>{}];
      return this->child(child).makeLocalView(ordering, push_back(prefix, child), sub_space_path);
    } else {
      auto unfold_children = [&](auto... i) {
        using LSpace =
          LocalNode<std::decay_t<decltype(*this->child(i).makeLocalView(
            ordering, push_back(prefix, i), sub_space_path))>...>;
        typename LSpace::NodeStorage storage{ this->child(i).makeLocalView(
          ordering, push_back(prefix, i), sub_space_path)... };
        return std::make_unique<LSpace>(std::move(storage));
      };
      auto indices = std::make_index_sequence<sizeof...(ChildOrdering)>{};
      return unpackIntegerSequence(unfold_children, indices);
    }
  }

private:
  template<class... C>
  using StaticSizes = decltype((C::size() + ...));

  template<class... T>
  using HasCommonType = std::common_type_t<T...>;

public:
  template<class Traits>
  static constexpr auto makeVectorContainer()
  {
    if constexpr (OrderingNode::containerBlocked()) {
      return Traits::template makeTuple<
        decltype(ChildOrdering::template makeVectorContainer<Traits>())...>();
    } else {
      using ChildContainers = std::tuple<
        decltype(ChildOrdering::template makeVectorContainer<Traits>())...>;
      return unpackIntegerSequence(
        [](auto... i) {
          static_assert(
            Std::is_detected<
              HasCommonType,
              typename Traits::template block_type<
                std::tuple_element_t<i, ChildContainers>>::type...>{},
            "Non-Blocked tuple children should yield a common container type");
          using ChildBlock =
            std::common_type_t<typename Traits::template block_type<
              std::tuple_element_t<i, ChildContainers>>::type...>;
          if constexpr (Std::is_detected<
                          StaticSizes,
                          std::tuple_element_t<i, ChildContainers>...>{}) {
            constexpr std::size_t aggregated_size =
              (std::tuple_element_t<i, ChildContainers>::size() + ...);
            return Traits::template makeArray<ChildBlock, aggregated_size>();
          } else
            return Traits::template makeVector<ChildBlock>();
        },
        std::index_sequence_for<ChildOrdering...>{});
    }
  }
};

} // namespace Dune::PDELab::inline Experimental::Impl

#endif // DUNE_PDELAB_BASIS_ORDERING_ENTITY_COMPOSITE_HH
