#ifndef DUNE_PDELAB_BASIS_ORDERING_LEXICOGRAPHIC_HH
#define DUNE_PDELAB_BASIS_ORDERING_LEXICOGRAPHIC_HH

#include <dune/pdelab/basis/prebasis/concept.hh>
#include <dune/pdelab/basis/merging_strategy.hh>

#include <dune/pdelab/common/multiindex.hh>
#include <dune/pdelab/common/tree_traversal.hh>

namespace Dune::PDELab::inline Experimental::Impl {

template<class Node, class MS>
class LexicographicOrderingNode
{
  using SizeType = typename MS::SizeType;

public:
  LexicographicOrderingNode() { update(); }

  LexicographicOrderingNode(const LexicographicOrderingNode&) = delete;
  LexicographicOrderingNode(LexicographicOrderingNode&&) = delete;

  static constexpr auto containerBlocked()
  {
    return std::integral_constant<bool, MS::Blocked>{};
  }

  // Gives the maximum size of a prefix produced by this ordering
  [[nodiscard]] static constexpr std::size_t maxContainerDepth()
  {
    static_assert(not Concept::LeafTreeNode<Node>);
    auto child_depth = [&]() {
      if constexpr (Concept::ArrayTreeNode<Node> ||
                    Concept::VectorTreeNode<Node>) {
        return Node::ChildType::maxContainerDepth();
      } else {
        static_assert(Node::isComposite);
        return unpackIntegerSequence(
          [](auto... i) {
            return std::max({
              TypeTree::template Child<Node, i>::maxContainerDepth()...});
          },
          std::make_index_sequence<Node::degree()>{});
      }
    }();
    if constexpr (containerBlocked())
      return child_depth + 1;
    else
      return child_depth;
  }

  [[nodiscard]] static constexpr auto prioryFixedSize()
  {
    static_assert(not Node::isLeaf);
    if constexpr (Node::isPower)
      return Node::ChildType::prioryFixedSize();
    else {
      static_assert(Node::isComposite);
      auto unfold_children = [&](auto... i) {
        constexpr bool all_fixed_size =
          (TypeTree::template Child<Node, i>::prioryFixedSize() && ...);
        return std::integral_constant<bool, all_fixed_size>{};
      };
      auto indices = std::make_index_sequence<Node::degree()>{};
      return unpackIntegerSequence(unfold_children, indices);
    }
  }

  [[nodiscard]] auto fixedSize() const
  {
    if constexpr (prioryFixedSize())
      return std::true_type{};
    else
      return _fixed_size;
  }

  [[nodiscard]] SizeType maxLocalCount() const
  {
    SizeType max_ls = 0;
    forEach(node(), [&](const auto& child) {
      using std::max;
      max_ls = max(max_ls, child.maxLocalCount());
    });
    return max_ls;
  }

  [[nodiscard]] SizeType dimension() const
  {
    SizeType dof_count = 0;
    forEach(node(), [&](const auto& child) { dof_count += child.dimension(); });
    return dof_count;
  }

  [[nodiscard]] auto blockCount() const
  {
    if constexpr (containerBlocked())
      return node().degree();
    else
      return _child_block_offsets.back();
  }

  // Check if a given codimension is mapped to multi-indices
  [[nodiscard]] bool containsCodim(SizeType codim) const noexcept
  {
    bool contained = false;
    forEach(node(), [&](const auto& child) { contained |= child.containsCodim(codim); });
    return contained;
  }

  template<Concept::MultiIndex CompositionSuffix>
  [[nodiscard]] Concept::MultiIndex auto firstContainerIndex(
    CompositionSuffix comp_suf,
    SizeType gt_index,
    SizeType entity_index) const noexcept
  {
    // Note: Multi-index is read Outer->Inner
    const auto child_index = front(comp_suf);
    // get container suffix
    auto ci = node()
                .child(child_index)
                .firstContainerIndex(pop_front(comp_suf), gt_index, entity_index);

    if constexpr (containerBlocked())
      return push_front(ci, child_index);
    else
      return accumulate_front(ci, _child_block_offsets[child_index]);
  }

  template<Concept::MultiIndex ContainerSuffix>
  std::size_t containerSize(const ContainerSuffix& cs) const
  {
    // Note: Multi-index is read Inner->Outer
    // suffix wants the size for this node
    if (cs.size() == 0)
      return node().blockCount();

    // helper to return from any child with a dynamic child index
    auto childContainerSize = [&](std::size_t child_i,
                                  auto next_suffix) -> SizeType {
      if constexpr (Concept::ArrayTreeNode<Node> ||
                    Concept::VectorTreeNode<Node>) {
        return node().child(child_i).containerSize(next_suffix);
      } else {
        static_assert(Concept::TupleTreeNode<Node>);
        // at this point we recoverd the index, but there is no way to
        // propagate its static information outside of this function (i.e. a
        // return type that depends on the child index)
        SizeType _size{0};
        // make a loop over all nodes and check which one matches the child
        // index
        forEach(node(), [&](auto& child, auto i) {
          if (i == child_i)
            _size = child.containerSize(next_suffix);
        });
        return _size;
      }
    };

    return applyToChild(cs, childContainerSize);
  }

  // gets a handle to lock an specific degree of freedom
  Concept::Lockable auto dofLockHandle(Concept::MultiIndex auto suffix) {
    // Note: Multi-index is read Inner->Outer

    // helper to return from any child with a dynamic child index
    auto childLockHandle = [&](std::size_t child_i,
                                  auto next_suffix) -> SizeType {
      if constexpr (Concept::ArrayTreeNode<Node> ||
                    Concept::VectorTreeNode<Node>) {
        return node().child(child_i).dofLockHandle(next_suffix);
      } else {
        static_assert(Concept::TupleTreeNode<Node>);
        // at this point we recoverd the index, but there is no way to
        // propagate its static information outside of this function (i.e. a
        // return type that depends on the child index)
        using Lock = decltype(node().child(Indices::_0).dofLockHandle(next_suffix));
        Lock _lock;
        // make a loop over all nodes and check which one matches the child
        // index
        forEach(node(), [&](auto& child, auto i) {
          using LockI = decltype(node().child(i).dofLockHandle(next_suffix));
          static_assert(std::assignable_from<Lock,LockI>, "Locks from ordering nodes are different!");
          if (i == child_i)
            _lock = child.dofLockHandle(next_suffix);
        });
        return _lock;
      }
    };

    return applyToChild(suffix, childLockHandle);
  }


  void update()
  {

    SizeType block_carry = 0;
    _fixed_size = prioryFixedSize();

    if constexpr (not containerBlocked()) {
      _child_block_offsets.resize(node().degree() + 1);
      std::fill(_child_block_offsets.begin(), _child_block_offsets.end(), 0);
    }

    forEach(node(), [&](auto& child, auto i) {
      child.update();

      if constexpr (not prioryFixedSize())
        _fixed_size &= child.fixedSize();

      block_carry += child.blockCount();

      if constexpr (not containerBlocked())
        _child_block_offsets[i+1] = block_carry;
    });
  }

private:

  template<Concept::MultiIndex ContainerSuffix>
  auto applyToChild(const ContainerSuffix& cs, auto fapply) const
  {
    // Note: Multi-index is read Inner->Outer
    // transform to reserved multi-index to avoid problems on pop back
    auto rcs = Dune::PDELab::MultiIndex<SizeType, maxContainerDepth()>{ cs };

    // the next index to find out its size
    auto back_index = back(rcs);
    // task: find child the child node for whom this index corresponds
    if constexpr (containerBlocked()) {
      // easy case, the back_index is exactly the index of the child node
      return fapply(back_index, pop_back(rcs));
    } else {
      // here we need to "recover" the child index that describes the
      // back_index (invers of firstContainerIndex operation)
      auto dof_begin = _child_block_offsets.begin();
      auto dof_end = _child_block_offsets.end();
      auto dof_it = std::upper_bound(dof_begin, dof_end, back_index);
      auto next = accumulate_back(rcs, SizeType{ 0 });
      if (dof_it != dof_begin) {
        std::advance(dof_it, -1);
        assert(back(cs) >= *dof_it);
        next = accumulate_back(next, -(*dof_it));
      }
      std::size_t child_index = std::distance(dof_begin, dof_it);
      assert(node().degree() > child_index);
      return fapply(child_index, next);
    }
  }


  //! Cast to node implementation (Barton–Nackman trick)
  const Node& node() const { return static_cast<const Node&>(*this); }
  Node& node() { return static_cast<Node&>(*this); }

  bool _fixed_size;
  std::vector<SizeType> _child_block_offsets;
};

template<class Node>
struct LexicographicLocalNodeBase {
  template<class T>
  void bind(const T& t) {
    forEach(node(),[&](auto& child){ child.bind(t); });
  }

  void unbind() {
    forEach(node(),[](auto& child){ child.unbind(); });
  }

  void lock() noexcept {
    forEach(node(),[](auto& child){ child.lock(); });
  };

  [[nodiscard]] bool try_lock() noexcept {
    static_assert(not Concept::LeafTreeNode<Node>);
    using namespace Dune::Indices;

    bool succed = true;
    forEach(node(), [&](auto& child, auto i){
      // try to lock every child on the node
      if (not child.try_lock()) {
        // already locked, we have to roll back
        Dune::Hybrid::forEach(Dune::range(i), [&](auto j){
          node().child( Dune::Hybrid::minus(i,j) ).unlock();
        });
        // ...and inform that we could not adquire the lock
        succed = false;
      }
    });
    return succed;
  };

  void unlock() noexcept  {
    forEach(node(),[](auto& child){ child.unlock(); });
  };
private:
  //! Cast to node implementation (Barton–Nackman trick)
  const Node& node() const { return static_cast<const Node&>(*this); }
  Node& node() { return static_cast<Node&>(*this); }
};

template<class MergingStrategy, class ChildOrdering, std::size_t degree>
class ArrayLexicographicOrdering
  : public TypeTree::PowerNode<ChildOrdering, degree>
  , public LexicographicOrderingNode<
      ArrayLexicographicOrdering<MergingStrategy, ChildOrdering, degree>,
      MergingStrategy>
{
  using TreeNode = TypeTree::PowerNode<ChildOrdering, degree>;
  using OrderingNode = LexicographicOrderingNode<
    ArrayLexicographicOrdering<MergingStrategy, ChildOrdering, degree>,
    MergingStrategy>;

  template<class ChildLocalNode>
  struct LocalNode
   : public TypeTree::PowerNode<ChildLocalNode, degree>
   , public LexicographicLocalNodeBase<LocalNode<ChildLocalNode>>
  {
    LocalNode(typename TypeTree::PowerNode<ChildLocalNode,degree>::NodeStorage&& storage)
      : TypeTree::PowerNode<ChildLocalNode, degree>{ std::move(storage) }
    {
    }
  };

public:
  ArrayLexicographicOrdering(typename TreeNode::NodeStorage&& storage)
    : TreeNode{ std::move(storage) }
    , OrderingNode{}
  {
  }

  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubPreBasisPath>
  auto makeLocalIndexSet(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubPreBasisPath& sub_pre_basis_path) const
  {
    if constexpr (Prefix::size() < SubPreBasisPath::size()) {
      auto child = sub_pre_basis_path[index_constant<Prefix::size()>{}];
      return this->child(child).makeLocalIndexSet(ordering, push_back(prefix, child), sub_pre_basis_path);
    } else {
      using Child = std::decay_t<decltype(*this->child(0).makeLocalIndexSet(ordering, push_back(prefix, 0), sub_pre_basis_path))>;
      typename LocalNode<Child>::NodeStorage storage;
      for (std::size_t i = 0; i < degree; ++i)
        storage[i] = this->child(i).makeLocalIndexSet(ordering, push_back(prefix, i), sub_pre_basis_path);
      return std::make_unique<LocalNode<Child>>(std::move(storage));
    }
  }

  //! Constructs an array of the local entity orderings from the entity
  //! orderings
  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubPreBasisPath>
  auto makeLocalView(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubPreBasisPath& sub_pre_basis_path) const
  {
    if constexpr (Prefix::size() < SubPreBasisPath::size()) {
      std::size_t child = sub_pre_basis_path[index_constant<Prefix::size()>{}];
      assert(child < this->degree());
      return this->child(child).makeLocalView(ordering, push_back(prefix, child), sub_pre_basis_path);
    } else {
      using Child = std::decay_t<decltype(*this->child(0).makeLocalView(ordering, push_back(prefix, 0), sub_pre_basis_path))>;
      typename LocalNode<Child>::NodeStorage storage;
      for (std::size_t i = 0; i < degree; ++i)
        storage[i] = this->child(i).makeLocalView(ordering, push_back(prefix, i), sub_pre_basis_path);
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

template<class MergingStrategy, class ChildOrdering>
class VectorLexicographicOrdering
  : public TypeTree::DynamicPowerNode<ChildOrdering>
  , public LexicographicOrderingNode<
      VectorLexicographicOrdering<MergingStrategy, ChildOrdering>,
      MergingStrategy>
{
  using TreeNode = TypeTree::DynamicPowerNode<ChildOrdering>;
  using OrderingNode = LexicographicOrderingNode<
    VectorLexicographicOrdering<MergingStrategy, ChildOrdering>,
    MergingStrategy>;

  template<class ChildLocalNode>
  struct LocalNode
  : public TypeTree::DynamicPowerNode<ChildLocalNode>
  , public LexicographicLocalNodeBase<LocalNode<ChildLocalNode>>
  {
    LocalNode(typename TypeTree::DynamicPowerNode<ChildLocalNode>::NodeStorage&& storage)
      : TypeTree::DynamicPowerNode<ChildLocalNode>{ std::move(storage) }
    {}
  };

public:
  VectorLexicographicOrdering(typename TreeNode::NodeStorage&& storage)
    : TreeNode{ std::move(storage) }
    , OrderingNode{}
  {
  }

  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubPreBasisPath>
  auto makeLocalIndexSet(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubPreBasisPath& sub_pre_basis_path) const
  {
    if constexpr (Prefix::size() < SubPreBasisPath::size()) {
      auto child = sub_pre_basis_path[index_constant<Prefix::size()>{}];
      assert(child < this->degree());
      return this->child(child).makeLocalIndexSet(ordering, push_back(prefix, child), sub_pre_basis_path);
    } else {
      using Child = std::decay_t<decltype(*this->child(0).makeLocalIndexSet(
        ordering, push_back(prefix, 0), sub_pre_basis_path))>;
      typename LocalNode<Child>::NodeStorage storage(this->degree());
      for (std::size_t i = 0; i < this->degree(); ++i)
        storage[i] =
          this->child(i).makeLocalIndexSet(ordering, push_back(prefix, i), sub_pre_basis_path);
      return std::make_unique<LocalNode<Child>>(std::move(storage));
    }
  }

  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubPreBasisPath>
  auto makeLocalView(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubPreBasisPath& sub_pre_basis_path) const
  {
    if constexpr (Prefix::size() < SubPreBasisPath::size()) {
      auto child = sub_pre_basis_path[index_constant<Prefix::size()>{}];
      assert(child < this->degree());
      return this->child(child).makeLocalView(ordering, push_back(prefix, child), sub_pre_basis_path);
    } else {
      using Child = std::decay_t<decltype(*this->child(0).makeLocalView(
        ordering, push_back(prefix, 0), sub_pre_basis_path))>;
      typename LocalNode<Child>::NodeStorage storage(this->degree());
      for (std::size_t i = 0; i < this->degree(); ++i)
        storage[i] = this->child(i).makeLocalView(ordering, push_back(prefix, i), sub_pre_basis_path);
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

template<class MergingStrategy, class... ChildOrdering>
class TupleLexicographicOrdering
  : public TypeTree::CompositeNode<ChildOrdering...>
  , public LexicographicOrderingNode<
      TupleLexicographicOrdering<MergingStrategy, ChildOrdering...>,
      MergingStrategy>
{
  static_assert(sizeof...(ChildOrdering) > 0);
  using TreeNode = TypeTree::CompositeNode<ChildOrdering...>;
  using OrderingNode = LexicographicOrderingNode<
    TupleLexicographicOrdering<MergingStrategy, ChildOrdering...>,
    MergingStrategy>;

  template<class... ChildLocalNode>
  struct LocalNode
    : public TypeTree::CompositeNode<ChildLocalNode...>
    , public LexicographicLocalNodeBase<LocalNode<ChildLocalNode...>>
  {
    LocalNode(typename TypeTree::CompositeNode<ChildLocalNode...>::NodeStorage&& storage)
      : TypeTree::CompositeNode<ChildLocalNode...>{ std::move(storage) }
    {}
  };

public:
  TupleLexicographicOrdering(typename TreeNode::NodeStorage&& storage)
    : TreeNode{ std::move(storage) }
    , OrderingNode{}
  {}

  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubPreBasisPath>
  auto makeLocalIndexSet(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubPreBasisPath& sub_pre_basis_path) const
  {
    if constexpr (Prefix::size() < SubPreBasisPath::size()) {
      auto child = sub_pre_basis_path[index_constant<Prefix::size()>{}];
      static_assert(child < this->degree());
      return this->child(child).makeLocalIndexSet(ordering, push_back(prefix, child), sub_pre_basis_path);
    } else {
      auto unfold_children = [&](auto... i) {
        using LPreBasis = LocalNode<std::decay_t<decltype(*this->child(i).makeLocalIndexSet(ordering, push_back(prefix, i), sub_pre_basis_path))>...>;
        typename LPreBasis::NodeStorage storage{ this->child(i).makeLocalIndexSet(ordering, push_back(prefix, i), sub_pre_basis_path)... };
        return std::make_unique<LPreBasis>(std::move(storage));
      };
      auto indices = std::make_index_sequence<sizeof...(ChildOrdering)>{};
      return unpackIntegerSequence(unfold_children, indices);
    }
  }

  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubPreBasisPath>
  auto makeLocalView(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubPreBasisPath& sub_pre_basis_path) const
  {
    if constexpr (Prefix::size() < SubPreBasisPath::size()) {
      auto child = sub_pre_basis_path[index_constant<Prefix::size()>{}];
      static_assert(child < this->degree());
      return this->child(child).makeLocalView(ordering, push_back(prefix, child), sub_pre_basis_path);
    } else {
      auto unfold_children = [&](auto... i) {
        using LPreBasis = LocalNode<std::decay_t<decltype(*this->child(i).makeLocalView(ordering, push_back(prefix, i), sub_pre_basis_path))>...>;
        typename LPreBasis::NodeStorage storage{ this->child(i).makeLocalView(ordering, push_back(prefix, i), sub_pre_basis_path)... };
        return std::make_unique<LPreBasis>(std::move(storage));
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
    using ChildContainers = std::tuple<decltype(ChildOrdering::template makeVectorContainer<Traits>())...>;

    auto common_child_container = unpackIntegerSequence(
      [](auto... i) {
        return Std::is_detected<HasCommonType, typename Traits::template block_type<std::tuple_element_t<i, ChildContainers>>::type...>{};
      },std::index_sequence_for<ChildOrdering...>{});

    if constexpr (OrderingNode::containerBlocked()) {
      if constexpr (common_child_container) {
          using ChildBlock = std::common_type_t<decltype(ChildOrdering::template makeVectorContainer<Traits>())...>;
          // return Traits::template makeArray<ChildBlock, sizeof...(ChildOrdering)>(); // help! FieldVector<BlockVector<>> is a bad idea for the solvers :()
          return Traits::template makeVector<ChildBlock>();
      } else {
        return Traits::template makeTuple<decltype(ChildOrdering::template makeVectorContainer<Traits>())...>();
      }
    } else {
      return unpackIntegerSequence(
        [=](auto... i) {
          static_assert(common_child_container, "Non-Blocked tuple children should yield a common container type");
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

template<Concept::Impl::PreBasisTree PreBasis>
auto
makeLexicographicOrdering(const PreBasis& pre_basis)
{

  using NodeTag = typename PreBasis::NodeTag;
  using MergingStrategy =
    typename PreBasis::Traits::MergingStrategy;

  static_assert(not std::is_same<NodeTag, TypeTree::LeafNodeTag>{},
    "Lexicographic merging strategy shall not be applied to leaf nodes. "
    "Consider using merging by entity."
  );

  auto make_ordering = [](const auto& pre_basis){
    return makeOrdering(pre_basis, pre_basis.mergingStrategy());
  };

  if constexpr (std::is_same<NodeTag, TypeTree::PowerNodeTag>{}) {
    constexpr std::size_t degree = PreBasis::degree();
    using Child = std::decay_t<decltype(*make_ordering(pre_basis.child(0)))>;
    using LexicographicOrdering = ArrayLexicographicOrdering<MergingStrategy, Child, degree>;
    std::array<std::shared_ptr<Child>, degree> storage;
    for (std::size_t i = 0; i < degree; ++i)
      storage[i] = make_ordering(pre_basis.child(i));
    return std::make_unique<LexicographicOrdering>( std::move(storage) );
  } else if constexpr (std::is_same<NodeTag, TypeTree::DynamicPowerNodeTag>{}) {
    std::size_t degree = pre_basis.degree();
    using Child = std::decay_t<decltype(*make_ordering(pre_basis.child(0)))>;
    using LexicographicOrdering =
      VectorLexicographicOrdering<MergingStrategy, Child>;
    std::vector<std::shared_ptr<Child>> storage(degree);
    for (std::size_t i = 0; i < degree; ++i)
      storage[i] = make_ordering(pre_basis.child(i));
    return std::make_unique<LexicographicOrdering>( std::move(storage) );
  } else {
    static_assert(std::is_same<NodeTag, TypeTree::CompositeNodeTag>{});
    auto unfold_children = [&](auto... i) {
      using LexicographicOrdering = TupleLexicographicOrdering<MergingStrategy, std::decay_t<decltype(*make_ordering(pre_basis.child(i)))>...>;
      auto storage = std::tuple{ make_ordering(pre_basis.child(i))... };
      return std::make_unique<LexicographicOrdering>( std::move(storage) );
    };
    auto indices = std::make_index_sequence<PreBasis::degree()>{};
    return unpackIntegerSequence(unfold_children, indices);
  }
}

template<Concept::Impl::PreBasisTree PreBasis>
auto
makeOrdering(const PreBasis& pre_basis, Strategy::FlatLexicographic)
{
  return makeLexicographicOrdering(pre_basis);
}

template<Concept::Impl::PreBasisTree PreBasis>
auto
makeOrdering(const PreBasis& pre_basis, Strategy::BlockedLexicographic)
{
  return makeLexicographicOrdering(pre_basis);
}

} // namespace Dune::PDELab::inline Experimental::Impl

#endif // DUNE_PDELAB_BASIS_ORDERING_LEXICOGRAPHIC_HH
