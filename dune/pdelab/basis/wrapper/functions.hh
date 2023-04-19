#ifndef DUNE_ASSEMBLER_DISCRETE_FUNCTION_SPACE_WRAPPER_FUNCTIONS_HH
#define DUNE_ASSEMBLER_DISCRETE_FUNCTION_SPACE_WRAPPER_FUNCTIONS_HH

#include <dune/assembler/concepts/space.hh>
#include <dune/assembler/concepts/indexable.hh>

#include <dune/assembler/common/reservedmultiindex.hh>
#include <dune/assembler/common/multiindex.hh>
#include <dune/assembler/common/entityset.hh>

#include <dune/assembler/space/constraints/container.hh>

#include <dune/typetree/treecontainer.hh>
#include <dune/typetree/leafnode.hh>
#include <dune/typetree/powernode.hh>
#include <dune/typetree/compositenode.hh>

#if !HAVE_DUNE_FUNCTIONS
#error "This header is only available if dune-functions headers are found"
#endif

#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/backends/istlvectorbackend.hh>

#ifndef DUNE_ASSEMBLER_ENABLE_DOUBLE_BIND
#define DUNE_ASSEMBLER_ENABLE_DOUBLE_BIND 1
#endif

namespace Dune::Assembler::Functions {

template<class LocalViewNode, class LocalView, class TreePath>
class LeafLocalSpace : public TypeTree::LeafNode
{

  using MI = typename LocalView::MultiIndex;

public:
  using size_type = typename LocalView::size_type;
  using Element = typename LocalViewNode::Element;
  using Path = TreePath;
  using FiniteElement = typename LocalViewNode::FiniteElement;
  using MultiIndex =
    ReservedMultiIndex<typename MI::value_type, MI::max_size()>;
  using ConstraintWeight = double;

  struct Traits
  {
    using FiniteElement = typename LocalViewNode::FiniteElement;
    using FiniteElementType = typename LocalViewNode::FiniteElement;
  };

  LeafLocalSpace(const LocalViewNode& node, const LocalView& local_view, Path path)
    : _node{ node }
    , _local_view{ local_view }
    , _path{ path }
  {
  }

  [[nodiscard]] MultiIndex index(size_type dof) const noexcept
  {
    MI source = _local_view.index(this->localIndex(dof));
    MultiIndex target;
    if constexpr (MI::max_size() == 1)
      target[0] = source[0];
    else {
      target.resize(source.size());
      std::copy(source.begin(), source.end(), target.begin());
    }
    return target;
  }

  [[nodiscard]] size_type localIndex(size_type dof) const noexcept { return _node.localIndex(dof); }

  [[nodiscard]] size_type size() const noexcept { return _node.size(); }

  [[nodiscard]] const Element& element() const noexcept { return _node.element(); }

  [[nodiscard]] const FiniteElement& finiteElement() const noexcept { return _node.finiteElement(); }


  [[nodiscard]] Path path() const noexcept { return _path; }

  [[nodiscard]] friend decltype(auto) localContainerEntry(auto& container, const LeafLocalSpace& lspace, size_type dof) noexcept {
    auto mi = lspace._local_view.index(lspace._node.localIndex(dof));
    return Dune::Functions::istlVectorBackend(container)[mi];
  }

private:
  const LocalViewNode& _node;
  const LocalView& _local_view;
  [[no_unique_address]] Path _path;
};

template<class LocaSpaceNode, std::size_t k>
class ArrayLocalSpace : public TypeTree::PowerNode<LocaSpaceNode, k>
{
  using TreeNode = TypeTree::PowerNode<LocaSpaceNode, k>;

public:
  ArrayLocalSpace(const typename TreeNode::NodeStorage& storage)
    : TreeNode{ storage }
  {
  }
};

template<class... LocaSpaceNode>
class TupleLocalSpace : public TypeTree::CompositeNode<LocaSpaceNode...>
{
  using TreeNode = TypeTree::CompositeNode<LocaSpaceNode...>;

public:
  TupleLocalSpace(const typename TreeNode::NodeStorage& storage)
    : TreeNode{ storage }
  {
  }
};

template<class ViewNode, class View, class Path>
auto
makeLocalSpaceTree(const ViewNode& node, const View& local_view, Path path)
{
  if constexpr (ViewNode::isLeaf) {
    using LocalView = LeafLocalSpace<ViewNode, View, Path>;
    return std::make_unique<LocalView>(node, local_view, path);
  } else if constexpr (ViewNode::isPower) {
    using ChidlNode =
      std::decay_t<decltype(*makeLocalSpaceTree(node.child(0), local_view, push_back(path,0)))>;
    using LocalView = ArrayLocalSpace<ChidlNode, ViewNode::degree()>;
    typename LocalView::NodeStorage storage;
    for (std::size_t i = 0; i < node.degree(); ++i)
      storage[i] = makeLocalSpaceTree(node.child(i), local_view, push_back(path,i));
    return std::make_unique<LocalView>(storage);
  } else {
    static_assert(ViewNode::isComposite);
    return unpackIntegerSequence(
      [&](auto... i) {
        using ChildNodes = std::tuple<std::decay_t<decltype(*makeLocalSpaceTree(
          node.child(i), local_view, push_back(path,i)))>...>;
        using LocalView =
          TupleLocalSpace<std::tuple_element_t<i, ChildNodes>...>;
        typename LocalView::NodeStorage storage{ makeLocalSpaceTree(
          node.child(i), local_view, push_back(path,i))... };
        return std::make_unique<LocalView>(storage);
      },
      std::make_index_sequence<ViewNode::degree()>{});
  }
}

template<class View>
auto
makeLocalSpaceTree(const View& local_view)
{
  return makeLocalSpaceTree(local_view.tree(), local_view, multiIndex());
}

template<class Basis, Concept::Tree ConstraintsTree>
class Space
{
  using MI = typename Basis::MultiIndex;
  using SP = typename Basis::SizePrefix;

public:
  using MultiIndex = ReservedMultiIndex<typename MI::value_type, MI::max_size()>;
  using SizePrefix = ReservedMultiIndex<typename SP::value_type, MI::max_size()>;
  using EntitySet = typename Basis::GridView;
  using size_type = typename Basis::size_type;

private:

  struct ConstraintsContainerGenerator {
    template<Concept::LeafTreeNode LeafNode, Concept::MultiIndex Path>
    auto operator()(const LeafNode&, Path) {
      using Constraints = TypeTree::ChildForTreePath<ConstraintsTree, Path>;
      using ConstrainsContainer = typename Constraints::template Container<MultiIndex,EntitySet>;
      return std::make_shared<ConstrainsContainer>(_entity_set);
    }

    EntitySet _entity_set;
  };

  using LocalSpaceTree = std::decay_t<decltype(*makeLocalSpaceTree(std::declval<typename Basis::LocalView>()))>;
  using RootConstraintsContainerStorage = decltype(makeConstraintsContainer(std::declval<const LocalSpaceTree&>(), std::declval<ConstraintsContainerGenerator>()));
  using RootConstraintsContainer = typename RootConstraintsContainerStorage::element_type;
public:

  using LocalConstraints = decltype(std::declval<RootConstraintsContainer>().localView(std::declval<LocalSpaceTree>(),multiIndex()));


  Space(std::shared_ptr<Basis> basis_ptr, const ConstraintsTree& constraints_tree)
    : _basis{ std::move(basis_ptr) }
    , _constraints_container{makeConstraintsContainer(*makeLocalSpaceTree(_basis->localView()), ConstraintsContainerGenerator{_basis->gridView()})}
    , _constraints_tree{std::make_shared<ConstraintsTree>(constraints_tree)}
  {
    updateConstraints();
  }

  Space(const Basis& basis, const ConstraintsTree& constraints_tree)
    : Space{std::make_shared<Basis>(basis), constraints_tree}
  {}

  class LocalView
  {
    using BasisLocalView = typename Basis::LocalView;
  public:
    using GlobalBasis = Space;
    using MultiIndex = ReservedMultiIndex<typename MI::value_type, MI::max_size()>;
    using Element = typename BasisLocalView::Element;
    using Tree = LocalSpaceTree;
    using size_type = typename Basis::size_type;

    LocalView(const Space& space)
      : _space{space}
      , _local_view{ std::make_unique<BasisLocalView>(_space._basis->localView()) }
      , _ltree_storage{ makeLocalSpaceTree(*_local_view) }
    {}

    LocalView(const LocalView& local_space) : LocalView{local_space._space} {}

    LocalView(LocalView&& other)
      : _space{ std::move(other._space) }
      , _local_view{ std::move(other._local_view) }
      , _ltree_storage{ std::move(other._ltree_storage) }
      , _ltree_view{ other._ltree_view }
    {}

    template<std::convertible_to<Element> E>
    LocalView& bind(E&& element) {
      bindElement(std::forward<E>(element));
      _local_view->bind(this->element());
      _ltree_view = _ltree_storage.get();
      _mem_region = _space.entitySet().memoryRegion(element);
      return *this;
    }

    LocalView& unbind() {
      _local_view->unbind();
      _entity_view = nullptr;
      _entity_storage = std::nullopt;
      return *this;
    }

    friend void unbind(LocalView& lspace0, auto& lspace1) {
      lspace1.unbind();
      lspace0.unbind();
    }


    template<class Element>
    friend void bind(Element&& element, LocalView& lspace0, auto& lspace1) {
      lspace0.bind(std::forward<Element>(element));
      lspace1.bind(lspace0.element());
    }

#if DUNE_ASSEMBLER_ENABLE_DOUBLE_BIND
    template<class Element>
    friend void bind(Element&& element, LocalView& lspace0, LocalView& lspace1) {
      lspace0.doubleBind(std::forward<Element>(element), lspace1);
    }

    friend void unbind(LocalView& lspace0, LocalView& lspace1) {
      lspace0.doubleUnbind(lspace1);
      lspace0._entity_view = lspace1._entity_view = nullptr;
      lspace0._entity_storage = std::nullopt;
      lspace1._entity_storage = std::nullopt;
    }
#endif


    [[nodiscard]] size_type size() const noexcept { return _local_view->size(); };

    [[nodiscard]] const Tree& tree() const noexcept { return *_ltree_view; }
    [[nodiscard]] const Element& element() const noexcept { return *_entity_view; }
    [[nodiscard]] size_type maxSize() const noexcept { return _local_view->maxSize(); }

    [[nodiscard]] const GlobalBasis& globalBasis() const noexcept { return _space; }


    [[nodiscard]] std::convertible_to<MemoryRegion> auto memoryRegion() const noexcept {
      return _mem_region; // TODO!
    }

    // Whether local index sets match in all processors
    [[nodiscard]] auto conforming() const noexcept {
      return std::false_type{}; // not conforming in general...
    }

    [[nodiscard]] MultiIndex index(size_type dof) const noexcept {
      return _local_view->index(dof);
    }

  private:

    template<class Element>
    void doubleBind(Element&& element, LocalView& other) {
      bind(std::forward<Element>(element));
      if (_space == other._space) {
        other._ltree_view = _ltree_view;
        other._mem_region = _mem_region;
      } else
        other.bind(this->element());
    }


    void doubleUnbind(LocalView& other) {
      unbind();
      if (_space == other._space)
        other._ltree_view = nullptr;
      else
        other.unbind();
    }

    void bindElement(Element&& element) {
      // the caller assigned the ownership to us
      _entity_storage.emplace(std::move(element));
      _entity_view = &(*_entity_storage);
    }

    void bindElement(const Element& element) {
        // ownership is managed elsewhere
      _entity_view = &element;
    }

    Space _space;
    std::unique_ptr<BasisLocalView> _local_view;
    std::unique_ptr<Tree> _ltree_storage;
    Tree const * _ltree_view = nullptr;
    Element const* _entity_view;
    std::optional<Element> _entity_storage;
    MemoryRegion _mem_region;
  };

  // local index set is not implemented at all for dune-functions because there
  // is no know way (to me) to extract the local indices of an entity without
  // reconstructing the whole map
  // we can reuse the local view to at least honor the space interface
  struct LocalIndexSet : public LocalView {
    LocalIndexSet() = delete;

    void bind(const Dune::Concept::Entity auto& entity) { DUNE_THROW(NotImplemented, ""); }

    friend void bind(const Dune::Concept::Entity auto& element, LocalView& lspace0, auto& lspace1) {
      DUNE_THROW(NotImplemented, "");
    }

    LocalIndexSet& bind(const Dune::Concept::Entity auto& element) {
      DUNE_THROW(NotImplemented, "");
    }

    LocalIndexSet& unbind() {
      DUNE_THROW(NotImplemented, "");
    }
  };

  [[nodiscard]] LocalView localView() const { return LocalView{ *this }; }

  [[nodiscard]] LocalIndexSet localIndexSet() const {
    DUNE_THROW(NotImplemented, "");
   }

  [[nodiscard]] LocalConstraints localConstraints() const {
    return _constraints_container->localView(localView().tree(),multiIndex());
  }

  [[nodiscard]] const EntitySet& entitySet() const noexcept { return _basis->gridView(); }

  [[nodiscard]] const size_type dimension() const noexcept { return _basis->dimension(); }

  [[nodiscard]] size_type size(const SizePrefix& prefix) const noexcept
  {
    return _basis->size(prefix);
  }

  [[nodiscard]] auto degree() const noexcept {
    return Basis::LocalView::Tree::degree();
  }

  [[nodiscard]] size_type size() const noexcept
  {
    return _basis->size();
  }

  [[nodiscard]] bool fixedSize(std::size_t dim, std::size_t codim) const {
    DUNE_THROW(NotImplemented, "");
   }

  [[nodiscard]] bool contains(std::size_t dim, std::size_t codim) const {
    DUNE_THROW(NotImplemented, "");
   }

  void update(const EntitySet& entity_set) {
    _basis->update(entity_set);
    updateConstraints();
  }

  [[nodiscard]] friend bool operator==(const Space&, const Space&) = default;
  [[nodiscard]] friend bool operator!=(const Space&, const Space&) = default;

private:
  void updateConstraints() {
    _constraints_container->assembleConstraints(*this, TypeTree::makeTreeContainer(*_constraints_tree, [](auto node){return node;}));
  }
  std::shared_ptr<Basis> _basis;
  std::shared_ptr<RootConstraintsContainer> _constraints_container;
  std::shared_ptr<ConstraintsTree> _constraints_tree;
};

} // Dune::Assembler::Functions

#endif // DUNE_ASSEMBLER_DISCRETE_FUNCTION_SPACE_WRAPPER_FUNCTIONS_HH
