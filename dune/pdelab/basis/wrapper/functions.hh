#ifndef DUNE_PDELAB_BASIS_WRAPPER_FUNCTIONS_HH
#define DUNE_PDELAB_BASIS_WRAPPER_FUNCTIONS_HH

#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/concepts/indexable.hh>

#include <dune/pdelab/common/multiindex.hh>

#include <dune/pdelab/basis/constraints/container.hh>

#include <dune/typetree/treecontainer.hh>
#include <dune/typetree/leafnode.hh>
#include <dune/typetree/powernode.hh>
#include <dune/typetree/compositenode.hh>
#include <dune/typetree/treepath.hh>

#if !HAVE_DUNE_FUNCTIONS
#error "This header is only available if dune-functions headers are found"
#endif

#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/backends/istlvectorbackend.hh>

#ifndef DUNE_PDELAB_ENABLE_DOUBLE_BIND
#define DUNE_PDELAB_ENABLE_DOUBLE_BIND 1
#endif

namespace Dune::PDELab::inline Experimental::Functions {

template<class LocalViewNode, class LocalView, class TreePath>
class LocalBasisLeaf : public TypeTree::LeafNode
{

  using MI = typename LocalView::MultiIndex;

public:
  using size_type = typename LocalView::size_type;
  using Element = typename LocalViewNode::Element;
  using Path = TreePath;
  using FiniteElement = typename LocalViewNode::FiniteElement;
  using MultiIndex = Dune::PDELab::MultiIndex<typename MI::value_type, MI::max_size()>;
  using ConstraintWeight = double;

  struct Traits
  {
    using FiniteElement = typename LocalViewNode::FiniteElement;
    using FiniteElementType = typename LocalViewNode::FiniteElement;
  };

  LocalBasisLeaf(const LocalViewNode& node, const LocalView& local_view, Path path)
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

  [[nodiscard]] friend decltype(auto) localContainerEntry(auto& container, const LocalBasisLeaf& lbasis, size_type dof) noexcept {
    auto mi = lbasis._local_view.index(lbasis._node.localIndex(dof));
    return Dune::Functions::istlVectorBackend(container)[mi];
  }

private:
  const LocalViewNode& _node;
  const LocalView& _local_view;
  [[no_unique_address]] Path _path;
};

template<class LocaBasisNode, std::size_t k>
class LocalBasisArray : public TypeTree::PowerNode<LocaBasisNode, k>
{
  using TreeNode = TypeTree::PowerNode<LocaBasisNode, k>;

public:
  LocalBasisArray(const typename TreeNode::NodeStorage& storage)
    : TreeNode{ storage }
  {
  }
};

template<class... LocaBasisNode>
class LocalBasisTuple : public TypeTree::CompositeNode<LocaBasisNode...>
{
  using TreeNode = TypeTree::CompositeNode<LocaBasisNode...>;

public:
  LocalBasisTuple(const typename TreeNode::NodeStorage& storage)
    : TreeNode{ storage }
  {
  }
};

template<class ViewNode, class View, class Path>
auto
makeLocalBasisTree(const ViewNode& node, const View& local_view, Path path)
{
  if constexpr (ViewNode::isLeaf) {
    using LocalView = LocalBasisLeaf<ViewNode, View, Path>;
    return std::make_unique<LocalView>(node, local_view, path);
  } else if constexpr (ViewNode::isPower) {
    using ChidlNode =
      std::decay_t<decltype(*makeLocalBasisTree(node.child(0), local_view, push_back(path,0)))>;
    using LocalView = LocalBasisArray<ChidlNode, ViewNode::degree()>;
    typename LocalView::NodeStorage storage;
    for (std::size_t i = 0; i < node.degree(); ++i)
      storage[i] = makeLocalBasisTree(node.child(i), local_view, push_back(path,i));
    return std::make_unique<LocalView>(storage);
  } else {
    static_assert(ViewNode::isComposite);
    return unpackIntegerSequence(
      [&](auto... i) {
        using ChildNodes = std::tuple<std::decay_t<decltype(*makeLocalBasisTree(
          node.child(i), local_view, push_back(path,i)))>...>;
        using LocalView =
          LocalBasisTuple<std::tuple_element_t<i, ChildNodes>...>;
        typename LocalView::NodeStorage storage{ makeLocalBasisTree(
          node.child(i), local_view, push_back(path,i))... };
        return std::make_unique<LocalView>(storage);
      },
      std::make_index_sequence<ViewNode::degree()>{});
  }
}

template<class View>
auto
makeLocalBasisTree(const View& local_view)
{
  return makeLocalBasisTree(local_view.tree(), local_view, TypeTree::treePath());
}

template<class WrappedBasis, Concept::Tree ConstraintsTree>
class Basis
{
  using MI = typename WrappedBasis::MultiIndex;
  using SP = typename WrappedBasis::SizePrefix;

public:
  using MultiIndex = Dune::PDELab::MultiIndex<typename MI::value_type, MI::max_size()>;
  using SizePrefix = Dune::PDELab::MultiIndex<typename SP::value_type, MI::max_size()>;
  using EntitySet = typename WrappedBasis::GridView;
  using size_type = typename WrappedBasis::size_type;

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

  using LocalBasisTree = std::decay_t<decltype(*makeLocalBasisTree(std::declval<typename WrappedBasis::LocalView>()))>;
  using RootConstraintsContainerStorage = decltype(makeConstraintsContainer(std::declval<const LocalBasisTree&>(), std::declval<ConstraintsContainerGenerator>()));
  using RootConstraintsContainer = typename RootConstraintsContainerStorage::element_type;
public:

  using LocalConstraints = decltype(std::declval<RootConstraintsContainer>().localView(std::declval<LocalBasisTree>(),TypeTree::treePath()));


  Basis(std::shared_ptr<WrappedBasis> basis_ptr, const ConstraintsTree& constraints_tree)
    : _wrapped_basis{ std::move(basis_ptr) }
    , _constraints_container{makeConstraintsContainer(*makeLocalBasisTree(_wrapped_basis->localView()), ConstraintsContainerGenerator{_wrapped_basis->gridView()})}
    , _constraints_tree{std::make_shared<ConstraintsTree>(constraints_tree)}
  {
    updateConstraints();
  }

  Basis(const WrappedBasis& basis, const ConstraintsTree& constraints_tree)
    : Basis{std::make_shared<WrappedBasis>(basis), constraints_tree}
  {}

  class LocalView
  {
    using BasisLocalView = typename WrappedBasis::LocalView;
  public:
    using GlobalBasis = Basis;
    using MultiIndex = Dune::PDELab::MultiIndex<typename MI::value_type, MI::max_size()>;
    using Element = typename BasisLocalView::Element;
    using Tree = LocalBasisTree;
    using size_type = typename WrappedBasis::size_type;

    LocalView(const Basis& basis)
      : _basis{basis}
      , _local_view{ std::make_unique<BasisLocalView>(_basis._wrapped_basis->localView()) }
      , _ltree_storage{ makeLocalBasisTree(*_local_view) }
    {}

    LocalView(const LocalView& local_basis) : LocalView{local_basis._basis} {}

    LocalView(LocalView&& other)
      : _basis{ std::move(other._basis) }
      , _local_view{ std::move(other._local_view) }
      , _ltree_storage{ std::move(other._ltree_storage) }
      , _ltree_view{ other._ltree_view }
    {}

    template<std::convertible_to<Element> E>
    LocalView& bind(E&& element) {
      bindElement(std::forward<E>(element));
      _local_view->bind(this->element());
      _ltree_view = _ltree_storage.get();
      // _mem_region = _basis.entitySet().memoryRegion(element);
      return *this;
    }

    LocalView& unbind() {
      _local_view->unbind();
      _entity_view = nullptr;
      _entity_storage = std::nullopt;
      return *this;
    }

    friend void unbind(LocalView& lbasis0, auto& lbasis1) {
      lbasis1.unbind();
      lbasis0.unbind();
    }


    template<class Element>
    friend void bind(Element&& element, LocalView& lbasis0, auto& lbasis1) {
      lbasis0.bind(std::forward<Element>(element));
      lbasis1.bind(lbasis0.element());
    }

#if DUNE_PDELAB_ENABLE_DOUBLE_BIND
    template<class Element>
    friend void bind(Element&& element, LocalView& lbasis0, LocalView& lbasis1) {
      lbasis0.doubleBind(std::forward<Element>(element), lbasis1);
    }

    friend void unbind(LocalView& lbasis0, LocalView& lbasis1) {
      lbasis0.doubleUnbind(lbasis1);
      lbasis0._entity_view = lbasis1._entity_view = nullptr;
      lbasis0._entity_storage = std::nullopt;
      lbasis1._entity_storage = std::nullopt;
    }
#endif


    [[nodiscard]] size_type size() const noexcept { return _local_view->size(); };

    [[nodiscard]] const Tree& tree() const noexcept { return *_ltree_view; }
    [[nodiscard]] const Element& element() const noexcept { return *_entity_view; }
    [[nodiscard]] size_type maxSize() const noexcept { return _local_view->maxSize(); }

    [[nodiscard]] const GlobalBasis& globalBasis() const noexcept { return _basis; }


    // [[nodiscard]] std::convertible_to<MemoryRegion> auto memoryRegion() const noexcept {
    //   return _mem_region; // TODO!
    // }

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
      if (_basis == other._basis) {
        other._ltree_view = _ltree_view;
        // other._mem_region = _mem_region;
      } else
        other.bind(this->element());
    }


    void doubleUnbind(LocalView& other) {
      unbind();
      if (_basis == other._basis)
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

    Basis _basis;
    std::unique_ptr<BasisLocalView> _local_view;
    std::unique_ptr<Tree> _ltree_storage;
    Tree const * _ltree_view = nullptr;
    Element const* _entity_view;
    std::optional<Element> _entity_storage;
    // MemoryRegion _mem_region;
  };

  // local index set is not implemented at all for dune-functions because there
  // is no know way (to me) to extract the local indices of an entity without
  // reconstructing the whole map
  // we can reuse the local view to at least honor the basis interface
  struct LocalIndexSet : public LocalView {
    LocalIndexSet() = delete;

    void bind(const Dune::Concept::Entity auto& entity) { DUNE_THROW(NotImplemented, ""); }

    friend void bind(const Dune::Concept::Entity auto& element, LocalView& lbasis0, auto& lbasis1) {
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
    return _constraints_container->localView(localView().tree(),TypeTree::treePath());
  }

  [[nodiscard]] const EntitySet& entitySet() const noexcept { return _wrapped_basis->gridView(); }

  [[nodiscard]] const size_type dimension() const noexcept { return _wrapped_basis->dimension(); }

  [[nodiscard]] size_type size(const SizePrefix& prefix) const noexcept
  {
    return _wrapped_basis->size(prefix);
  }

  [[nodiscard]] auto degree() const noexcept {
    return WrappedBasis::LocalView::Tree::degree();
  }

  [[nodiscard]] size_type size() const noexcept
  {
    return _wrapped_basis->size();
  }

  [[nodiscard]] bool fixedSize(std::size_t dim, std::size_t codim) const {
    DUNE_THROW(NotImplemented, "");
   }

  [[nodiscard]] bool contains(std::size_t dim, std::size_t codim) const {
    DUNE_THROW(NotImplemented, "");
   }

  void update(const EntitySet& entity_set) {
    _wrapped_basis->update(entity_set);
    updateConstraints();
  }

  [[nodiscard]] friend bool operator==(const Basis&, const Basis&) = default;
  [[nodiscard]] friend bool operator!=(const Basis&, const Basis&) = default;

private:
  void updateConstraints() {
    _constraints_container->assembleConstraints(*this, TypeTree::makeTreeContainer(*_constraints_tree, [](auto node){return node;}));
  }
  std::shared_ptr<WrappedBasis> _wrapped_basis;
  std::shared_ptr<RootConstraintsContainer> _constraints_container;
  std::shared_ptr<ConstraintsTree> _constraints_tree;
};

} // Dune::PDELab::inline Experimental::Functions

#endif // DUNE_PDELAB_BASIS_WRAPPER_FUNCTIONS_HH
