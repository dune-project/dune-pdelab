#ifndef DUNE_PDELAB_BASIS_WRAPPER_GRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_BASIS_WRAPPER_GRIDFUNCTIONSPACE_HH

#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/concepts/indexable.hh>

#include <dune/pdelab/common/multiindex.hh>
// #include <dune/pdelab/common/communication/entity_data_handler.hh>

#include <dune/pdelab/basis/constraints/container.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>

#include <dune/typetree/treecontainer.hh>
#include <dune/typetree/treepath.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <array>

#ifndef DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE_LEGACY
#define DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE_LEGACY 0
#endif

#ifndef DUNE_PDELAB_ENABLE_DOUBLE_BIND
#define DUNE_PDELAB_ENABLE_DOUBLE_BIND 1
#endif


namespace Dune::PDELab::inline Experimental::Legacy {

template<class Cache, class TreePath, bool fast>
class LeafLocaIndexSet : public TypeTree::LeafNode
{
public:
  LeafLocaIndexSet(TreePath path)
    : _path{ path }
  {}

  using Path = TreePath;
  using size_type = typename Cache::size_type;
  using MultiIndex = typename Cache::ContainerIndex;

  [[nodiscard]] size_type localIndex(size_type i) const { return _cache_view->offsets()[_tree_index] + i; }

  [[nodiscard]] size_type size() const { return _cache_view->offsets()[_tree_index+1] - _cache_view->offsets()[_tree_index]; }

  [[nodiscard]] Path path() const {
    return _path;
  }

  void setTreeOffset(std::size_t tree_offset, std::size_t tree_index) {
    _tree_index = tree_index;
    assert(tree_offset == _cache_view->offsets()[_tree_index]);
  }

  [[nodiscard]] MultiIndex index(size_type dof) const
  {
    MultiIndex ci = _cache_view->containerIndex(this->localIndex(fast ? 0 : dof));
    if constexpr (fast) ci.front() += dof;
    // dune-assembler assumes a reversed order of indices wrt dune-pdelab.
    // Reversing them is quite inefficient. That's why we provide a friend
    // function to index the container for then native handler.
    if constexpr (MultiIndex::max_size() > 1)
      return ci;
    else
      return reverse(ci);
  }

  [[nodiscard]] friend decltype(auto) localContainerEntry(
                           auto& container,
                           const LeafLocaIndexSet& lbasis,
                           size_type dof)
  {
    MultiIndex ci = lbasis._cache_view->containerIndex(lbasis.localIndex(fast ? 0 : dof));
    if constexpr (fast) ci.front() += dof;
    if constexpr (requires {container[ci];}) // PDELab backend
      return container[ci];
    else // standard indexable container
      return localContainerEntryImpl(container, ci);
  }

  void setCacheView(Cache const * cache_view) {
    _cache_view = cache_view;
  }

private:

  static decltype(auto) localContainerEntryImpl(auto&& container, MultiIndex prefix)
  {
    if constexpr (requires {container[back(prefix)];}) {
      assert(not prefix.empty());
      auto index = back(prefix);
      return localContainerEntryImpl(container[index], pop_back(prefix));
    } else {
      return container;
    }
  }

  std::size_t _tree_index;
  Cache const * _cache_view = nullptr;
  [[no_unique_address]] Path _path;
};


template<class Cache, class TreePath, bool fast>
class LocalBasisLeaf : public TypeTree::LeafNode
{
  using LFS = typename Cache::LocalFunctionSpace;
  using LeafLFS = TypeTree::ChildForTreePath<LFS, TreePath>;

public:
  LocalBasisLeaf(TreePath path)
    : _path{ path }
  {}

  using Path = TreePath;
  using Element = typename LeafLFS::Traits::Element;
  using FiniteElement = typename LeafLFS::Traits::FiniteElement;
  using size_type = typename LeafLFS::Traits::SizeType;
  using MultiIndex = Cache::ContainerIndex;

  struct Traits
  {
    using FiniteElement = typename LeafLFS::Traits::FiniteElement;
    using FiniteElementType = typename LeafLFS::Traits::FiniteElement;
  };

  void setCacheView(Cache const * cache_view) {
    _cache_view = cache_view;
    _leaf_lfs_view = &TypeTree::child(_cache_view->localFunctionSpace(), _path);
  }

  void bind(Element const* element) { _entity_view = element; }

  [[nodiscard]] const FiniteElement& finiteElement() const noexcept
  {
    return _leaf_lfs_view->finiteElement();
  }

  [[nodiscard]] const Element& element() const noexcept
  {
    assert(_entity_view);
    return *_entity_view;
  }

  [[nodiscard]] size_type localIndex(size_type i) const noexcept { return _tree_offset + i; }

  [[nodiscard]] size_type size() const noexcept { return _leaf_lfs_view->size(); }

  [[nodiscard]] Path path() const noexcept {
    return _path;
  }

  void setTreeOffset(std::size_t tree_offset, std::size_t tree_index) { _tree_offset = tree_offset; }

  [[nodiscard]] MultiIndex index(size_type dof) const noexcept
  {
    MultiIndex ci;
    if constexpr (fast) {
      ci = _cache_view->containerIndex(this->localIndex(0));
      ci.front() += dof;
    } else {
      ci =_cache_view->containerIndex(this->localIndex(dof));
    }
    // dune-assembler assumes a reversed order of indices wrt dune-pdelab.
    // Reversing them is quite inefficient. That's why we provide a friend
    // function to index the container for then native handler.
    if constexpr (MultiIndex::max_size() == 1)
      return ci;
    else
      return reverse(ci);
  }

  [[nodiscard]] friend decltype(auto) localContainerEntry(
                           auto& container,
                           const LocalBasisLeaf& lbasis,
                           size_type dof) noexcept
  {
    MultiIndex ci;
    if constexpr (fast) {
      ci = lbasis._cache_view->containerIndex(lbasis.localIndex(0));
      ci.front() += dof;
    } else {
      ci = lbasis._cache_view->containerIndex(lbasis.localIndex(dof));
    }
    if constexpr (requires {container[ci];}) // PDELab backend
      return container[ci];
    else // standard indexable container
      return localContainerEntryImpl(container, ci, ci.size()-1 );
  }

private:

  static decltype(auto) localContainerEntryImpl(auto&& container, const MultiIndex& prefix, std::size_t i) noexcept
  {
    if constexpr (requires {container[prefix[i]];}) {
      return localContainerEntryImpl(container[prefix[i]], prefix, i-1);
    } else {
      return container;
    }
  }

  std::size_t _tree_offset;
  Cache const * _cache_view = nullptr;
  LeafLFS const * _leaf_lfs_view = nullptr;
  Element const* _entity_view = nullptr;
  [[no_unique_address]] Path _path;
};

template<class T, std::size_t k>
class LocalBasisArray : public TypeTree::PowerNode<T, k>
{
  using TreeNode = TypeTree::PowerNode<T, k>;

public:
  LocalBasisArray(const typename TreeNode::NodeStorage& storage)
    : TreeNode{ storage }
  {
  }
};

template<class... T>
class LocalBasisTuple : public TypeTree::CompositeNode<T...>
{
  using TreeNode = TypeTree::CompositeNode<T...>;

public:
  LocalBasisTuple(const typename TreeNode::NodeStorage& storage)
    : TreeNode{ storage }
  {
  }
};

template<bool fast, class Cache, class GFSNode, class TreePath>
auto makeLocalBasisTree(const GFSNode& gfs_node, TreePath path)
{
  if constexpr (GFSNode::isLeaf) {
    using LocalView = LocalBasisLeaf<Cache, TreePath, fast>;
    return std::make_unique<LocalView>(path);
  } else if constexpr (GFSNode::isPower) {
    using ChidlNode = std::decay_t<decltype(*makeLocalBasisTree<fast,Cache>(
      gfs_node.child(0), push_back(path, 0)))>;
    using LocalView = LocalBasisArray<ChidlNode, GFSNode::degree()>;
    typename LocalView::NodeStorage storage;
    for (std::size_t i = 0; i < gfs_node.degree(); ++i)
      storage[i] =
        makeLocalBasisTree<fast,Cache>(gfs_node.child(i), push_back(path, i));
    return std::make_unique<LocalView>(storage);
  } else {
    static_assert(GFSNode::isComposite);
    return unpackIntegerSequence(
      [&](auto... i) {
        using ChildNodes = std::tuple<std::decay_t<decltype(*makeLocalBasisTree<fast,Cache>(
          gfs_node.child(i), push_back(path, i)))>...>;
        using LocalView =
          LocalBasisTuple<std::tuple_element_t<i, ChildNodes>...>;
        typename LocalView::NodeStorage storage{ makeLocalBasisTree<fast, Cache>(
          gfs_node.child(i), push_back(path, i))... };
        return std::make_unique<LocalView>(storage);
      },
      std::make_index_sequence<GFSNode::degree()>{});
  }
}

template<bool fast, class Cache, class GFSNode, class TreePath>
auto makeLocalIndexSetTree(const GFSNode& gfs_node, TreePath path)
{
  if constexpr (GFSNode::isLeaf) {
    using LocalView = LeafLocaIndexSet<Cache, TreePath, fast>;
    return std::make_unique<LocalView>(path);
  } else if constexpr (GFSNode::isPower) {
    using ChidlNode = std::decay_t<decltype(*makeLocalIndexSetTree<fast, Cache>(
      gfs_node.child(0), push_back(path, 0)))>;
    using LocalView = LocalBasisArray<ChidlNode, GFSNode::degree()>;
    typename LocalView::NodeStorage storage;
    for (std::size_t i = 0; i < gfs_node.degree(); ++i)
      storage[i] =
        makeLocalIndexSetTree<fast, Cache>(gfs_node.child(i), push_back(path, i));
    return std::make_unique<LocalView>(storage);
  } else {
    static_assert(GFSNode::isComposite);
    return unpackIntegerSequence(
      [&](auto... i) {
        using ChildNodes = std::tuple<std::decay_t<decltype(*makeLocalIndexSetTree<fast, Cache>(
          gfs_node.child(i), push_back(path, i)))>...>;
        using LocalView =
          LocalBasisTuple<std::tuple_element_t<i, ChildNodes>...>;
        typename LocalView::NodeStorage storage{ makeLocalIndexSetTree<fast, Cache>(
          gfs_node.child(i), push_back(path, i))... };
        return std::make_unique<LocalView>(storage);
      },
      std::make_index_sequence<GFSNode::degree()>{});
  }
}

template<class GridFunctionSpace, Concept::Tree ConstraintsTree>
class Basis
{

  template<class Node>
  constexpr static auto implFastDG() {
    // In case of fast DG, we asume no interleaving or permutation on the front
    if constexpr (Concept::LeafTreeNode<Node>) {
      constexpr int dim = Node::Traits::EntitySet::dimension;
      if (not Node::Traits::FiniteElementMap::hasDOFs(0))
        return false;
      bool other_codims = false;
      for (int codim = 1; codim <= dim; ++codim)
        other_codims |= Node::Traits::FiniteElementMap::hasDOFs(codim);
      return not other_codims;
    } else if constexpr (Concept::ArrayTreeNode<Node> ||
                         Concept::VectorTreeNode<Node>) {
      return implFastDG<typename Node::ChildType>();
    } else if constexpr (Node::isComposite) {
      auto unfold_children = [&](auto... i) {
        constexpr bool fast_dg =
          (implFastDG<TypeTree::template Child<Node, i>>() || ...);
        return std::integral_constant<bool, fast_dg>{};
      };
      auto indices = std::make_index_sequence<Node::degree()>{};
      return unpackIntegerSequence(unfold_children, indices);
    } else {
      static_assert(Dune::AlwaysFalse<Node>{}, "Not known Node Type");
    }
  }

  constexpr static auto fast = std::integral_constant<bool,implFastDG<GridFunctionSpace>()>{};

  using LFS = Dune::PDELab::LocalFunctionSpace<GridFunctionSpace>;
  using LFSCache = Dune::PDELab::LFSIndexCache<LFS,Dune::PDELab::EmptyTransformation,fast>;
  using LISCache = Dune::PDELab::EntityIndexCache<GridFunctionSpace, false>;

  using LocalBasisTree = std::decay_t<decltype(*makeLocalBasisTree<fast,LFSCache>(std::declval<GridFunctionSpace>(), TypeTree::treePath()))>;
  using LocalIndexSetTree = std::decay_t<decltype(*makeLocalIndexSetTree<fast,LISCache>(std::declval<GridFunctionSpace>(), TypeTree::treePath()))>;

public:
  using MultiIndex = typename GridFunctionSpace::Ordering::Traits::ContainerIndex;
  using SizePrefix = MultiIndex;
  using EntitySet = typename GridFunctionSpace::Traits::EntitySet;
  using size_type = typename GridFunctionSpace::Traits::SizeType;

private:
  static constexpr auto constraints_container_generator = []<Concept::LeafTreeNode LeafBasis, Concept::MultiIndex Path>(const LeafBasis& leaf_basis, Path path){
    using Constraints = TypeTree::ChildForTreePath<ConstraintsTree, Path>;
    using ConstrainsContainer = typename Constraints::template Container<MultiIndex, typename LeafBasis::Traits::EntitySet>;
    return std::make_shared<ConstrainsContainer>(leaf_basis.entitySet());
  };

  using RootConstraintsContainerStorage = decltype(makeConstraintsContainer(std::declval<const GridFunctionSpace&>(), constraints_container_generator));
  using RootConstraintsContainer = typename RootConstraintsContainerStorage::element_type;
public:

  using LocalConstraints = decltype(std::declval<RootConstraintsContainer>().localView(std::declval<LocalBasisTree>(),TypeTree::treePath()));

  Basis(std::shared_ptr<GridFunctionSpace> gfs_ptr, const ConstraintsTree& constraints_tree)
    : _gfs{ std::move(gfs_ptr) }
    , _constraints_container{makeConstraintsContainer(*_gfs, constraints_container_generator)}
    , _constraints_tree{std::make_shared<ConstraintsTree>(constraints_tree)}
    , _conforming_local_index_set{std::make_shared<bool>(false)}
  {
    updateConstraints();
  }

private:

  template<class LocalTree>
  struct LocalIndexSetBase {

    using MultiIndex = typename GridFunctionSpace::Ordering::Traits::ContainerIndex;
    using Tree = LocalTree;
    using size_type = typename LISCache::size_type;
    using GlobalBasis = Basis;

    LocalIndexSetBase(const Basis& basis, std::unique_ptr<LocalTree> ltree_storage)
      : _basis{basis}
      , _ltree_storage{move(ltree_storage)}
    {
        if constexpr (DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE_LEGACY)
          _indices = std::make_unique<MultiIndex[]>(maxSize());
    }

    LocalIndexSetBase(const LocalIndexSetBase&) = delete;

    LocalIndexSetBase(LocalIndexSetBase&& other)
      : _basis{ std::move(other._basis) }
      , _ltree_storage{ move(other._ltree_storage) }
      , _ltree_view{ other._ltree_view }
      // , _mem_region{ other._mem_region }
    {
      if constexpr (fast and DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE_LEGACY)
        _indices = move(other._indices);
    }

    void bind(const Dune::Concept::Entity auto& entity) {
      std::size_t tree_index = 0;
      std::size_t tree_offset = 0;
      PDELab::forEachLeafNode(*_ltree_storage, [&](auto& leaf){
        leaf.setTreeOffset(tree_offset, tree_index);
        if constexpr (fast and DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE_LEGACY)
          for(size_type dof = 0; dof != leaf.size(); ++dof)
            _indices[tree_offset + dof] = leaf.index(dof);
        tree_index  += 1;
        tree_offset += leaf.size();
      });
      assert(tree_offset <= maxSize());
      if constexpr (fast and DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE_LEGACY)
        _indices_view = _indices.get();
      _ltree_view = _ltree_storage.get();
      // _mem_region = _basis.entitySet().memoryRegion(entity);
    }

    void unbind() {
      _ltree_view = nullptr;
      _indices_view = nullptr;
    }

    // [[nodiscard]] std::convertible_to<MemoryRegion> auto memoryRegion() const noexcept {
    //   return _mem_region; // TODO!
    // }

    [[nodiscard]] size_type maxSize() const noexcept {
      return _basis._gfs->maxLocalSize();
    };

    [[nodiscard]] const Tree& tree() const noexcept {
      return *_ltree_view;
    }

    [[nodiscard]] const GlobalBasis& globalBasis() const noexcept {
      return _basis;
    }

    // Whether local index sets match in all processors
    [[nodiscard]] auto conforming() const noexcept {
      return _basis.isLocalIndexSetConforming();
    }

  protected:

    void doubleBind(const Dune::Concept::Entity auto& element, LocalIndexSetBase& other) {
      bind(element);
      if (_basis == other._basis) {
        other._ltree_view = _ltree_view;
        other._indices_view = _indices_view;
        // other._mem_region = _mem_region;
      } else {
        other.bind(element);
      }
    }

    void doubleUnbind(LocalIndexSetBase& other) {
      unbind();
      if (_basis == other._basis) {
        other._ltree_view = nullptr;
        other._indices_view = nullptr;
      } else {
        other.unbind();
      }
    }

    GlobalBasis _basis;
    std::unique_ptr<LocalTree> _ltree_storage;
    LocalTree const * _ltree_view = nullptr;
    std::unique_ptr<MultiIndex[]> _indices;
    MultiIndex const * _indices_view = nullptr;
    // MemoryRegion _mem_region;
  };


public:

  class LocalIndexSet : public LocalIndexSetBase<LocalIndexSetTree> {
    using Base = LocalIndexSetBase<LocalIndexSetTree>;
  public:

    LocalIndexSet(const Basis& basis)
      : Base{basis, makeLocalIndexSetTree<fast,LISCache>(*basis._gfs, TypeTree::treePath())}
      , _cache{ std::make_unique<LISCache>(*basis._gfs) }
    {
      PDELab::forEachLeafNode(*this->_ltree_storage, [&](auto& leaf){
        leaf.setCacheView(_cache.get());
      });
    }

    LocalIndexSet(const LocalIndexSet& localview)
      : LocalIndexSet( localview._basis )
    {}

    LocalIndexSet(LocalIndexSet&& other)
      : Base(std::move(other))
      , _cache( move(other._cache) )
    {}

    LocalIndexSet& bind(const Dune::Concept::Entity auto& element) {
      _cache->update(element);
      Base::bind(element);
      return *this;
    }

    LocalIndexSet& unbind() {
      return *this;
    }

    friend void bind(const Dune::Concept::Entity auto& element, LocalIndexSet& lbasis0, auto& lbasis1) {
      lbasis1.bind(element);
      lbasis0.bind(element);
    }

#if DUNE_PDELAB_ENABLE_DOUBLE_BIND
    friend void bind(const Dune::Concept::Entity auto& element, LocalIndexSet& lbasis0, LocalIndexSet& lbasis1) {
      lbasis0._cache->update(element);
      lbasis0.doubleBind(element, lbasis1);
    }

    friend void unbind(LocalIndexSet& lbasis0, LocalIndexSet& lbasis1) {
      lbasis0.doubleUnbind(lbasis1);
    }
#endif

    [[nodiscard]] size_type size() const noexcept {
      return _cache->size();
    }

    [[nodiscard]] MultiIndex index(size_type dof) const noexcept {
      if constexpr (fast and not DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE_LEGACY)
        DUNE_THROW(NotImplemented, "To enable this feature complie with `-DDUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE_LEGACY`");
      else if constexpr (fast)
        return this->_indices_view[dof];
      else
        return _cache->containerIndex(dof);
    }

  private:
    std::unique_ptr<LISCache> _cache;
  };

  class LocalView : public LocalIndexSetBase<LocalBasisTree> {
    using Base = LocalIndexSetBase<LocalBasisTree>;
  public:

    using Element = typename LFS::Traits::Element;

    LocalView(const Basis& basis)
      : Base{basis, makeLocalBasisTree<fast, LFSCache>(*basis._gfs, TypeTree::treePath())}
      , _lfs{ std::make_unique<LFS>(*basis._gfs) }
      , _cache{ std::make_unique<LFSCache>(*_lfs) }
    {
      PDELab::forEachLeafNode(*this->_ltree_storage, [&](auto& leaf){
        leaf.setCacheView(_cache.get());
      });
    }

    LocalView(const LocalView& localview)
      : LocalView( localview.globalBasis() )
    {}

    LocalView(LocalView&& other)
      : Base( std::move(other) )
      , _lfs(move(other._lfs))
      , _cache(move(other._cache))
      , _entity_view(other._entity_view)
      , _entity_storage(move(other._entity_storage))
    {}

    [[nodiscard]] size_type size() const {
      return _lfs->size();
    };


    [[nodiscard]] MultiIndex index(size_type dof) const noexcept {
      if constexpr (fast and not DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE_LEGACY)
        DUNE_THROW(NotImplemented, "To enable this feature complie with `-DDUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE_LEGACY`");
      else if constexpr (fast)
        return this->_indices_view[dof];
      else
        return _cache->containerIndex(dof);
    }

    template<std::convertible_to<Element> E>
    LocalView& bind(E&& element) {
      bindElement(std::forward<E>(element));
      _lfs->bind(this->element());
      _cache->update();
      Base::bind(this->element());
      return *this;
    }

    LocalView& unbind() {
      Base::unbind();
      _entity_view = nullptr;
      _entity_storage = std::nullopt;
      return *this;
    }

    template<std::convertible_to<Element> E>
    friend void bind(E&& element, LocalView& lbasis0, auto& lbasis1) {
      lbasis0.bind(std::forward<E>(element));
      lbasis1.bind(lbasis0.element());
    }

    friend void unbind(LocalView& lbasis0, auto& lbasis1) {
      lbasis1.unbind();
      lbasis0.unbind();
    }

#if DUNE_PDELAB_ENABLE_DOUBLE_BIND
    template<std::convertible_to<Element> E>
    friend void bind(E&& element, LocalView& lbasis0, LocalView& lbasis1) {
      lbasis0.bindElement(std::forward<E>(element));
      lbasis0._lfs->bind(lbasis0.element());
      lbasis0._cache->update();
      lbasis1.bindElement(lbasis0.element());
      lbasis0.doubleBind(lbasis0.element(), lbasis1);
    }

    friend void unbind(LocalView& lbasis0, LocalView& lbasis1) {
      lbasis0.doubleUnbind(lbasis1);
      lbasis0._entity_view = lbasis1._entity_view = nullptr;
      lbasis0._entity_storage = std::nullopt;
      lbasis1._entity_storage = std::nullopt;
    }
#endif

    [[nodiscard]] const Element& element() const noexcept {
      assert(_entity_view);
      return *_entity_view;
    };

  private:
    void bindElement(Element&& element) {
      // the caller assigned the ownership to us
      _entity_storage.emplace(std::move(element));
      _entity_view = &(*_entity_storage);
    }

    void bindElement(const Element& element) {
        // ownership is managed elsewhere
      _entity_view = &element;
    }

    std::unique_ptr<LFS> _lfs;
    std::unique_ptr<LFSCache> _cache;
    Element const* _entity_view;
    std::optional<Element> _entity_storage;
  };


  [[nodiscard]] LocalView localView() const { return LocalView{ *this }; }

  [[nodiscard]] LocalIndexSet localIndexSet() const { return LocalIndexSet{ *this }; }

  [[nodiscard]] LocalConstraints localConstraints() const {
    return _constraints_container->localView(localView().tree(),TypeTree::treePath());
  }

  [[nodiscard]] const EntitySet& entitySet() const noexcept { return _gfs->entitySet(); }

  [[nodiscard]] const size_type dimension() const noexcept { return _gfs->size(); }

  [[nodiscard]] auto degree() const noexcept {
    return _gfs->degree();
  }

  [[nodiscard]] size_type size(const SizePrefix& prefix) const noexcept
  {
    return _gfs->ordering().size(reverse(prefix));
  }

  [[nodiscard]] size_type size() const noexcept { return _gfs->size(); }

  [[nodiscard]] bool fixedSize(std::size_t dim, std::size_t codim) const noexcept { return _gfs->dataHandleFixedSize(codim); }

  [[nodiscard]] bool contains(std::size_t dim, std::size_t codim) const noexcept { return _gfs->dataHandleContains(codim); }

  // Whether local index sets match in all processors (constraints may still differ!sett)
  [[nodiscard]] auto isLocalIndexSetConforming() const noexcept {
    return *_conforming_local_index_set;
  }

  void update(const EntitySet& entity_set) {
    _gfs->update(entity_set);
    updateConstraints();

    MultipleCodimMultipleGeomTypeMapper mapper{entitySet(), [](auto, auto){return 1;}};

    // technically we don't need the container, but it's easier to reuse the entity data handle
    std::vector<std::size_t> mismatching_sizes(mapper.size(), 0);
    // notice that we always calculate this with the root node.
    // this allows us to share and reuse the result to any sub-basis
    auto local_index_set = localIndexSet();

    // // communicate local index set sizes and compare them at receiving end
    // applyToDistributedEntities(
    //   mapper,
    //   mismatching_sizes,
    //   Dune::InterfaceType::All_All_Interface,
    //   [&](auto phase,
    //       const auto& entity,
    //       auto& remote_value,
    //       auto& local_value) {
    //     local_index_set.bind(entity);
    //     if constexpr (decltype(phase)::value == Communication::gather)
    //       remote_value = local_index_set.size();
    //     else
    //       local_value = (local_index_set.size() != remote_value);
    //   });

    // // accumulate number of entities that mismatch
    // auto missmatching =
    //   std::accumulate(begin(mismatching_sizes), end(mismatching_sizes), 0);
    // *_conforming_local_index_set = (missmatching == 0);
  }


  [[nodiscard]] friend bool operator==(const Basis&, const Basis&) = default;
  [[nodiscard]] friend bool operator!=(const Basis&, const Basis&) = default;

private:
  void updateConstraints() {
    _constraints_container->assembleConstraints(*this, TypeTree::makeTreeContainer(*_constraints_tree, [](auto node){return node;}));
  }

  std::shared_ptr<GridFunctionSpace> _gfs;
  std::shared_ptr<RootConstraintsContainer> _constraints_container;
  std::shared_ptr<ConstraintsTree> _constraints_tree;
  std::shared_ptr<bool> _conforming_local_index_set;
};

} // namespace Dune::PDELab::inline Experimental::Legacy

#endif // DUNE_PDELAB_BASIS_WRAPPER_GRIDFUNCTIONSPACE_HH
