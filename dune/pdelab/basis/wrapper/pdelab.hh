#ifndef DUNE_ASSEMBLER_DISCRETE_FUNCTION_SPACE_WRAPPER_PDELAB_HH
#define DUNE_ASSEMBLER_DISCRETE_FUNCTION_SPACE_WRAPPER_PDELAB_HH

#include <dune/assembler/concepts/space.hh>
#include <dune/assembler/concepts/indexable.hh>

#include <dune/assembler/common/multiindex.hh>
#include <dune/assembler/common/reservedmultiindex.hh>
#include <dune/assembler/common/entityset.hh>
#include <dune/assembler/common/communication/entity_data_handler.hh>

#include <dune/assembler/space/constraints/container.hh>

#if !HAVE_DUNE_PDELAB
#error "This header is only available if PDELab headers are found"
#endif


#include <dune/pdelab/gridfunctionspace/gridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>

#include <dune/typetree/treecontainer.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <array>

#ifndef DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE_PDELAB
#define DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE_PDELAB 0
#endif

#ifndef DUNE_ASSEMBLER_ENABLE_DOUBLE_BIND
#define DUNE_ASSEMBLER_ENABLE_DOUBLE_BIND 1
#endif


namespace Dune::Assembler::PDELab {

template<class Cache, class TreePath, bool fast>
class LeafLocaIndexSet : public TypeTree::LeafNode
{
  using CI = typename Cache::ContainerIndex;

public:
  LeafLocaIndexSet(TreePath path)
    : _path{ path }
  {}

  using Path = TreePath;
  using size_type = typename Cache::size_type;
  using MultiIndex = ReservedMultiIndex<typename CI::value_type, CI::max_size()>;

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
    CI source = _cache_view->containerIndex(this->localIndex(fast ? 0 : dof));
    if constexpr (fast) source.front() += dof;
    // dune-assembler assumes a reversed order of indices wrt dune-pdelab.
    // Reversing them is quite inefficient. That's why we provide a friend
    // function to index the container for then native handler.
    MultiIndex target;
    if constexpr (CI::max_size() == 1)
      target[0] = source[0];
    else {
      target.resize(source.size());
      std::reverse_copy(source.begin(), source.end(), target.begin());
    }
    return target;
  }

  [[nodiscard]] friend decltype(auto) localContainerEntry(
                           auto& container,
                           const LeafLocaIndexSet& lspace,
                           size_type dof)
  {
    CI ci = lspace._cache_view->containerIndex(lspace.localIndex(fast ? 0 : dof));
    if constexpr (fast) ci.front() += dof;
    if constexpr (requires {container[ci];}) // PDELab backend
      return container[ci];
    else // standard indexable container
      return localContainerEntryImpl(container, MultiIndex{ci});
  }

  void setCacheView(Cache const * cache_view) {
    _cache_view = cache_view;
  }

private:

  static decltype(auto) localContainerEntryImpl(auto&& container, Concept::MultiIndex auto prefix)
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
class LeafLocalSpace : public TypeTree::LeafNode
{
  using LFS = typename Cache::LocalFunctionSpace;
  using LeafLFS = TypeTree::ChildForTreePath<LFS, TreePath>;
  using CI = typename Cache::ContainerIndex;

public:
  LeafLocalSpace(TreePath path)
    : _path{ path }
  {}

  using Path = TreePath;
  using Element = typename LeafLFS::Traits::Element;
  using FiniteElement = typename LeafLFS::Traits::FiniteElement;
  using size_type = typename LeafLFS::Traits::SizeType;
  using MultiIndex = ReservedMultiIndex<typename CI::value_type, CI::max_size()>;

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
    CI source;
    if constexpr (fast) {
      source = _cache_view->containerIndex(this->localIndex(0));
      source.front() += dof;
    } else {
      source =_cache_view->containerIndex(this->localIndex(dof));
    }
    // dune-assembler assumes a reversed order of indices wrt dune-pdelab.
    // Reversing them is quite inefficient. That's why we provide a friend
    // function to index the container for then native handler.
    MultiIndex target;
    if constexpr (CI::max_size() == 1)
      target[0] = source[0];
    else {
      target.resize(source.size());
      std::reverse_copy(source.begin(), source.end(), target.begin());
    }
    return target;
  }

  [[nodiscard]] friend decltype(auto) localContainerEntry(
                           auto& container,
                           const LeafLocalSpace& lspace,
                           size_type dof) noexcept
  {
    CI ci;
    if constexpr (fast) {
      ci = lspace._cache_view->containerIndex(lspace.localIndex(0));
      ci.front() += dof;
    } else {
      ci = lspace._cache_view->containerIndex(lspace.localIndex(dof));
    }
    if constexpr (requires {container[ci];}) // PDELab backend
      return container[ci];
    else // standard indexable container
      return localContainerEntryImpl(container, ci, ci.size()-1 );
  }

private:

  static decltype(auto) localContainerEntryImpl(auto&& container, const auto& prefix, std::size_t i) noexcept
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
class ArrayLocalSpace : public TypeTree::PowerNode<T, k>
{
  using TreeNode = TypeTree::PowerNode<T, k>;

public:
  ArrayLocalSpace(const typename TreeNode::NodeStorage& storage)
    : TreeNode{ storage }
  {
  }
};

template<class... T>
class TupleLocalSpace : public TypeTree::CompositeNode<T...>
{
  using TreeNode = TypeTree::CompositeNode<T...>;

public:
  TupleLocalSpace(const typename TreeNode::NodeStorage& storage)
    : TreeNode{ storage }
  {
  }
};

template<bool fast, class Cache, class GFSNode, class TreePath>
auto makeLocalSpaceTree(const GFSNode& gfs_node, TreePath path)
{
  if constexpr (GFSNode::isLeaf) {
    using LocalView = LeafLocalSpace<Cache, TreePath, fast>;
    return std::make_unique<LocalView>(path);
  } else if constexpr (GFSNode::isPower) {
    using ChidlNode = std::decay_t<decltype(*makeLocalSpaceTree<fast,Cache>(
      gfs_node.child(0), push_back(path, 0)))>;
    using LocalView = ArrayLocalSpace<ChidlNode, GFSNode::degree()>;
    typename LocalView::NodeStorage storage;
    for (std::size_t i = 0; i < gfs_node.degree(); ++i)
      storage[i] =
        makeLocalSpaceTree<fast,Cache>(gfs_node.child(i), push_back(path, i));
    return std::make_unique<LocalView>(storage);
  } else {
    static_assert(GFSNode::isComposite);
    return unpackIntegerSequence(
      [&](auto... i) {
        using ChildNodes = std::tuple<std::decay_t<decltype(*makeLocalSpaceTree<fast,Cache>(
          gfs_node.child(i), push_back(path, i)))>...>;
        using LocalView =
          TupleLocalSpace<std::tuple_element_t<i, ChildNodes>...>;
        typename LocalView::NodeStorage storage{ makeLocalSpaceTree<fast, Cache>(
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
    using LocalView = ArrayLocalSpace<ChidlNode, GFSNode::degree()>;
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
          TupleLocalSpace<std::tuple_element_t<i, ChildNodes>...>;
        typename LocalView::NodeStorage storage{ makeLocalIndexSetTree<fast, Cache>(
          gfs_node.child(i), push_back(path, i))... };
        return std::make_unique<LocalView>(storage);
      },
      std::make_index_sequence<GFSNode::degree()>{});
  }
}

template<class GridFunctionSpace, Concept::Tree ConstraintsTree>
class Space
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

  using LocalSpaceTree = std::decay_t<decltype(*makeLocalSpaceTree<fast,LFSCache>(std::declval<GridFunctionSpace>(), multiIndex()))>;
  using LocalIndexSetTree = std::decay_t<decltype(*makeLocalIndexSetTree<fast,LISCache>(std::declval<GridFunctionSpace>(), multiIndex()))>;

  using MI = typename GridFunctionSpace::Ordering::Traits::ContainerIndex;
  using SP = typename GridFunctionSpace::Ordering::Traits::ContainerIndex;

public:
  using MultiIndex = ReservedMultiIndex<typename MI::value_type, MI::max_size()>;
  using SizePrefix = ReservedMultiIndex<typename SP::value_type, SP::max_size()>;
  using EntitySet = typename GridFunctionSpace::Traits::EntitySet;
  using size_type = typename GridFunctionSpace::Traits::SizeType;

private:
  static constexpr auto constraints_container_generator = []<Concept::LeafTreeNode LeafSpace, Concept::MultiIndex Path>(const LeafSpace& leaf_space, Path path){
    using Constraints = TypeTree::ChildForTreePath<ConstraintsTree, Path>;
    using ConstrainsContainer = typename Constraints::template Container<MultiIndex, typename LeafSpace::Traits::EntitySet>;
    return std::make_shared<ConstrainsContainer>(leaf_space.entitySet());
  };

  using RootConstraintsContainerStorage = decltype(makeConstraintsContainer(std::declval<const GridFunctionSpace&>(), constraints_container_generator));
  using RootConstraintsContainer = typename RootConstraintsContainerStorage::element_type;
public:

  using LocalConstraints = decltype(std::declval<RootConstraintsContainer>().localView(std::declval<LocalSpaceTree>(),multiIndex()));

  Space(std::shared_ptr<GridFunctionSpace> gfs_ptr, const ConstraintsTree& constraints_tree)
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

    using MultiIndex = ReservedMultiIndex<typename MI::value_type, MI::max_size()>;
    using Tree = LocalTree;
    using size_type = typename LISCache::size_type;
    using GlobalBasis = Space;

    LocalIndexSetBase(const Space& space, std::unique_ptr<LocalTree> ltree_storage)
      : _space{space}
      , _ltree_storage{move(ltree_storage)}
    {
        if constexpr (DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE_PDELAB)
          _indices = std::make_unique<MultiIndex[]>(maxSize());
    }

    LocalIndexSetBase(const LocalIndexSetBase&) = delete;

    LocalIndexSetBase(LocalIndexSetBase&& other)
      : _space{ std::move(other._space) }
      , _ltree_storage{ move(other._ltree_storage) }
      , _ltree_view{ other._ltree_view }
      , _mem_region{ other._mem_region }
    {
      if constexpr (fast and DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE_PDELAB)
        _indices = move(other._indices);
    }

    void bind(const Dune::Concept::Entity auto& entity) {
      std::size_t tree_index = 0;
      std::size_t tree_offset = 0;
      Assembler::forEachLeafNode(*_ltree_storage, [&](auto& leaf){
        leaf.setTreeOffset(tree_offset, tree_index);
        if constexpr (fast and DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE_PDELAB)
          for(size_type dof = 0; dof != leaf.size(); ++dof)
            _indices[tree_offset + dof] = leaf.index(dof);
        tree_index  += 1;
        tree_offset += leaf.size();
      });
      assert(tree_offset <= maxSize());
      if constexpr (fast and DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE_PDELAB)
        _indices_view = _indices.get();
      _ltree_view = _ltree_storage.get();
      _mem_region = _space.entitySet().memoryRegion(entity);
    }

    void unbind() {
      _ltree_view = nullptr;
      _indices_view = nullptr;
    }

    [[nodiscard]] std::convertible_to<MemoryRegion> auto memoryRegion() const noexcept {
      return _mem_region; // TODO!
    }

    [[nodiscard]] size_type maxSize() const noexcept {
      return _space._gfs->maxLocalSize();
    };

    [[nodiscard]] const Tree& tree() const noexcept {
      return *_ltree_view;
    }

    [[nodiscard]] const GlobalBasis& globalBasis() const noexcept {
      return _space;
    }

    // Whether local index sets match in all processors
    [[nodiscard]] auto conforming() const noexcept {
      return _space.isLocalIndexSetConforming();
    }

  protected:

    void doubleBind(const Dune::Concept::Entity auto& element, LocalIndexSetBase& other) {
      bind(element);
      if (_space == other._space) {
        other._ltree_view = _ltree_view;
        other._indices_view = _indices_view;
        other._mem_region = _mem_region;
      } else {
        other.bind(element);
      }
    }

    void doubleUnbind(LocalIndexSetBase& other) {
      unbind();
      if (_space == other._space) {
        other._ltree_view = nullptr;
        other._indices_view = nullptr;
      } else {
        other.unbind();
      }
    }

    GlobalBasis _space;
    std::unique_ptr<LocalTree> _ltree_storage;
    LocalTree const * _ltree_view = nullptr;
    std::unique_ptr<MultiIndex[]> _indices;
    MultiIndex const * _indices_view = nullptr;
    MemoryRegion _mem_region;
  };


public:

  class LocalIndexSet : public LocalIndexSetBase<LocalIndexSetTree> {
    using Base = LocalIndexSetBase<LocalIndexSetTree>;
  public:

    LocalIndexSet(const Space& space)
      : Base{space, makeLocalIndexSetTree<fast,LISCache>(*space._gfs, multiIndex())}
      , _cache{ std::make_unique<LISCache>(*space._gfs) }
    {
      Assembler::forEachLeafNode(*this->_ltree_storage, [&](auto& leaf){
        leaf.setCacheView(_cache.get());
      });
    }

    LocalIndexSet(const LocalIndexSet& localview)
      : LocalIndexSet( localview._space )
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

    friend void bind(const Dune::Concept::Entity auto& element, LocalIndexSet& lspace0, auto& lspace1) {
      lspace1.bind(element);
      lspace0.bind(element);
    }

#if DUNE_ASSEMBLER_ENABLE_DOUBLE_BIND
    friend void bind(const Dune::Concept::Entity auto& element, LocalIndexSet& lspace0, LocalIndexSet& lspace1) {
      lspace0._cache->update(element);
      lspace0.doubleBind(element, lspace1);
    }

    friend void unbind(LocalIndexSet& lspace0, LocalIndexSet& lspace1) {
      lspace0.doubleUnbind(lspace1);
    }
#endif

    [[nodiscard]] size_type size() const noexcept {
      return _cache->size();
    }

    [[nodiscard]] MultiIndex index(size_type dof) const noexcept {
      if constexpr (fast and not DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE_PDELAB)
        DUNE_THROW(NotImplemented, "To enable this feature complie with `-DDUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE_PDELAB`");
      else if constexpr (fast)
        return this->_indices_view[dof];
      else
        return _cache->containerIndex(dof);
    }

  private:
    std::unique_ptr<LISCache> _cache;
  };

  class LocalView : public LocalIndexSetBase<LocalSpaceTree> {
    using Base = LocalIndexSetBase<LocalSpaceTree>;
  public:

    using Element = typename LFS::Traits::Element;

    LocalView(const Space& space)
      : Base{space, makeLocalSpaceTree<fast, LFSCache>(*space._gfs, multiIndex())}
      , _lfs{ std::make_unique<LFS>(*space._gfs) }
      , _cache{ std::make_unique<LFSCache>(*_lfs) }
    {
      Assembler::forEachLeafNode(*this->_ltree_storage, [&](auto& leaf){
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
      if constexpr (fast and not DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE_PDELAB)
        DUNE_THROW(NotImplemented, "To enable this feature complie with `-DDUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE_PDELAB`");
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
    friend void bind(E&& element, LocalView& lspace0, auto& lspace1) {
      lspace0.bind(std::forward<E>(element));
      lspace1.bind(lspace0.element());
    }

    friend void unbind(LocalView& lspace0, auto& lspace1) {
      lspace1.unbind();
      lspace0.unbind();
    }

#if DUNE_ASSEMBLER_ENABLE_DOUBLE_BIND
    template<std::convertible_to<Element> E>
    friend void bind(E&& element, LocalView& lspace0, LocalView& lspace1) {
      lspace0.bindElement(std::forward<E>(element));
      lspace0._lfs->bind(lspace0.element());
      lspace0._cache->update();
      lspace1.bindElement(lspace0.element());
      lspace0.doubleBind(lspace0.element(), lspace1);
    }

    friend void unbind(LocalView& lspace0, LocalView& lspace1) {
      lspace0.doubleUnbind(lspace1);
      lspace0._entity_view = lspace1._entity_view = nullptr;
      lspace0._entity_storage = std::nullopt;
      lspace1._entity_storage = std::nullopt;
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
    return _constraints_container->localView(localView().tree(),multiIndex());
  }

  [[nodiscard]] const EntitySet& entitySet() const noexcept { return _gfs->entitySet(); }

  [[nodiscard]] const size_type dimension() const noexcept { return _gfs->size(); }

  [[nodiscard]] auto degree() const noexcept {
    return _gfs->degree();
  }

  [[nodiscard]] size_type size(const SizePrefix& prefix) const noexcept
  {
    MI pdelab_suffix;
    pdelab_suffix.resize(prefix.size());
    std::reverse_copy(prefix.begin(), prefix.end(),pdelab_suffix.begin());
    return _gfs->ordering().size(pdelab_suffix);
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
    // this allows us to share and reuse the result to any sub-space
    auto local_index_set = localIndexSet();

    // communicate local index set sizes and compare them at receiving end
    applyToDistributedEntities(
      mapper,
      mismatching_sizes,
      Dune::InterfaceType::All_All_Interface,
      [&](auto phase,
          const auto& entity,
          auto& remote_value,
          auto& local_value) {
        local_index_set.bind(entity);
        if constexpr (decltype(phase)::value == Communication::gather)
          remote_value = local_index_set.size();
        else
          local_value = (local_index_set.size() != remote_value);
      });

    // accumulate number of entities that mismatch
    auto missmatching =
      std::accumulate(begin(mismatching_sizes), end(mismatching_sizes), 0);
    *_conforming_local_index_set = (missmatching == 0);
  }


  [[nodiscard]] friend bool operator==(const Space&, const Space&) = default;
  [[nodiscard]] friend bool operator!=(const Space&, const Space&) = default;

private:
  void updateConstraints() {
    _constraints_container->assembleConstraints(*this, TypeTree::makeTreeContainer(*_constraints_tree, [](auto node){return node;}));
  }

  std::shared_ptr<GridFunctionSpace> _gfs;
  std::shared_ptr<RootConstraintsContainer> _constraints_container;
  std::shared_ptr<ConstraintsTree> _constraints_tree;
  std::shared_ptr<bool> _conforming_local_index_set;
};

} // namespace Dune::Assembler::PDELab

#endif // DUNE_ASSEMBLER_DISCRETE_FUNCTION_SPACE_WRAPPER_PDELAB_HH
