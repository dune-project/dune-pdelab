#ifndef DUNE_PDELAB_BASIS_ORDERING_ENTITY_LEAF_HH
#define DUNE_PDELAB_BASIS_ORDERING_ENTITY_LEAF_HH

#include <dune/pdelab/basis/prebasis/concept.hh>
#include <dune/pdelab/basis/ordering/entity_base.hh>

#include <dune/pdelab/common/multiindex.hh>
#include <dune/pdelab/common/tree_traversal.hh>
#include <dune/pdelab/common/container_entry.hh>

#include <dune/pdelab/concepts/multiindex.hh>
#include <dune/pdelab/concepts/treenode.hh>

#include <dune/typetree/leafnode.hh>
#include <dune/typetree/dynamicpowernode.hh>
#include <dune/typetree/treepath.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>
#include <dune/geometry/referenceelements.hh>

#include <vector>
#include <optional>
#include <numeric>

#ifndef DUNE_PDELAB_FAST_DG_OPTIMIZATION
#define DUNE_PDELAB_FAST_DG_OPTIMIZATION 1
#endif

namespace Dune::PDELab::inline Experimental::Impl {

/**
 * @brief Leaf implementation of an EntityOrderingNode
 * @note This is the first ordering constructed on every ordering tree
 *
 * @tparam PreBasisLeaf   Type containing a finite element map and a
 *                        merging strategy
 */
template<Concept::Impl::PreBasisLeaf PreBasisLeaf>
class LeafEntityOrdering
  : public TypeTree::LeafNode
  , public EntityOrderingNode<LeafEntityOrdering<PreBasisLeaf>,
                              typename PreBasisLeaf::Traits::MergingStrategy>
{
  using TreeNode = TypeTree::LeafNode;
  using OrderingNode =
    EntityOrderingNode<LeafEntityOrdering<PreBasisLeaf>,
                       typename PreBasisLeaf::Traits::MergingStrategy>;
  using FE = typename PreBasisLeaf::Traits::FiniteElementMap::Traits::FiniteElement;
  using ES = typename PreBasisLeaf::Traits::MergingStrategy::EntitySet;
  static constexpr std::size_t fem_dim = FE::Traits::LocalBasisType::Traits::dimDomain;
  static constexpr std::size_t fem_codim = ES::dimension - fem_dim;
public:
  using Space = PreBasisLeaf;
  using SizeType = typename PreBasisLeaf::Traits::MergingStrategy::SizeType;

  //! Construct an ordering based on the templated discrete-function-space
  LeafEntityOrdering(const Space& space)
    : TreeNode{}
    , OrderingNode{ space.mergingStrategy() }
    , _space{ space }
  {
  }

  LeafEntityOrdering(const LeafEntityOrdering&) = delete;
  LeafEntityOrdering(LeafEntityOrdering&&) = default;

  LeafEntityOrdering& operator=(const LeafEntityOrdering&) = delete;
  LeafEntityOrdering& operator=(LeafEntityOrdering&&) = default;

  //! Access to the underlying discrete-function-space
  const Space& space() const { return _space; }

  //! Access to the underlying discrete-function-space
  Space& space() { return _space; }

  [[nodiscard]] int maxSubEntities() const { return _max_sub_entities; }
  void setMaxSubEntities(int s) { _max_sub_entities = s; }

  /**
   * @brief Coefficient vector factory
   *
   * @tparam BackendTraits    A traits class providing container types
   * @return constexpr auto   a suitable container type for the underlying
   * ordering
   */
  template<class BackendTraits>
  static constexpr auto makeVectorContainer()
  {
    // common size for all geometry types
    constexpr auto gt_common_size = commonSizePerGeometryType();

    auto field = BackendTraits::template makeField<Space>();
    using Field = std::decay_t<decltype(field)>;
    if constexpr (gt_common_size)
      return BackendTraits::template makeArray<Field, gt_common_size.value()>();
    else
      return BackendTraits::template makeVector<Field>();
  }

  //! gets a common size for all active geometry types if available at compile time
  // this needs the the Space::Traits::FiniteElementMap::size(...) is static constexpr.
  static consteval std::optional<std::size_t> commonSizePerGeometryType() {
    using FEM = typename Space::Traits::FiniteElementMap;
    std::size_t size = 0;
    const auto fem_size = [](auto gt) consteval {
      if constexpr (requires { { FEM::size(gt) } -> std::convertible_to<std::size_t>; })
        return FEM::size(gt);
      else
        return 0;
    };
    // iterate over all possible geometry types and find out if all share the same size
    for (std::size_t dim = 0; dim <= FEM::dimension ; ++dim) {
      std::size_t gt_size = fem_size(GeometryTypes::none(dim));
      if (gt_size > 0) {
        if (size > 0 and size != gt_size)
          return std::nullopt;
        else
          size = gt_size;
      }
      for (std::size_t topology_id = 0 ; topology_id < (std::size_t{1} << dim) ; ++topology_id) {
        std::size_t gt_size = fem_size(GeometryType(topology_id,dim));
        if (gt_size > 0) {
          if (size > 0 and size != gt_size)
            return std::nullopt;
          else
            size = gt_size;
        }
      }
    }
    if (size == 0)
      return std::nullopt;
    else
      return size;
  }

private:

  template<Concept::MultiIndex ViewPath, Concept::MultiIndex ContainerIndex>
  class LocalIndexSet : public TypeTree::LeafNode
  {
    using DisjointCodimClosure = decltype(std::declval<LeafEntityOrdering>().disjointCodimClosure());
  public:
    using size_type = SizeType;
    using MultiIndex = ContainerIndex;
    using Path = ViewPath;

    //! Class constructing the local space
    LocalIndexSet(const Path& view_path, const DisjointCodimClosure& disjoint_codim_closure)
      : TypeTree::LeafNode{}
      , _indices{}
      , _tree_offset{ 0 }
      , _view_path{view_path}
      , _disjoint_codim_closure{disjoint_codim_closure}
    {
    }

    LocalIndexSet(const LocalIndexSet&) = delete;
    LocalIndexSet(LocalIndexSet&&) = default;

    LocalIndexSet& operator=(const LocalIndexSet&) = delete;
    LocalIndexSet& operator=(LocalIndexSet&&) = default;

    //! Grant mutable access to the underlying local container indices
    std::vector<MultiIndex>& indices() noexcept { return _indices; }

    //! Obtain a container index (Inner2Outer) to the i-th local degree of freedom
    [[nodiscard]] MultiIndex index(size_type dof) const noexcept
    {
      assert(dof < size());
      if (optimizeFastDG() and disjointCodimClosure()) {
        // In this case, we know that indices are contiguous and can be
        //   reconstructed from the local index so there is no need to fill
        //   the whole vector.
        return accumulate_back(_indices[0], dof);
      } else {
        return _indices[dof];
      }
    }


    void unbind() noexcept {
      _indices.clear();
#ifndef NDEBUG
      _tree_offset = std::numeric_limits<std::size_t>::max()/2;
#endif
    }

    // returns a local view path: join(orderingViewPath(), subEntityPath())
    [[nodiscard]] Path path() const noexcept {
      return _view_path;
    }

    [[nodiscard]] auto orderingViewPath() const noexcept {
      if constexpr (fem_dim != ES::dimension)
        return pop_back(_view_path);
      else
        return _view_path;
    }

    [[nodiscard]] auto subEntityPath() const noexcept {
      if constexpr (fem_dim != ES::dimension)
        return back(_view_path);
      else
        return TypeTree::treePath();
    }

    //! Sets an offest for the local degrees of freedom
    void setTreeOffset(size_type tree_offset) noexcept { _tree_offset = tree_offset; }

    //! Returns unique index whitin the tree dofs
    [[nodiscard]] size_type localIndex(size_type dof) const noexcept
    {
      return _tree_offset + dof;
    }

    //! Amount of degrees of freedom for this node
    [[nodiscard]] size_type size() const noexcept { return _indices.size(); }

    //! Bool whether dofs are shared to other entity local spaces
    [[nodiscard]] auto disjointCodimClosure() const noexcept
    {
      return _disjoint_codim_closure;
    }

    [[nodiscard]] static constexpr bool optimizeFastDG() noexcept
    {
      return DUNE_PDELAB_FAST_DG_OPTIMIZATION;
    }

    [[nodiscard]] friend constexpr decltype(auto) localContainerEntry(auto& container, const LocalIndexSet& node, size_type dof) noexcept {
      return containerEntry(container, node.index(dof));
    }

    template<class OtherLocalIndexSet>
    [[nodiscard]] friend decltype(auto) localContainerEntry(auto&& container,
                                         const LocalIndexSet& test_lis,
                                         size_type test_dof,
                                         const OtherLocalIndexSet& trial_lis,
                                         size_type trial_dof) noexcept
    {
      static_assert(MultiIndex::size() == OtherLocalIndexSet::MultiIndex::size());
      MultiIndex i = test_lis.index(test_dof);
      // TODO check that `trial_lis` has the same multi-index direction semantics
      Concept::MultiIndex auto j = trial_lis.index(trial_dof);
      Concept::MultiIndex auto k = containerIndexMerge(container, j, i);
      return containerEntry(container, k);
    }


  private:

    std::vector<MultiIndex> _indices;
    size_type _tree_offset;
    // the following objects may be compile-time or run-time values, in the
    // first case (i.e. MultiIndex<> and std::bool_constant<...>), they may be
    // completely empty thus we do not rely on their addresses and try make data
    // layout of this class smaller.
    [[no_unique_address]] const Path _view_path;
    [[no_unique_address]] const DisjointCodimClosure _disjoint_codim_closure;
  };

  template<Concept::MultiIndex ViewPath, Concept::MultiIndex ContainerIndex>
  class SubEntityLocalIndexSet : public TypeTree::DynamicPowerNode<LocalIndexSet<ViewPath,ContainerIndex>>
  {
    using Base = TypeTree::DynamicPowerNode<LocalIndexSet<ViewPath,ContainerIndex>>;
  public:
    using size_type = SizeType;
    using MultiIndex = ContainerIndex;
    SubEntityLocalIndexSet(ViewPath view_path, typename Base::NodeStorage&& storage)
      : Base{ std::move(storage) }
      , _view_path{ view_path }
    {
    }

    SubEntityLocalIndexSet(const SubEntityLocalIndexSet&) = delete;
    SubEntityLocalIndexSet(SubEntityLocalIndexSet&&) = default;

    [[nodiscard]] auto orderingViewPath() const noexcept {
        return _view_path;
    }

  private:
    [[no_unique_address]] const ViewPath _view_path;
  };


  /**
   * @brief A realization of a local space node
   * This class provides an ordering to the space induced by the
   * discrete-function-space but restricted to a particular entity. Countrary to
   * LeafEntityOrdering, this class is only instantiated after the whole
   * ordering tree has been constructed.
   *
   * @note This class fullfils the Concept::LocalViewLeaf, thus, is the
   * facade class for the user.
   * @note ContainerIndex usually differs between different nodes of the same
   * ordering. In particular, it is an hybrid multi-index where known parent
   * indices for this node may be pre-filled at compile-time.
   *
   * @tparam ContainerIndex  The type to store container indices
   */
  template<Concept::MultiIndex ViewPath, Concept::MultiIndex ContainerIndex>
  class LocalView : public LocalIndexSet<ViewPath, ContainerIndex>
  {
    using DisjointCodimClosure = decltype(std::declval<LeafEntityOrdering>().disjointCodimClosure());
  public:
    using FiniteElement = FE;
    using EntitySet = ES;
    using Element = typename EntitySet::template Codim<fem_codim>::Entity;
    using size_type = SizeType;
    using MultiIndex = ContainerIndex;
    using Path = ViewPath;
    static constexpr std::size_t dimDomain = fem_dim;

    //! Local space traits (compatiblility layer for PDELab)
    struct Traits
    {
      using FiniteElementType /*[[deprecated]]*/ = typename Space::Traits::
        FiniteElementMap::Traits::FiniteElement;
      using FiniteElement /*[[deprecated]]*/ = typename Space::Traits::
        FiniteElementMap::Traits::FiniteElement;
      using SizeType /*[[deprecated]]*/ = typename PreBasisLeaf::Traits::MergingStrategy::SizeType;
    };

    //! Class constructing the local space
    LocalView(const Path& view_path, const DisjointCodimClosure& disjoint_codim_closure)
      : LocalIndexSet<ViewPath, ContainerIndex>{ view_path, disjoint_codim_closure }
      , _entity_view{ nullptr }
      , _fe_view{ nullptr }
    {
    }

    LocalView(const LocalView&) = delete;
    LocalView(LocalView&&) = default;

    LocalView& operator=(const LocalView&) = delete;
    LocalView& operator=(LocalView&&) = default;


    //! Binds a finite element: rvalues are stored, lvalues are referenced
    template<class FE>
    void bindFiniteElement(FE&& finite_element) noexcept
    {
      static_assert(std::same_as<std::decay_t<FE>, FiniteElement>);
      if constexpr (std::is_rvalue_reference_v<FE&&>) {
        static_assert(std::move_constructible<FiniteElement>);
        static_assert(std::is_move_assignable_v<FiniteElement>);
        // we will pay the price for unique_ptr ctor on our first bind
        if (_fe_store) [[likely]]
          (*_fe_store) = std::move(finite_element);
        else
          _fe_store =
            std::make_unique<FiniteElement>(std::move(finite_element));
        _fe_view = _fe_store.get();
      } else {
        _fe_view = &finite_element;
      }
    }

    //! Binds a view on the entity. Internally, we keep a reference the object
    void bindElement(const Element* entity) noexcept {
      _entity_view = entity;
    }

    void unbind() noexcept {
      LocalIndexSet<ViewPath, ContainerIndex>::unbind();
      _fe_view = nullptr;
      _entity_view = nullptr;
    }

    //! Returns a view on the local finite element
    [[nodiscard]] const FiniteElement& finiteElement() const noexcept
    {
      assert(_fe_view);
      return *_fe_view;
    }

    //! Returns a view on the bound entity
    [[nodiscard]] const Element& element() const noexcept
    {
      assert(boundElement() && "Entity is not bound: local function is not bound or it has no support on the bound entity");
      return *_entity_view;
    }

    [[nodiscard]] bool boundElement() const noexcept {
      return _entity_view != nullptr;
    }

  private:
    std::unique_ptr<FiniteElement> _fe_store;
    Element const* _entity_view;
    FiniteElement const* _fe_view;
  };

  template<Concept::MultiIndex ViewPath, Concept::MultiIndex ContainerIndex>
  class SubEntityLocalView : public TypeTree::DynamicPowerNode<LocalView<ViewPath, ContainerIndex>> {
    using LV = LocalView<ViewPath, ContainerIndex>;
    using Base = TypeTree::DynamicPowerNode<LocalView<ViewPath, ContainerIndex>>;
  public:

    using size_type = SizeType;
    using MultiIndex = ContainerIndex;

    SubEntityLocalView(ViewPath view_path,
      const typename Base::NodeStorage&& storage)
      : Base{ std::move(storage) }
      , _view_path{ view_path }
    {
    }

    void bindElement(const auto& entity) noexcept {
      constexpr std::size_t fem_dim = LV::FiniteElement::Traits::LocalBasisType::Traits::dimDomain;
      constexpr std::size_t fem_codim = LV::EntitySet::dimension - fem_dim;
      _sub_entities.resize(entity->subEntities(fem_codim));
      // use emplace in order to use the move constructor (move assignement did not work on multidomainggrid+UGGrid)
      for (std::size_t s = 0; s != _sub_entities.size(); ++s) {
        _sub_entities[s].emplace(entity->template subEntity<fem_codim>(s));
        this->child(s).bindElement(&_sub_entities[s].value());
      }
      assert(_sub_entities.size() <= Base::degree());
    }

    [[nodiscard]] size_type size() const noexcept {
      size_type _size = 0;
      for (std::size_t i = 0; i != _sub_entities.size(); ++i)
        _size += this->child(i).size();
      return _size;
    }

    [[nodiscard]] auto orderingViewPath() const noexcept {
        return _view_path;
    }

  private:
    std::vector<std::optional<typename LV::Element>> _sub_entities;
    [[no_unique_address]] const ViewPath _view_path;

  };

public:
  /**
   * @brief Creates a local ordering pointer based on a root global ordering
   *
   * @tparam Ordering   The root entity set ordering (it contains this object as
   * a node)
   * @tparam Prefix     A prefix container tree to reach this node from the root
   * ordering
   * @param ordering    The root entity set ordering (it contains this object as
   * a node)
   * @param prefix      A prefix container tree to reach this node from the root
   * ordering
   * @return auto       A shared pointer for a local ordering
   */
  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubSpacePath>
  auto makeLocalIndexSet(const std::shared_ptr<Ordering>& ordering,
                         const Prefix& prefix, const SubSpacePath& sub_space_path) const
  {
    static_assert(SubSpacePath::size() <= Prefix::size(), "Sub space path cannot be larger than the total path from root");

    constexpr std::size_t sub_space_depth = SubSpacePath::size();
    constexpr std::size_t view_depth = Prefix::size() - sub_space_depth;
    auto view_path = unpackIntegerSequence([&](auto... i){
      return TypeTree::treePath(prefix[index_constant<sub_space_depth+i>{}] ...);
    }, std::make_index_sequence<view_depth>{});

    assert(&child(*ordering, prefix) == this &&
           "This class should be a node of the root ordering");
    // From the ordering, we extract the type of the container index for this node
    using ContainerIndex =
      decltype(ordering->firstContainerIndex(prefix, SizeType{}, SizeType{}));
    // and use it to instantiate the local ordering
    using FEM = typename Space::Traits::FiniteElementMap;
    using EntitySet = typename Space::Traits::MergingStrategy::EntitySet;
    constexpr std::size_t fem_dim = FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits::dimDomain;
    if constexpr (fem_dim == EntitySet::dimension) {
      return std::make_unique<LocalIndexSet<decltype(view_path), ContainerIndex>>(view_path, this->disjointCodimClosure());
    } else {
      using LIS = LocalIndexSet<decltype(push_back(view_path,0)), ContainerIndex>;
      using SLIS = SubEntityLocalIndexSet<decltype(push_back(view_path,0)), ContainerIndex>;
      typename SLIS::NodeStorage storage(_max_sub_entities);
      for (std::size_t sub_entity = 0; sub_entity != storage.size(); ++sub_entity)
        storage[sub_entity] = std::make_unique<LIS>(push_back(view_path,sub_entity), this->disjointCodimClosure());
      return std::make_unique<SLIS>(view_path, std::move(storage));
    }
  }

  template<class Ordering, Concept::MultiIndex Prefix, Concept::MultiIndex SubSpacePath>
  auto makeLocalView(const std::shared_ptr<Ordering>& ordering,
                     const Prefix& prefix, SubSpacePath) const
  {
    static_assert(SubSpacePath::size() <= Prefix::size(), "Sub space path cannot be larger than the total path from root");

    constexpr std::size_t sub_space_depth = SubSpacePath::size();
    constexpr std::size_t view_depth = Prefix::size() - sub_space_depth;
    auto view_path = unpackIntegerSequence([&](auto... i){
      return TypeTree::treePath(prefix[index_constant<sub_space_depth+i>{}] ...);
    }, std::make_index_sequence<view_depth>{});

    assert(&child(*ordering, prefix) == this &&
           "This class should be a node of the root ordering");
    // From the ordering, we extract the type of the container index for this
    // node
    using ContainerIndex = decltype(ordering->firstContainerIndex(prefix, SizeType{}, SizeType{}));
    using FEM = typename Space::Traits::FiniteElementMap;
    using EntitySet = typename Space::Traits::MergingStrategy::EntitySet;
    constexpr std::size_t fem_dim = FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits::dimDomain;
    // and use it to instantiate the local ordering
    if constexpr (fem_dim == EntitySet::dimension) {
      return std::make_unique<LocalView<decltype(view_path), ContainerIndex>>(view_path, this->disjointCodimClosure());
    } else {
      using LV = LocalView<decltype(push_back(view_path,0)), ContainerIndex>;
      using SLV = SubEntityLocalView<decltype(push_back(view_path,0)), ContainerIndex>;
      typename SLV::NodeStorage storage(_max_sub_entities);
      for (std::size_t sub_entity = 0; sub_entity != storage.size(); ++sub_entity)
        storage[sub_entity] = std::make_unique<LV>(push_back(view_path,sub_entity), this->disjointCodimClosure());
      return std::make_unique<SLV>(view_path, std::move(storage));
    }
  }

private:
  Space _space;
  int _max_sub_entities = 0;
};

} // namespace Dune::PDELab::inline Experimental::Impl

#endif // DUNE_PDELAB_BASIS_ORDERING_ENTITY_LEAF_HH
