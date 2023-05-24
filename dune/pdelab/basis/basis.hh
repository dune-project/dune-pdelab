#ifndef DUNE_PDELAB_BASIS_BASIS_HH
#define DUNE_PDELAB_BASIS_BASIS_HH

#include <dune/pdelab/concepts/multiindex.hh>
#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/concepts/local_basis.hh>

// #include <dune/pdelab/common/partition/identity.hh>
#include <dune/pdelab/common/multiindex.hh>
#include <dune/pdelab/common/container_entry.hh>
// #include <dune/pdelab/common/container_resize.hh>
// #include <dune/pdelab/common/communication/entity_data_handler.hh>

#include <dune/pdelab/basis/prebasis/concept.hh>
#include <dune/pdelab/basis/ordering.hh>
#include <dune/pdelab/basis/merging_strategy.hh>

#include <dune/pdelab/basis/constraints/container.hh>

#include <dune/pdelab/basis/prebasis/leaf.hh>
#include <dune/pdelab/basis/prebasis/composite.hh>

#include <dune/grid/common/mcmgmapper.hh>

#ifndef DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE
#define DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE 0
#endif

#ifndef DUNE_PDELAB_ENABLE_DOUBLE_BIND
#define DUNE_PDELAB_ENABLE_DOUBLE_BIND 1
#endif

namespace Dune::PDELab::inline Experimental {

  template<Dune::Concept::GridView ES, Concept::Impl::PreBasisTree PB, Concept::FixedSizeMultiIndex SubBasisPath = TypeTree::HybridTreePath<> >
  class Basis {
    // We need some types sub-basis independent so that we can shared them between templated sub-basiss
    using Ordering = std::decay_t<decltype(*Impl::makeOrdering(std::declval<PB>()))>;
    using LocalViewTree = std::decay_t<decltype(*std::declval<Ordering>().makeLocalView(std::shared_ptr<Ordering>{}, TypeTree::treePath(), SubBasisPath{}))>;
    using LocalIndexSetTree = std::decay_t<decltype(*std::declval<Ordering>().makeLocalIndexSet(std::shared_ptr<Ordering>{}, TypeTree::treePath(), SubBasisPath{}))>;
    static constexpr std::size_t ContainerDepth = Ordering::maxContainerDepth();
    static_assert(ContainerDepth > 0);

    using RootLocalViewTree = std::decay_t<decltype(*std::declval<Ordering>().makeLocalView(std::shared_ptr<Ordering>{}, TypeTree::treePath(), TypeTree::treePath()))>;

    struct ConstraintsContainerGenerator {
      template<Concept::Impl::PreBasisLeaf PreBasisLeaf, Concept::FixedSizeMultiIndex Path>
      auto operator()(const PreBasisLeaf& pre_basis_leaf, Path) const {
        using EntitySet = PreBasisLeaf::Traits::MergingStrategy::EntitySet;
        using MultiIndex = typename TypeTree::ChildForTreePath<RootLocalViewTree,Path>::MultiIndex;
        using ConstraintsContainer = PreBasisLeaf::Traits::ConstraintsOperator::template Container<MultiIndex, EntitySet>;
        return std::make_shared<ConstraintsContainer>(pre_basis_leaf.mergingStrategy().entitySet());
      }
    };

    using RootConstraintsContainer = std::decay_t<decltype(*makeConstraintsContainer(std::declval<const PB&>(), ConstraintsContainerGenerator{}))>;

  public:

    template<class Backend>
    using Container = std::decay_t<decltype(Ordering::template makeVectorContainer<Backend>())>;

    using size_type = std::size_t;
    // using EntitySetPartition = ESP;
    using EntitySet = ES;
    using MultiIndex = Dune::PDELab::MultiIndex<size_type,ContainerDepth>;
    using SizePrefix = Dune::PDELab::MultiIndex<size_type,ContainerDepth>;
    using PreBasis = PB;
    using LocalConstraints = decltype(std::declval<RootConstraintsContainer>().localView(std::declval<LocalViewTree>(),SubBasisPath{}));

    Basis(const EntitySet& entity_set, const PreBasis& pre_basis)
      : _pre_basis{pre_basis}
      , _entity_set{entity_set}
      , _ordering{Impl::makeOrdering(_pre_basis)}
      , _conforming_local_index_set{ std::make_shared<bool>(false) }
      , _sub_basis_path{TypeTree::treePath()}
    {
      update(_entity_set);
    }

    template<Dune::Concept::GridView _EntitySet, Concept::Impl::PreBasisTree _PreBasis, Concept::FixedSizeMultiIndex _SubBasisPath>
    friend class Basis;

    template<class... Args>
    Basis(const Basis<Args...>& other_basis, const EntitySet& entity_set, SubBasisPath sub_basis_path)
      : _pre_basis{other_basis._pre_basis}
      , _entity_set{entity_set}
      , _ordering{other_basis._ordering}
      , _constraints_container{other_basis._constraints_container}
      , _conforming_local_index_set{other_basis._conforming_local_index_set}
      , _sub_basis_path{sub_basis_path}
    {}

  private:
    template<class LocalTree>
    struct LocalIndexSetBase {
      using size_type = std::size_t;
      using MultiIndex = Dune::PDELab::MultiIndex<size_type,ContainerDepth>;
      using GlobalBasis = Basis;
      using Tree = LocalTree;

      LocalIndexSetBase(const Basis& basis, std::unique_ptr<LocalTree> ltree_storage)
        : _basis{basis}
        , _ltree_storage{move(ltree_storage)}
        , _ltree_view{_ltree_storage.get()}
      {
        if constexpr (DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE)
          _indices = std::make_unique<MultiIndex[]>(maxSize());
      }

      LocalIndexSetBase(const LocalIndexSetBase&) = delete;
      LocalIndexSetBase(LocalIndexSetBase&&) = default;

      LocalIndexSetBase& operator=(const LocalIndexSetBase&) = delete;
      LocalIndexSetBase& operator=(LocalIndexSetBase&&) = default;

      void bind(const Dune::Concept::Entity auto& element) {
        _ltree_storage->bind(element);
        assert(_size == 0);
        TypeTree::forEachLeafNode(*_ltree_storage, [&](auto& leaf, auto& path){
          leaf.setTreeOffset(_size);
          if constexpr (DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE)
            for(size_type dof = 0; dof != leaf.size(); ++dof)
              _indices[_size + dof] = leaf.index(dof);
          _size += leaf.size();
        });
        assert(_size <= maxSize());
        if constexpr (DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE)
          _indices_view = _indices.get();
        _ltree_view = _ltree_storage.get();
        // _mem_region = _basis.entitySetPartition().memoryRegion(element);
      }

      void unbind() {
        _ltree_storage->unbind();
        _ltree_view = _ltree_storage.get();
        _indices_view = nullptr;
        _size = 0;
      }

      friend void bind(const Dune::Concept::Entity auto& element, LocalIndexSetBase& lbasis0, auto& lbasis1) {
        lbasis0.bind(element);
        lbasis1.bind(element);
      }

      friend void unbind(LocalIndexSetBase& lbasis0, auto& lbasis1) {
        lbasis1.unbind();
        lbasis0.unbind();
      }

#if DUNE_PDELAB_ENABLE_DOUBLE_BIND
      friend void bind(const Dune::Concept::Entity auto& element, LocalIndexSetBase& lbasis0, LocalIndexSetBase& lbasis1) {
        lbasis0.doubleBind(element, lbasis1);
      }

      friend void unbind(LocalIndexSetBase& lbasis0, LocalIndexSetBase& lbasis1) {
        lbasis0.doubleUnbind(lbasis1);
      }
#endif

      [[nodiscard]] size_type size() const noexcept {
        return _size;
      }

      [[nodiscard]] const GlobalBasis& globalBasis() const noexcept {
        return _basis;
      }

      // Maximum number of coefficients that may be associated to a local view
      [[nodiscard]] size_type maxSize() const noexcept {
        return _basis._ordering->maxLocalCount();
      }

      // Whether local index sets match in all processors
      [[nodiscard]] auto conforming() const noexcept {
        return _basis.isLocalIndexSetConforming();
      }

      [[nodiscard]] MultiIndex index(size_type dof) const noexcept {
        if constexpr (not DUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE)
          DUNE_THROW(NotImplemented, "To enable this feature complie with `-DDUNE_PDELAB_ENABLE_INDEX_ON_ROOT_NODE`");
        else
          return _indices_view[dof];
      }

      [[nodiscard]] const Tree& tree() const noexcept {
        assert(_ltree_view);
        return *_ltree_view;
      }

      // [[nodiscard]] std::convertible_to<MemoryRegion> auto memoryRegion() const noexcept {
      //   return _mem_region; //TODO!!
      // }

      void lock() noexcept {
        assert(_ltree_view);
        _ltree_view->lock();
      }

      [[nodiscard]] bool try_lock() noexcept {
        assert(_ltree_view);
        return _ltree_view->try_lock();
      }

      void unlock() noexcept {
        assert(_ltree_view);
        _ltree_view->unlock();
      }

    protected:

      void doubleBind(const Dune::Concept::Entity auto& element, LocalIndexSetBase& other) {
        bind(element);
        if (_basis == other._basis) {
          // notice that we only share a (bound) view of our local tree
          other._ltree_view = _ltree_view;
          other._indices_view = _indices_view;
          other._size = _size;
          // other._mem_region = _mem_region;
        } else
          other.bind(element);
      }

      void doubleUnbind(LocalIndexSetBase& other) {
        if (_ltree_view == other._ltree_view) {
          other._ltree_view = nullptr;
          other._indices_view = nullptr;
          other._size = 0;
        } else
          other.unbind();
        unbind();
      }

      Basis _basis;
      size_type _size = 0;
      std::unique_ptr<LocalTree> _ltree_storage;
      LocalTree * _ltree_view = nullptr;
      std::unique_ptr<MultiIndex[]> _indices;
      MultiIndex const * _indices_view = nullptr;
      // MemoryRegion _mem_region;
    };

  public:

    class LocalIndexSet : public LocalIndexSetBase<LocalIndexSetTree> {
      using Base = LocalIndexSetBase<LocalIndexSetTree>;
    public:
      LocalIndexSet(const Basis& basis)
        : Base{basis, basis._ordering->makeLocalIndexSet(basis._ordering, TypeTree::treePath(), basis._sub_basis_path)}
      {}

      LocalIndexSet(const LocalIndexSet& other)
        : LocalIndexSet{other.globalBasis()}
      {}

      LocalIndexSet& bind(const Dune::Concept::Entity auto& element) {
        Base::bind(element);
        return *this;
      }

      LocalIndexSet& unbind() {
        Base::unbind();
        return *this;
      }
    };

    class LocalView : public LocalIndexSetBase<LocalViewTree> {
      using Base = LocalIndexSetBase<LocalViewTree>;
    public:
      using Element = typename EntitySet::template Codim<0>::Entity;

      LocalView(const Basis& basis)
        : Base{basis, basis._ordering->makeLocalView(basis._ordering, TypeTree::treePath(), basis._sub_basis_path)}
      {}

      LocalView(const LocalView& other) {
        (*this) = other;
      }

      LocalView& operator=(const LocalView& other) {
        return (*this) = LocalView{other.globalBasis()};
      }

      LocalView(LocalView&&) = default;
      LocalView& operator=(LocalView&&) = default;

      template<std::convertible_to<Element> E>
      LocalView& bind(E&& element) {
        bindElement(std::forward<E>(element));
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

      Element const * _entity_view = nullptr;
      std::optional<const Element> _entity_storage;
    };


    [[nodiscard]] size_type size(const SizePrefix& prefix) const noexcept {
      return _ordering->containerSize(reverse(prefix));
    }

    [[nodiscard]] auto dofLockHandle(const SizePrefix& multiindex) noexcept {
      return _ordering->dofLockHandle(reverse(multiindex));
    }

    [[nodiscard]] size_type size() const noexcept {
      return dimension();
    }

    [[nodiscard]] size_type dimension() const noexcept {
      return PDELab::containerEntry(*_ordering, _sub_basis_path).dimension();
    }

//     [[nodiscard]] const EntitySetPartition& entitySetPartition() const noexcept {
//       return _entity_set;
//     }

    [[nodiscard]] EntitySet entitySet() const noexcept {
      return _entity_set;
    }


    [[nodiscard]] std::string name() const {
      return PDELab::containerEntry(_pre_basis, _sub_basis_path).name();
    }

    void name(std::string new_name) {
      return PDELab::containerEntry(_pre_basis, _sub_basis_path).name(new_name);
    }


    [[nodiscard]] auto degree() const {
      return PDELab::containerEntry(_pre_basis, _sub_basis_path).degree();
    }


    template<class = void>
    void update(const EntitySet& entity_set) {
      static_assert(SubBasisPath::size() == 0, "Update of function basiss only be called on root basis");
      _entity_set = entity_set;
      _ordering->update(); // TODO: how to inform entity ordering about updated multidomaingrid?

      _constraints_container = makeConstraintsContainer(_pre_basis, ConstraintsContainerGenerator{});
      auto constraints_ops = TypeTree::makeTreeContainer(_pre_basis, [](auto& pre_basis_node){
        return pre_basis_node.constraintsOperator();
      });
      _constraints_container->assembleConstraints(*this, constraints_ops);

      // we need to know if other processors also have the same size information at the border
      // fixed-size orderings will always have the same local index set on each (sub)-entity
      // this also holds if there are more geometry types at the each codimension
      if (_ordering->fixedSize())
        return;

      // in sequential case we know the answer...
      if (entitySet().comm().size() == 0) {
        *_conforming_local_index_set = true;
        return;
      }

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


    [[nodiscard]] LocalView localView() const {
      return LocalView{*this};
    }

    [[nodiscard]] LocalConstraints localConstraints() const {
      return _constraints_container->localView(localView().tree(), _sub_basis_path);
    }

    [[nodiscard]] LocalIndexSet localIndexSet() const {
      return LocalIndexSet{*this};
    }

    // Whether local index sets match in all processors (constraints may still differ!sett)
    [[nodiscard]] auto isLocalIndexSetConforming() const noexcept {
      return *_conforming_local_index_set;
    }

    [[nodiscard]] auto fixedSize(std::size_t dim, std::size_t codim) const noexcept{
      return _ordering->fixedSize();
    }

    [[nodiscard]] bool contains(std::size_t dim, std::size_t codim) const noexcept {
      assert(dim == EntitySet::dimension); // not sure what is dim here
      return PDELab::containerEntry(*_ordering, _sub_basis_path).containsCodim(codim);
    }

    template<class Backend>
    [[nodiscard]] Container<Backend> makeContainer(Backend) const {
      Container<Backend> container{};
      containerResize(container, *this);
      forEachContainerEntry(container, []<class T>(T& v){v = T{0};});
      return container;
    }

    [[nodiscard]] friend bool operator==(const Basis& lhs, const Basis& rhs) noexcept {
      // // entity sets may be equally comparable
      // if constexpr (std::equality_comparable<EntitySetPartition>) {
      //   if (lhs._entity_set != rhs._entity_set)
      //     return false;
      // }

      if (lhs._ordering != rhs._ordering)
        return false;

      bool same_leafs = true;
      forEachLeafNode(*lhs._ordering, [&](const auto& lhs_ordering, auto path){
        same_leafs &= (&lhs_ordering == &PDELab::containerEntry(*rhs._ordering, path));
      });
      return same_leafs;
    }

    [[nodiscard]] friend bool operator!=(const Basis& lhs, const Basis& rhs) {
      return !(lhs == rhs);
    }

    // template<Concept::GridViewPartition OtherEntitySetPartition>
    // [[nodiscard]] auto subBasis(const OtherEntitySetPartition& partition, Concept::FixedSizeMultiIndex auto sub_basis_path) const {
    //   auto joined_sub_basis_path = join(this->_sub_basis_path, sub_basis_path);
    //   using JoinedSubBasisPath = decltype(joined_sub_basis_path);
    //   static_assert(requires {
    //     { PDELab::containerEntry(_pre_basis, joined_sub_basis_path).name() } -> std::convertible_to<std::string>; },
    //     "Child Pre-Basis node for the requested sub basis does not exist");
    //   return Basis<PreBasis, OtherEntitySetPartition, JoinedSubBasisPath>{*this, partition, joined_sub_basis_path};
    // }

    // [[nodiscard]] auto subBasis(Concept::FixedSizeMultiIndex auto sub_basis_path) const {
    //   return this->subBasis(this->entitySetPartition(), sub_basis_path);
    // }

  private:
    auto rootBasis() const
    {
      return Basis<EntitySet, PreBasis>{*this, this->entitySetPartition(), TypeTree::treePath()};
    }

    PreBasis _pre_basis;
    EntitySet _entity_set;
    std::shared_ptr<Ordering> _ordering;
    std::shared_ptr<RootConstraintsContainer> _constraints_container;
    std::shared_ptr<bool> _conforming_local_index_set;
    SubBasisPath _sub_basis_path;
  };


  // template<Concept::Impl::PreBasisTree PreBasis, Concept::GridViewPartition EntitySetPartition>
  // [[nodiscard]] Concept::Basis auto makeOrderedBasis(const PreBasis& pre_basis, const EntitySetPartition& partition)
  // {
  //   return Basis{pre_basis, partition};
  // }

  template<Concept::Impl::PreBasisTree PreBasis, Dune::Concept::GridView EntitySet>
  [[nodiscard]] /*Concept::Basis*/ auto makeBasis(const EntitySet& entity_set, const PreBasis& pre_basis)
  {
    return Basis{entity_set, pre_basis};
  }

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_BASIS_BASIS_HH
