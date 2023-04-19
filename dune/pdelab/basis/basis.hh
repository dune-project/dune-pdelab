#ifndef DUNE_PDELAB_BASIS_BASIS_HH
#define DUNE_PDELAB_BASIS_BASIS_HH

#include <dune/pdelab/concepts/multiindex.hh>
// #include <dune/assembler/concepts/space.hh>
// #include <dune/assembler/concepts/local_space.hh>

// #include <dune/assembler/common/partition/identity.hh>
// #include <dune/assembler/common/reservedmultiindex.hh>
// #include <dune/assembler/common/container_entry.hh>
// #include <dune/assembler/common/container_resize.hh>
// #include <dune/assembler/common/communication/entity_data_handler.hh>

// #include <dune/assembler/space/concept.hh>
// #include <dune/assembler/space/ordering.hh>

// #include <dune/assembler/space/constraints/container.hh>

// #include <dune/assembler/space/leaf.hh>
// #include <dune/assembler/space/composite.hh>

// #include <dune/grid/common/mcmgmapper.hh>

// #ifndef DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE
// #define DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE 0
// #endif

// #ifndef DUNE_ASSEMBLER_ENABLE_DOUBLE_BIND
// #define DUNE_ASSEMBLER_ENABLE_DOUBLE_BIND 1
// #endif

// namespace Dune::Assembler {

//   template<Concept::Impl::SpaceTree PS, Concept::EntitySetPartition ESP, Concept::MultiIndex SubSpacePath = MultiIndex<> >
//   class Space {
//     // We need some types sub-space independent so that we can shared them between templated sub-spaces
//     using Ordering = std::decay_t<decltype(*Impl::makeOrdering(std::declval<PS>()))>;
//     using LocalViewTree = std::decay_t<decltype(*std::declval<Ordering>().makeLocalView(std::shared_ptr<Ordering>{}, multiIndex(), SubSpacePath{}))>;
//     using LocalIndexSetTree = std::decay_t<decltype(*std::declval<Ordering>().makeLocalIndexSet(std::shared_ptr<Ordering>{}, multiIndex(), SubSpacePath{}))>;
//     static constexpr std::size_t ContainerDepth = Ordering::maxContainerDepth();
//     static_assert(ContainerDepth > 0);

//     using RootLocalViewTree = std::decay_t<decltype(*std::declval<Ordering>().makeLocalView(std::shared_ptr<Ordering>{}, multiIndex(), multiIndex()))>;

//     struct ConstraintsContainerGenerator {
//       template<Concept::Impl::SpaceLeaf PreSpaceLeaf, Concept::MultiIndex Path>
//       auto operator()(const PreSpaceLeaf& pre_space_leaf, Path) const {
//         using EntitySet = PreSpaceLeaf::Traits::MergingStrategy::EntitySet;
//         using MultiIndex = typename TypeTree::ChildForTreePath<RootLocalViewTree,Path>::MultiIndex;
//         using ConstraintsContainer = PreSpaceLeaf::Traits::ConstraintsOperator::template Container<MultiIndex, EntitySet>;
//         return std::make_shared<ConstraintsContainer>(pre_space_leaf.mergingStrategy().gridView());
//       }
//     };

//     using RootConstraintsContainer = std::decay_t<decltype(*makeConstraintsContainer(std::declval<const PS&>(), ConstraintsContainerGenerator{}))>;

//   public:

//     template<class Backend>
//     using Container = std::decay_t<decltype(Ordering::template makeVectorContainer<Backend>())>;

//     using size_type = std::size_t;
//     using EntitySetPartition = ESP;
//     using GridView = typename EntitySetPartition::GridView;
//     using MultiIndex = Dune::Assembler::ReservedMultiIndex<size_type,ContainerDepth>;
//     using SizePrefix = Dune::Assembler::ReservedMultiIndex<size_type,ContainerDepth>;
//     using PreSpace = PS;
//     using LocalConstraints = decltype(std::declval<RootConstraintsContainer>().localView(std::declval<LocalViewTree>(),SubSpacePath{}));

//     Space(const PreSpace& pre_space, const EntitySetPartition& partition)
//       : _pre_space{pre_space}
//       , _partition{partition}
//       , _ordering{Impl::makeOrdering(_pre_space)}
//       , _conforming_local_index_set{ std::make_shared<bool>(false) }
//       , _sub_space_path{multiIndex()}
//     {
//       update(_partition);
//     }

//     template<Concept::Impl::SpaceTree _PreSpace, Concept::EntitySetPartition _EntitySet, Concept::MultiIndex _SubSpacePath>
//     friend class Space;

//     template<class... Args>
//     Space(const Space<Args...>& other_space, const EntitySetPartition& partition, SubSpacePath sub_space_path)
//       : _pre_space{other_space._pre_space}
//       , _partition{partition}
//       , _ordering{other_space._ordering}
//       , _constraints_container{other_space._constraints_container}
//       , _conforming_local_index_set{other_space._conforming_local_index_set}
//       , _sub_space_path{sub_space_path}
//     {}

//   private:
//     template<class LocalTree>
//     struct LocalIndexSetBase {
//       using size_type = std::size_t;
//       using MultiIndex = Dune::Assembler::ReservedMultiIndex<size_type,ContainerDepth>;
//       using GlobalBasis = Space;
//       using Tree = LocalTree;

//       LocalIndexSetBase(const Space& space, std::unique_ptr<LocalTree> ltree_storage)
//         : _space{space}
//         , _ltree_storage{move(ltree_storage)}
//         , _ltree_view{_ltree_storage.get()}
//       {
//         if constexpr (DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE)
//           _indices = std::make_unique<MultiIndex[]>(maxSize());
//       }

//       LocalIndexSetBase(const LocalIndexSetBase&) = delete;

//       LocalIndexSetBase(LocalIndexSetBase&& other)
//         : _space{move(other._space)}
//         , _ltree_storage{move(other._ltree_storage)}
//         , _ltree_view{other._ltree_view}
//         , _indices_view{other._indices_view}
//       {
//         if constexpr (DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE)
//           _indices = move(other._indices);
//       }

//       void bind(const Dune::Concept::Entity auto& element) {
//         _ltree_storage->bind(element);
//         assert(_size == 0);
//         TypeTree::forEachLeafNode(*_ltree_storage, [&](auto& leaf, auto& path){
//           leaf.setTreeOffset(_size);
//           if constexpr (DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE)
//             for(size_type dof = 0; dof != leaf.size(); ++dof)
//               _indices[_size + dof] = leaf.index(dof);
//           _size += leaf.size();
//         });
//         assert(_size <= maxSize());
//         if constexpr (DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE)
//           _indices_view = _indices.get();
//         _ltree_view = _ltree_storage.get();
//         _mem_region = _space.entitySetPartition().memoryRegion(element);
//       }

//       void unbind() {
//         _ltree_storage->unbind();
//         _ltree_view = _ltree_storage.get();
//         _indices_view = nullptr;
//         _size = 0;
//       }

//       friend void bind(const Dune::Concept::Entity auto& element, LocalIndexSetBase& lspace0, auto& lspace1) {
//         lspace0.bind(element);
//         lspace1.bind(element);
//       }

//       friend void unbind(LocalIndexSetBase& lspace0, auto& lspace1) {
//         lspace1.unbind();
//         lspace0.unbind();
//       }

// #if DUNE_ASSEMBLER_ENABLE_DOUBLE_BIND
//       friend void bind(const Dune::Concept::Entity auto& element, LocalIndexSetBase& lspace0, LocalIndexSetBase& lspace1) {
//         lspace0.doubleBind(element, lspace1);
//       }

//       friend void unbind(LocalIndexSetBase& lspace0, LocalIndexSetBase& lspace1) {
//         lspace0.doubleUnbind(lspace1);
//       }
// #endif

//       [[nodiscard]] size_type size() const noexcept {
//         return _size;
//       }

//       [[nodiscard]] const GlobalBasis& globalBasis() const noexcept {
//         return _space;
//       }

//       // Maximum number of coefficients that may be associated to a local view
//       [[nodiscard]] size_type maxSize() const noexcept {
//         return _space._ordering->maxLocalCount();
//       }

//       // Whether local index sets match in all processors
//       [[nodiscard]] auto conforming() const noexcept {
//         return _space.isLocalIndexSetConforming();
//       }

//       [[nodiscard]] MultiIndex index(size_type dof) const noexcept {
//         if constexpr (not DUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE)
//           DUNE_THROW(NotImplemented, "To enable this feature complie with `-DDUNE_ASSEMBLER_ENABLE_INDEX_ON_ROOT_NODE`");
//         else
//           return _indices_view[dof];
//       }

//       [[nodiscard]] const Tree& tree() const noexcept {
//         assert(_ltree_view);
//         return *_ltree_view;
//       }

//       [[nodiscard]] std::convertible_to<MemoryRegion> auto memoryRegion() const noexcept {
//         return _mem_region; //TODO!!
//       }

//       void lock() noexcept {
//         assert(_ltree_view);
//         _ltree_view->lock();
//       }

//       [[nodiscard]] bool try_lock() noexcept {
//         assert(_ltree_view);
//         return _ltree_view->try_lock();
//       }

//       void unlock() noexcept {
//         assert(_ltree_view);
//         _ltree_view->unlock();
//       }

//     protected:

//       void doubleBind(const Dune::Concept::Entity auto& element, LocalIndexSetBase& other) {
//         bind(element);
//         if (_space == other._space) {
//           // notice that we only share a (bound) view of our local tree
//           other._ltree_view = _ltree_view;
//           other._indices_view = _indices_view;
//           other._size = _size;
//           other._mem_region = _mem_region;
//         } else
//           other.bind(element);
//       }

//       void doubleUnbind(LocalIndexSetBase& other) {
//         unbind();
//         if (_space == other._space) {
//           other._ltree_view = nullptr;
//           other._indices_view = nullptr;
//           other._size = 0;
//         } else
//           other.unbind();
//       }

//       Space _space;
//       size_type _size = 0;
//       std::unique_ptr<LocalTree> _ltree_storage;
//       LocalTree * _ltree_view = nullptr;
//       std::unique_ptr<MultiIndex[]> _indices;
//       MultiIndex const * _indices_view = nullptr;
//       MemoryRegion _mem_region;
//     };

//   public:

//     class LocalIndexSet : public LocalIndexSetBase<LocalIndexSetTree> {
//       using Base = LocalIndexSetBase<LocalIndexSetTree>;
//     public:
//       LocalIndexSet(const Space& space)
//         : Base{space, space._ordering->makeLocalIndexSet(space._ordering, multiIndex(), space._sub_space_path)}
//       {}

//       LocalIndexSet(const LocalIndexSet& other)
//         : LocalIndexSet{other.globalBasis()}
//       {}

//       LocalIndexSet& bind(const Dune::Concept::Entity auto& element) {
//         Base::bind(element);
//         return *this;
//       }

//       LocalIndexSet& unbind() {
//         Base::unbind();
//         return *this;
//       }
//     };

//     class LocalView : public LocalIndexSetBase<LocalViewTree> {
//       using Base = LocalIndexSetBase<LocalViewTree>;
//     public:
//       using Element = typename EntitySetPartition::Entity;

//       LocalView(const Space& space)
//         : Base{space, space._ordering->makeLocalView(space._ordering, multiIndex(), space._sub_space_path)}
//       {}

//       LocalView(const LocalView& other)
//         : LocalView{other.globalBasis()}
//       {}

//       template<std::convertible_to<Element> E>
//       LocalView& bind(E&& element) {
//         bindElement(std::forward<E>(element));
//         Base::bind(this->element());
//         return *this;
//       }

//       LocalView& unbind() {
//         Base::unbind();
//         _entity_view = nullptr;
//         _entity_storage = std::nullopt;
//         return *this;
//       }

//       template<std::convertible_to<Element> E>
//       friend void bind(E&& element, LocalView& lspace0, auto& lspace1) {
//         lspace0.bind(std::forward<E>(element));
//         lspace1.bind(lspace0.element());
//       }

//       friend void unbind(LocalView& lspace0, auto& lspace1) {
//         lspace1.unbind();
//         lspace0.unbind();
//       }

// #if DUNE_ASSEMBLER_ENABLE_DOUBLE_BIND
//       template<std::convertible_to<Element> E>
//       friend void bind(E&& element, LocalView& lspace0, LocalView& lspace1) {
//         lspace0.bindElement(std::forward<E>(element));
//         lspace1.bindElement(lspace0.element());
//         lspace0.doubleBind(lspace0.element(), lspace1);
//       }

//       friend void unbind(LocalView& lspace0, LocalView& lspace1) {
//         lspace0.doubleUnbind(lspace1);
//         lspace0._entity_view = lspace1._entity_view = nullptr;
//         lspace0._entity_storage = std::nullopt;
//         lspace1._entity_storage = std::nullopt;
//       }
// #endif

//       [[nodiscard]] const Element& element() const noexcept {
//         assert(_entity_view);
//         return *_entity_view;
//       };

//     private:
//       void bindElement(Element&& element) {
//         // the caller assigned the ownership to us
//         _entity_storage.emplace(std::move(element));
//         _entity_view = &(*_entity_storage);
//       }

//       void bindElement(const Element& element) {
//           // ownership is managed elsewhere
//         _entity_view = &element;
//       }

//       Element const * _entity_view = nullptr;
//       std::optional<const Element> _entity_storage;
//     };


//     [[nodiscard]] size_type size(const SizePrefix& prefix) const noexcept {
//       return _ordering->containerSize(reverse(prefix));
//     }

//     [[nodiscard]] auto dofLockHandle(const SizePrefix& multiindex) noexcept {
//       return _ordering->dofLockHandle(reverse(multiindex));
//     }

//     [[nodiscard]] size_type size() const noexcept {
//       return dimension();
//     }

//     [[nodiscard]] size_type dimension() const noexcept {
//       return TypeTree::child(*_ordering, _sub_space_path).dimension();
//     }

//     [[nodiscard]] const EntitySetPartition& entitySetPartition() const noexcept {
//       return _partition;
//     }

//     [[nodiscard]] GridView gridView() const noexcept {
//       return _partition.gridView();
//     }

//     [[nodiscard]] std::string name() const {
//       return TypeTree::child(_pre_space, _sub_space_path).name();
//     }

//     void name(std::string new_name) {
//       return TypeTree::child(_pre_space, _sub_space_path).name(new_name);
//     }


//     [[nodiscard]] auto degree() const {
//       return TypeTree::child(_pre_space, _sub_space_path).degree();
//     }


//     template<class = void>
//     void update(const EntitySetPartition& partition) {
//       static_assert(SubSpacePath::size() == 0, "Update of function spaces only be called on root space");
//       _partition = partition;
//       _ordering->update(); // TODO: how to inform entity ordering about updated multidomaingrid?

//       _constraints_container = makeConstraintsContainer(_pre_space, ConstraintsContainerGenerator{});
//       auto constraints_ops = TypeTree::makeTreeContainer(_pre_space, [](auto& pre_space_node){
//         return pre_space_node.constraintsOperator();
//       });
//       _constraints_container->assembleConstraints(*this, constraints_ops);

//       // we need to know if other processors also have the same size information at the border
//       // fixed-size orderings will always have the same local index set on each (sub)-entity
//       // this also holds if there are more geometry types at the each codimension
//       if (_ordering->fixedSize())
//         return;

//       // in sequential case we know the answer...
//       if (gridView().comm().size() == 0) {
//         *_conforming_local_index_set = true;
//         return;
//       }

//       MultipleCodimMultipleGeomTypeMapper mapper{gridView(), [](auto, auto){return 1;}};

//       // technically we don't need the container, but it's easier to reuse the entity data handle
//       std::vector<std::size_t> mismatching_sizes(mapper.size(), 0);
//       // notice that we always calculate this with the root node.
//       // this allows us to share and reuse the result to any sub-space
//       auto local_index_set = localIndexSet();

//       // communicate local index set sizes and compare them at receiving end
//       applyToDistributedEntities(
//         mapper,
//         mismatching_sizes,
//         Dune::InterfaceType::All_All_Interface,
//         [&](auto phase,
//             const auto& entity,
//             auto& remote_value,
//             auto& local_value) {
//           local_index_set.bind(entity);
//           if constexpr (decltype(phase)::value == Communication::gather)
//             remote_value = local_index_set.size();
//           else
//             local_value = (local_index_set.size() != remote_value);
//         });

//       // accumulate number of entities that mismatch
//       auto missmatching =
//         std::accumulate(begin(mismatching_sizes), end(mismatching_sizes), 0);
//       *_conforming_local_index_set = (missmatching == 0);
//     }


//     [[nodiscard]] LocalView localView() const {
//       return LocalView{*this};
//     }

//     [[nodiscard]] LocalConstraints localConstraints() const {
//       return _constraints_container->localView(localView().tree(), _sub_space_path);
//     }

//     [[nodiscard]] LocalIndexSet localIndexSet() const {
//       return LocalIndexSet{*this};
//     }

//     // Whether local index sets match in all processors (constraints may still differ!sett)
//     [[nodiscard]] auto isLocalIndexSetConforming() const noexcept {
//       return *_conforming_local_index_set;
//     }

//     [[nodiscard]] auto fixedSize(std::size_t dim, std::size_t codim) const noexcept{
//       return _ordering->fixedSize();
//     }

//     [[nodiscard]] bool contains(std::size_t dim, std::size_t codim) const noexcept {
//       assert(dim == GridView::dimension); // not sure what is dim here
//       return TypeTree::child(*_ordering, _sub_space_path).containsCodim(codim);
//     }

//     template<class Backend>
//     [[nodiscard]] Container<Backend> makeContainer(Backend) const {
//       Container<Backend> container{};
//       containerResize(container, *this);
//       forEachContainerEntry(container, []<class T>(T& v){v = T{0};});
//       return container;
//     }

//     [[nodiscard]] friend bool operator==(const Space& lhs, const Space& rhs) noexcept {
//       // entity sets may be equally comparable
//       if constexpr (std::equality_comparable<EntitySetPartition>) {
//         if (lhs._partition != rhs._partition)
//           return false;
//       }

//       if (lhs._ordering != rhs._ordering)
//         return false;

//       bool same_leafs = true;
//       forEachLeafNode(*lhs._ordering, [&](const auto& lhs_ordering, auto path){
//         same_leafs &= (&lhs_ordering == &TypeTree::child(*rhs._ordering, path));
//       });
//       return same_leafs;
//     }

//     [[nodiscard]] friend bool operator!=(const Space& lhs, const Space& rhs) {
//       return !(lhs == rhs);
//     }

//     template<Concept::EntitySetPartition OtherEntitySetPartition>
//     [[nodiscard]] auto subSpace(const OtherEntitySetPartition& partition, Concept::MultiIndex auto sub_space_path) const {
//       auto joined_sub_space_path = join(this->_sub_space_path, sub_space_path);
//       using JoinedSubSpacePath = decltype(joined_sub_space_path);
//       static_assert(requires {
//         { TypeTree::child(_pre_space, joined_sub_space_path).name() } -> std::convertible_to<std::string>; },
//         "Child Pre-Space node for the requested sub space does not exist");
//       return Space<PreSpace, OtherEntitySetPartition, JoinedSubSpacePath>{*this, partition, joined_sub_space_path};
//     }

//     [[nodiscard]] auto subSpace(Concept::MultiIndex auto sub_space_path) const {
//       return this->subSpace(this->entitySetPartition(), sub_space_path);
//     }

//   private:
//     auto rootBasis() const
//     {
//       return Space<PreSpace, EntitySetPartition>{*this, this->entitySetPartition(), multiIndex()};
//     }

//     PreSpace _pre_space;
//     EntitySetPartition _partition;
//     std::shared_ptr<Ordering> _ordering;
//     std::shared_ptr<RootConstraintsContainer> _constraints_container;
//     std::shared_ptr<bool> _conforming_local_index_set;
//     SubSpacePath _sub_space_path;
//   };


//   template<Concept::Impl::SpaceTree PreSpace, Concept::EntitySetPartition EntitySetPartition>
//   [[nodiscard]] Concept::Space auto makeOrderedSpace(const PreSpace& pre_space, const EntitySetPartition& partition)
//   {
//     return Space{pre_space, partition};
//   }

//   template<Concept::Impl::SpaceTree PreSpace, Dune::Concept::GridView GridView>
//   [[nodiscard]] Concept::Space auto makeOrderedSpace(const PreSpace& pre_space, const GridView& grid_view)
//   {
//     return Space{pre_space, EntitySetPartition::Identity{grid_view}};
//   }

// } // namespace Dune::Assembler

#endif // DUNE_PDELAB_BASIS_BASIS_HH
