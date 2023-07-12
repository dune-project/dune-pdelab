#ifndef DUNE_PDELAB_BASIS_ORDERING_ENTITY_SET_HH
#define DUNE_PDELAB_BASIS_ORDERING_ENTITY_SET_HH

#include <dune/pdelab/basis/prebasis/concept.hh>

// #include <dune/pdelab/common/concurrency/vector_adaptivelock.hh>
// #include <dune/pdelab/common/concurrency/vector_spinlock.hh>
#include <dune/pdelab/common/concurrency/vector_bitspinlock.hh>
#include <dune/pdelab/common/entity_cast.hh>

#include <dune/pdelab/concepts/treenode.hh>
#include <dune/pdelab/concepts/lockable.hh>
#include <dune/pdelab/concepts/multiindex.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/typetree/traversal.hh>
#include <dune/typetree/treecontainer.hh>

#include <dune/grid/concepts/gridview.hh>

#include <memory>
#include <optional>
#include <thread>
#include <mutex>

#ifndef DUNE_PDELAB_VECTOR_LOCK_TYPE
#define DUNE_PDELAB_VECTOR_LOCK_TYPE VectorBitSpinLock // VectorSpinLock, VectorBitSpinLock, VectorAdaptiveLock
#endif

namespace Dune::PDELab::inline Experimental::Impl {

    using VectorLockType = DUNE_PDELAB_VECTOR_LOCK_TYPE;

  template<Concept::TreeNode EntityOrdering, class MS>
  class EntitySetOrdering : public EntityOrdering {
    using SizeType = typename MS::SizeType;
    using EntitySet = typename MS::EntitySet;
    using MergingStrategy = MS;
  public:

    EntitySetOrdering(EntityOrdering&& entity_ordering)
      : EntityOrdering{std::move(entity_ordering)}
    {}

    EntitySetOrdering(const EntitySetOrdering&) = delete;
    EntitySetOrdering(EntitySetOrdering&&) = delete;

    EntitySetOrdering& operator=(const EntitySetOrdering&) = delete;
    EntitySetOrdering& operator=(EntitySetOrdering&&) = default;

    static constexpr auto containerBlocked() {
      return std::integral_constant<bool,MergingStrategy::Blocked>{};
    }

    static constexpr auto maxContainerDepth() {
      if constexpr (containerBlocked())
        return EntityOrdering::maxContainerDepth()+1;
      else
        return EntityOrdering::maxContainerDepth();
    }

    // slice this as base class
    const EntityOrdering& entityOrdering() const { return static_cast<const EntityOrdering&>(*this); }
    EntityOrdering& entityOrdering() { return static_cast<EntityOrdering&>(*this); }

    using EntityOrdering::entitySet;
    using EntityOrdering::fixedSize;
    using EntityOrdering::maxLocalCount;

    template<Concept::MultiIndex CompositionSuffix>
    Concept::MultiIndex auto firstContainerIndex(CompositionSuffix comp_suf, const SizeType& gt_index, const SizeType& e_index) const noexcept {
      // Note: Multi-index is read Outer->Inner
      // get local container suffix
      assert(this->containsGeometry(gt_index));
      const auto lcont_prefix = entityOrdering().firstContainerIndex(comp_suf, gt_index, e_index);

      // add blocking if necessary
      const auto bcont_prefix = [&](){
        if constexpr (containerBlocked())
          return push_front(lcont_prefix, Indices::_0);
        else
          return lcont_prefix;
      }();

      // accumulate offset for this entity
      return accumulate_front(bcont_prefix, blockOffset(gt_index, e_index));
    }

    template<Concept::MultiIndex ContainerSuffix>
    std::size_t containerSize(ContainerSuffix suffix) const noexcept
    {
      // Note: Multi-index is read Inner->Outer
      if (suffix.size() == 0)
        return blockCount();
      else {
        // the next index to find out its size
        auto back_index = suffix.back();

        // we first have to figure out the entity index
        const auto [gt_index, entity_index] = getGeometryEntityIndex(back_index);

        // then, the entity ordering knows the size for a given entity.
        if constexpr (containerBlocked())
          return entityOrdering().containerSize(pop_back(suffix), gt_index, entity_index);
        else {
          // remove offset introduced by the entity_index
          suffix.back() -= blockOffset(gt_index, entity_index);
          return entityOrdering().containerSize(suffix, gt_index, entity_index);
        }
      }
    }


    void update()
    {
      entityOrdering().update();

      if (fixedSize())
        updateFixedSizeOrdering();
      else
        updateVariableSizeOrdering();

      updateEntitySetMutex();

      // TODO: attach handle to entity set to inform us whenever it changes!
    }

    using EntityOrdering::blockCount;

    [[nodiscard]] SizeType blockCount() const noexcept {
      if (fixedSize())
        return _gt_dof_offsets->back();
      else
        return _entity_dof_offsets->back();
    }

    [[nodiscard]] Concept::Lockable auto entityLockHandle(std::size_t gt_index, std::size_t entity_index) noexcept {
      return _gt_mutex[gt_index][entity_index];
    }

    // gets a handle to lock an specific degree of freedom
    [[nodiscard]] Concept::Lockable auto dofLockHandle(Concept::MultiIndex auto suffix) {
      // Note: Multi-index is read Inner->Outer
      // we first have to figure out the entity index
      const auto [gt_index, entity_index] = getGeometryEntityIndex(suffix.back());
      return entityLockHandle(gt_index,entity_index);
    }

  private:

    auto getGeometryEntityIndex(SizeType i) const {
      // we just need to make the inverse computation of the mapIndex funtion to find the entity index
      SizeType gt_index, entity_index;
      auto dof_begin = begin(fixedSize() ? *_gt_dof_offsets : *_entity_dof_offsets);
      auto dof_end = end(fixedSize() ? *_gt_dof_offsets : *_entity_dof_offsets);
      auto dof_it = std::prev(std::upper_bound(dof_begin, dof_end, i));
      auto dof_dist = std::distance(dof_begin, dof_it);

      if (fixedSize()) {
        gt_index = dof_dist;
        // assert(i >= *dof_it);
        entity_index = i - *dof_it;
        entity_index /= entityOrdering().blockCount(gt_index);
      } else {
        auto gt_begin = begin(*_gt_entity_offsets);
        auto gt_end = end(*_gt_entity_offsets);
        auto gt_it = std::prev(std::upper_bound(gt_begin, gt_end, dof_dist));
        gt_index = std::distance(gt_begin, gt_it);
        // assert(dof_dist >= *gt_it);
        entity_index = dof_dist - *gt_it;
      }
      return std::make_pair(gt_index, entity_index);
    }

    SizeType blockOffset(SizeType geometry_index, SizeType entity_index) const {
      if (fixedSize()) {
        const auto& gt_offset = (*_gt_dof_offsets)[geometry_index];
        if constexpr (containerBlocked())
          return gt_offset + entity_index;
        else
          return gt_offset + entity_index * entityOrdering().blockCount(geometry_index);
      } else {
        const auto& gt_offset = (*_gt_entity_offsets)[geometry_index];
        return (*_entity_dof_offsets)[gt_offset + entity_index];
      }
    }

    void updateFixedSizeOrdering() {
      assert(fixedSize());

      if (not _gt_dof_offsets)
        _gt_dof_offsets = std::make_unique<std::vector<SizeType>>();
      const SizeType gt_count = GlobalGeometryTypeIndex::size(EntitySet::dimension);
      _gt_dof_offsets->assign(gt_count + 1,0);

      for (std::size_t codim = 0; codim <= EntitySet::dimension; ++codim) {
        for (const auto& gt : entitySet().indexSet().types(codim)) {
          const SizeType gt_index = GlobalGeometryTypeIndex::index(gt);
          if (not entityOrdering().containsGeometry(gt_index))
            continue;
          SizeType gt_size = entityOrdering().blockCount(gt_index);
          const SizeType gt_entity_count = entitySet().indexSet().size(gt);
          if (containerBlocked())
            gt_size = (gt_size > 0);
          (*_gt_dof_offsets)[gt_index + 1] = gt_size * gt_entity_count;
        }
      }

      std::partial_sum(begin(*_gt_dof_offsets),end(*_gt_dof_offsets),begin(*_gt_dof_offsets));
    }

    void updateVariableSizeOrdering() {
      if constexpr (true or containerBlocked()) { // TODO!!!
        // in this case, we need to count the number of _used_ sub blocks
        assert(not fixedSize());
        if (not _gt_entity_offsets)
          _gt_entity_offsets = std::make_unique<std::vector<SizeType>>();
        const SizeType gt_count = GlobalGeometryTypeIndex::size(EntitySet::dimension);
        _gt_entity_offsets->assign(gt_count + 1,0);

        for (std::size_t codim = 0; codim <= EntitySet::dimension; ++codim) {
          for (const auto& gt : entitySet().indexSet().types(codim)) {
            if (!entityOrdering().containsGeometry(gt))
              continue;
            const SizeType gt_index = GlobalGeometryTypeIndex::index(gt);
            (*_gt_entity_offsets)[gt_index + 1] = entitySet().indexSet().size(gt);
          }
        }

        std::partial_sum(begin(*_gt_entity_offsets),end(*_gt_entity_offsets),begin(*_gt_entity_offsets));
        if (not _entity_dof_offsets)
          _entity_dof_offsets = std::make_unique<std::vector<SizeType>>();
        _entity_dof_offsets->assign(_gt_entity_offsets->back()+1,0);

        SizeType carry_block = 0;
        SizeType index = 0;
        for (SizeType gt_index = 0; gt_index < gt_count; ++gt_index) {
          if (not entityOrdering().containsGeometry(gt_index))
            continue;
          const SizeType entity_count = (*_gt_entity_offsets)[gt_index + 1] - (*_gt_entity_offsets)[gt_index];
          for (SizeType entity_index = 0; entity_index < entity_count; ++entity_index) {
            const SizeType size = entityOrdering().blockCount(gt_index, entity_index);
            carry_block += (containerBlocked() ? (size > 0) : size);
            (*_entity_dof_offsets)[++index] = carry_block;
          }
        }
      } else {
        // in this case, we need to count the number of sub blocks, this is
        // exactly the same as what the entity ordering does so we reuse them
        _gt_entity_offsets = EntityOrdering::_gt_entity_offsets;
        _entity_dof_offsets = EntityOrdering::_entity_dof_offsets;
      }
    }

    //! @todo Use parent entity set in the loop
    void updateEntitySetMutex() {
      const SizeType gt_count = GlobalGeometryTypeIndex::size(EntitySet::dimension);
      _gt_mutex = std::vector<VectorLockType>(gt_count);
      for (std::size_t codim = 0; codim <= EntitySet::dimension; ++codim) {
        for (const auto& gt : entitySet().indexSet().types(codim)) {
          const SizeType gt_index = GlobalGeometryTypeIndex::index(gt);
          const SizeType gt_entity_count = entitySet().indexSet().size(gt);
          _gt_mutex[gt_index].resize(gt_entity_count);
        }
      }
    }

    template<Concept::TreeNode EntityLocalIndexSet, Concept::TreeNode GlobalOrdering, Concept::MultiIndex SubSpacePath>
    class LocalIndexSet : public EntityLocalIndexSet {

    public:

      LocalIndexSet(EntityLocalIndexSet&& elbasis, const std::shared_ptr<GlobalOrdering>& ordering, EntitySetOrdering& _entity_set_ordering, SubSpacePath sub_space_path)
        : EntityLocalIndexSet{std::move(elbasis)}
        , _ordering{ordering}
        , _entity_set_ordering{_entity_set_ordering}
        , _sub_space_path{sub_space_path}
      {}

      LocalIndexSet(const LocalIndexSet&) = delete;
      LocalIndexSet(LocalIndexSet&&) = default;

      LocalIndexSet& operator=(const LocalIndexSet&) = delete;
      LocalIndexSet& operator=(LocalIndexSet&&) = default;


      EntityLocalIndexSet& entityLocalIndexSet() {
        return static_cast<EntityLocalIndexSet&>(*this);
      }

      template<class AssemblyEntity>
      void bind(const AssemblyEntity& entity) noexcept {
        const auto& entity_set = _entity_set_ordering.entitySet();
        if (not entity_set.contains(entity)) {
          forEachLeafNode(entityLocalIndexSet(), [&](auto& leaf_local_index_set, auto) {
            leaf_local_index_set.indices().clear();
          });
          return;
        }

        constexpr auto codim = std::integral_constant<std::size_t, AssemblyEntity::codimension>{};
        if (not (EntitySetOrdering::mayContainCodim(codim) and _entity_set_ordering.containsCodim(codim))) {
          forEachLeafNode(entityLocalIndexSet(), [&](auto& leaf_local_index_set, auto) {
            leaf_local_index_set.indices().clear();
          });
          return;
        }

        _gt_index = GlobalGeometryTypeIndex::index(entity.type());
        _entity_index = entity_set.indexSet().index(entity);

        forEachLeafNode(entityLocalIndexSet(), [&]<class LeafLocalIndex>(LeafLocalIndex& leaf_local_index_set, const auto& suffix) {
          auto ordering_path = join(_sub_space_path, leaf_local_index_set.orderingViewPath());
          const auto& leaf_ordering = PDELab::containerEntry(*_ordering, ordering_path);

          auto& indices = leaf_local_index_set.indices();
          indices.resize(leaf_ordering.blockCount(_gt_index, _entity_index));
          if (indices.empty()) return;

          indices[0] = _ordering->firstContainerIndex(ordering_path, _gt_index, _entity_index);
          for (std::size_t dof = 1; dof < indices.size(); ++dof)
            indices[dof] = accumulate_back(indices[0], dof);
        });
      }

      void unbind() {}

      void lock() noexcept {
        mutex().lock();
      }

      bool try_lock() noexcept {
        return mutex().try_lock();
      }

      void unlock() noexcept {
        mutex().unlock();
      }

    private:
      auto mutex() {
        return _entity_set_ordering.entityLockHandle(_gt_index,_entity_index);
      }

      std::shared_ptr<GlobalOrdering> _ordering;
      EntitySetOrdering& _entity_set_ordering;
      SubSpacePath _sub_space_path;
      typename MS::SizeType _gt_index, _entity_index;
    };

    template<Concept::TreeNode EntityLocalView, Concept::TreeNode GlobalOrdering, Concept::MultiIndex SubSpacePath>
    class LocalView : public EntityLocalView {
    public:
      using Element = typename EntitySet::template Codim<0>::Entity;

      LocalView(EntityLocalView&& elbasis, const std::shared_ptr<GlobalOrdering>& ordering, EntitySetOrdering& entity_set_ordering, SubSpacePath sub_space_path)
        : EntityLocalView{std::move(elbasis)}
        , _index_cache(1) // DG-FV branch
        , _ordering{ordering}
        , _entity_set_ordering{entity_set_ordering}
        , _sub_space_path{sub_space_path}
      {}

      LocalView(const LocalView&) = delete;
      LocalView(LocalView&&) = default;

      LocalView& operator=(const LocalView&) = delete;
      LocalView& operator=(LocalView&&) = default;

      EntityLocalView& entityLocalView() {
        return static_cast<EntityLocalView&>(*this);
      }

      // notice that entity needs to be stored somewhere. If the outer scope
      // cannot ensure that the entity will be alive while using the view
      // it should move the entity here otherwise the contents will dangle
      template<class AssemblyEntity>
      void bind(AssemblyEntity&& entity) noexcept {

        const auto& entity_set = _entity_set_ordering.entitySet();
        if (not entity_set.contains(entity)) {
          forEachLeafNode(entityLocalView(), [&](auto& leaf_local_space, auto) {
            leaf_local_space.indices().resize(0);
          });
          return;
        }

        // Here we transform (if necessary) the assembly entity into the entity
        // used for the ordering. This allows orderings to have a different
        // entity types than the grid used for assembly
        // (e.g. dune-multidomaingrid or dune-grid-glue) Both types, however,
        // still need to have the same dimension
        decltype(auto) casted_entity = entityCast(entity_set, std::forward<AssemblyEntity>(entity));
        if constexpr (std::is_rvalue_reference_v<decltype(casted_entity)&&>) {
          // the caller/cast assigned the ownership to us
          _entity_store.emplace(std::move(casted_entity));
          _entity_view = &(_entity_store.value());
        } else {
          // ownership is managed elsewhere
          _entity_view = &casted_entity;
        }

        // cache gt_index and entity_index for used sub-entities (unrolled loop)
        // we cache them for several reasons:
        //  * reuse them when the tree below is more than one node
        //  * unrolling the codim loop allows to have a constexpr codim (countrary to local key's codim)
        //  * access them in a more regular pattern (local kyes don't offer any guarantee)
        //  * lock/unlock sub-entity mutexes
        if (_entity_set_ordering.disjointCodimClosure()) {
          // DG/FV branch
          SizeType gt_index = GlobalGeometryTypeIndex::index(_entity_view->type());
          SizeType entity_index = entity_set.indexSet().index(*_entity_view);
          _index_cache[0] = {gt_index, entity_index};
        } else {
          std::fill(begin(_codim_offsets), end(_codim_offsets), 0);
          _index_cache.clear();
          Dune::Hybrid::forEach(std::make_index_sequence<Element::dimension+1>{}, [&](auto codim) {
            if (codim == 0 or (EntitySetOrdering::mayContainCodim(codim) and _entity_set_ordering.containsCodim(codim))) {
              std::size_t sub_entities = _entity_view->subEntities(codim);
              _index_cache.resize(_index_cache.size() + sub_entities);
              std::size_t codim_offset = _codim_offsets[codim];
              for (std::size_t s = 0; s != sub_entities; ++s) {
                const auto& sub_entity = _entity_view->template subEntity<codim>(s);
                SizeType gt_index = GlobalGeometryTypeIndex::index(sub_entity.type());
                SizeType entity_index = entity_set.indexSet().index(sub_entity);
                _index_cache[codim_offset + s] = {gt_index, entity_index};
              }
            }
            _codim_offsets[codim+1] = _index_cache.size();
          });
        }

        entityLocalView().bindElement(_entity_view);


        forEachLeafNode(entityLocalView(), [&]<class LeafLocalView>(LeafLocalView& leaf_local_view, const auto& suffix) {
          // Here is the heart of binding:
          //   we need to fill the indices for the local space
          auto& indices = leaf_local_view.indices();
          if (not leaf_local_view.boundElement()) {
            indices.clear();
            return;
          }

          constexpr std::size_t fem_codim = EntitySet::dimension - LeafLocalView::dimDomain;
          auto ordering_path = join(_sub_space_path, leaf_local_view.orderingViewPath());
          const auto& leaf_ordering = PDELab::containerEntry(*_ordering, ordering_path);
          const auto& fem = leaf_ordering.space().finiteElementMap();
           // TODO: document that finite element map must be thread-safe!
          leaf_local_view.bindFiniteElement(fem.find(leaf_local_view.element()));
          const auto& fe = leaf_local_view.finiteElement();
          assert(fe.type() == leaf_local_view.element().type());

          indices.resize(fe.size());

          using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LeafLocalView::FiniteElement>;
          if (leaf_local_view.disjointCodimClosure()) {
            // DG/FV branch
            const auto [gt_index, entity_index] = _index_cache[0];
            assert(fe.size() == leaf_ordering.blockCount(gt_index, entity_index));
            // In this case, we know that indices are contiguous and can be
            //   reconstructed from the local index. Thus, there is no need to
            //   fill the whole vector.
            if (not indices.empty())
              indices[0] = _ordering->firstContainerIndex(ordering_path, gt_index, entity_index);
            // notice that the other indices are garbage!
            // but we need them for the size to be correct
            if constexpr (not std::decay_t<LeafLocalView>::optimizeFastDG()) {
              // This the observed behavior for non-optimizaed cases:
              const auto& coeffs = FESwitch::coefficients(fe);
              assert(coeffs.localKey(0).index() == 0 && "Indices in DG/FV elements shall be ordered");
              for (std::size_t dof = 1; dof < indices.size(); ++dof) {
                assert(coeffs.localKey(dof).index() == dof && "Indices in DG/FV elements shall be ordered");
                indices[dof] = accumulate_back(indices[0], dof);
              }
            }
          } else {
            const auto& coeffs = FESwitch::coefficients(fe);
            // for each local dof, we need to figure out its assigned container index
            for (std::size_t dof = 0; dof < coeffs.size(); ++dof) {
              const auto& key = coeffs.localKey(dof);
              assert(_entity_set_ordering.containsCodim(fem_codim + key.codim()));
              const auto [gt_index, entity_index] = [&]{
                if constexpr (fem_codim == 0) {
                  auto offset = _codim_offsets[key.codim()];
                  offset += key.subEntity();
                  assert(offset < _index_cache.size());
                  return _index_cache[offset];
                } else {
                  // in this case the cache does not work because we do not know the embedding of the sub entity
                  auto fem_sub_entity = back(suffix);
                  const auto& sub_entity = _entity_view->template subEntity<fem_codim>(fem_sub_entity);
                  const auto& ref_el = referenceElement(sub_entity.geometry());
                  auto gt_index = GlobalGeometryTypeIndex::index(ref_el.type(key.subEntity(), key.codim()));
                  auto index = entity_set.indexSet().subIndex(sub_entity, key.subEntity(), fem_codim + key.codim());
                  return std::make_tuple(gt_index, index);
                }
              }();
              auto base_index = _ordering->firstContainerIndex(ordering_path, gt_index, entity_index);
              assert(key.index() < leaf_ordering.blockCount(gt_index, entity_index));
              indices[dof] = accumulate_back(base_index, key.index());
            }
          }
        });
      }

      void unbind() noexcept {
        _entity_view = nullptr;
        forEachLeafNode(entityLocalView(), [](auto& leaf_local_view) {
          leaf_local_view.unbind();
        });
      }

      void lock() noexcept {
        lock(_index_cache);
      }

      bool try_lock() noexcept {
        return try_lock(_index_cache);
      }

      void unlock() noexcept {
        unlock(_index_cache);
      }

    private:

      void lock(const std::vector<std::array<SizeType,2>>& indices) noexcept {
        // spin lock until we adquire all of the mutexes
        auto [gt_index0, entity_index0] = indices[0];
        if (_entity_set_ordering.disjointCodimClosure()) {
          _entity_set_ordering.entityLockHandle(gt_index0,entity_index0).lock();
        } else {
          while (not tryLockCodimClosure(indices)) {
            _entity_set_ordering.entityLockHandle(gt_index0,entity_index0).wait();
          }
        }
      }

      bool try_lock(const std::vector<std::array<SizeType,2>>& indices) noexcept {
        if (_entity_set_ordering.disjointCodimClosure()) {
          auto [gt_index, entity_index] = indices[0];
          return _entity_set_ordering.entityLockHandle(gt_index,entity_index).try_lock();
        } else {
          return tryLockCodimClosure(indices);
        }
      }

      void unlock(const std::vector<std::array<SizeType,2>>& indices) noexcept {
        if (_entity_set_ordering.disjointCodimClosure()) {
          auto [gt_index, entity_index] = indices[0];
          _entity_set_ordering.entityLockHandle(gt_index,entity_index).unlock();
        } else {
          for (auto [gt_index, entity_index] : indices) {
            _entity_set_ordering.entityLockHandle(gt_index,entity_index).unlock();
          }
        }
      }

      [[nodiscard]] bool tryLockCodimClosure(const std::vector<std::array<SizeType,2>>& indices) noexcept {
        for (std::size_t i = 0; i < indices.size(); ++i) {
          // try to lock every index on the vector
          if (not _entity_set_ordering.entityLockHandle(indices[i][0],indices[i][1]).try_lock()) [[unlikely]] {
            // entity was already locked, we have to roll back
            for (std::size_t j = i; j != 0; --j)
              _entity_set_ordering.entityLockHandle(indices[j-1][0],indices[j-1][1]).unlock();
            // ...and inform that we could not adquire the lock
            return false;
          }
        }
        // if all entities were locked, we succeded
        return true;
      }

      std::array<std::size_t,Element::dimension+2> _codim_offsets;
      std::vector<std::array<SizeType,2>> _index_cache;
      std::shared_ptr<GlobalOrdering> _ordering;
      std::optional<const Element> _entity_store;
      Element const * _entity_view;
      EntitySetOrdering& _entity_set_ordering;
      SubSpacePath _sub_space_path;
    };

  public:

    template<Concept::TreeNode GlobalOrdering, Concept::MultiIndex Prefix, Concept::MultiIndex SubSpacePath>
    auto makeLocalIndexSet(const std::shared_ptr<GlobalOrdering>& ordering, const Prefix& prefix, const SubSpacePath& sub_space_path) const {
      auto elis = entityOrdering().makeLocalIndexSet(ordering, prefix, sub_space_path);
      using EntityLocalIndexSet = std::decay_t<decltype(*elis)>;
      using LIndexSet = LocalIndexSet<EntityLocalIndexSet, GlobalOrdering, SubSpacePath>;
      EntitySetOrdering& entity_set_ordering = PDELab::containerEntry(*ordering, prefix);
      assert(&entity_set_ordering == this);
      return std::make_unique<LIndexSet>(std::move(*elis), ordering, entity_set_ordering, sub_space_path);
    }

    template<Concept::TreeNode GlobalOrdering, Concept::MultiIndex Prefix, Concept::MultiIndex SubSpacePath>
    auto makeLocalView(const std::shared_ptr<GlobalOrdering>& ordering, const Prefix& prefix, const SubSpacePath& sub_space_path) const {
      auto elv = entityOrdering().makeLocalView(ordering, prefix, sub_space_path);
      using EntityLocalView = std::decay_t<decltype(*elv)>;
      using LSpace = LocalView<EntityLocalView, GlobalOrdering, SubSpacePath>;
      EntitySetOrdering& entity_set_ordering = PDELab::containerEntry(*ordering, prefix);
      assert(&entity_set_ordering == this);
      return std::make_unique<LSpace>(std::move(*elv), ordering, entity_set_ordering, sub_space_path);
    }

    template<class Traits>
    static constexpr auto makeVectorContainer() {
      using EntityOrderingContainer = decltype(EntityOrdering::template makeVectorContainer<Traits>());
      if constexpr (containerBlocked()) {
        return Traits::template makeVector<EntityOrderingContainer>();
      } else {
        using MergedBlock = typename Traits::template block_type<EntityOrderingContainer>::type;
        return Traits::template makeVector<MergedBlock>();
      }
    }

  private:
    // These may be borrowed from the local ordering (base class). But notice
    // that all nodes are destroyed top to bottom, thus, we do not have the risk
    // of a memory leak.
    std::shared_ptr<std::vector<SizeType>> _gt_dof_offsets;
    std::shared_ptr<std::vector<SizeType>> _gt_entity_offsets;
    std::shared_ptr<std::vector<SizeType>> _entity_dof_offsets;
    std::vector<VectorLockType> _gt_mutex;
  };

  template<Concept::Impl::PreBasisTree PreBasis>
  auto makeEntitySetOrdering(const PreBasis& pre_basis) {
    auto entity_ordering = pre_basis.makeLocalOrdering();
    using LocalOrdering = std::decay_t<decltype(*entity_ordering)>;
    // TODO: assert every leaf has same entity set
    using Ordering = EntitySetOrdering<LocalOrdering,typename PreBasis::Traits::MergingStrategy>;
    return std::make_unique<Ordering>(std::move(*entity_ordering));
  }

} // namespace Dune::PDELab::inline Experimental::Impl

#endif // DUNE_PDELAB_BASIS_ORDERING_ENTITY_SET_HH
