#ifndef DUNE_PDELAB_COMMON_LOCAL_CONTAINER_HH
#define DUNE_PDELAB_COMMON_LOCAL_CONTAINER_HH

#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/concepts/local_index_set.hh>
#include <dune/pdelab/concepts/container.hh>

#include <dune/pdelab/common/local_container_entry.hh>

#include <dune/typetree/treecontainer.hh>
#include <dune/typetree/childextraction.hh>

#include <cassert>
#include <functional>
#include <utility>
#include <vector>
#include <mutex>
#include <atomic>
#if !__cpp_lib_atomic_ref && __has_include(<boost/atomic/atomic_ref.hpp>)
#include <boost/atomic/atomic_ref.hpp>
#endif

namespace Dune::PDELab::inline Experimental {

// this class provides a buffer to store intermediate results of local computations
// and provides the means to gather and scatter to the global container without data races

// TODO: this interface does not let you know if a vector was not loaded!
template<Concept::Basis Basis_, Concept::Container<Basis_> Container>
class LocalContainerBuffer
{
  template<class LeafNode>
  struct LeafToField {
    using Field = std::decay_t<decltype(localContainerEntry(std::declval<Container&>(), std::declval<LeafNode>(), 0))>;
  };

  using LocalTreeContainer = TypeTree::UniformTreeContainer<std::size_t, typename Basis_::LocalView::Tree>;
  using LocalConstraints = typename Basis_::LocalConstraints;
public:

  using Basis = Basis_;
  using Weight = int;

  LocalContainerBuffer(const Basis& space, Container& container)
    : _lconstraints{ space.localConstraints() }
    , _container{ std::ref(container) }
  {
    auto it = space.entitySet().template begin<0>();
    if (space.entitySet().size(0) == 0)
      return;
    auto lspace = space.localView();
    lspace.bind(*it);
    reserve(lspace);
  }

  LocalContainerBuffer(const Concept::LocalIndexSet auto& lspace, Container& container)
    : _lconstraints{ lspace.globalBasis().localConstraints() }
    , _container{ std::ref(container) }
  {
    reserve(lspace);
  }

  LocalContainerBuffer(const LocalContainerBuffer& other)
  : _buffer{other._buffer}
  , _offsets{other._offsets}
  , _lconstraints{other._lconstraints}
  , _container{other._container}
  {}

  LocalContainerBuffer& operator=(const LocalContainerBuffer& other) {
    _buffer = other._buffer;
    _offsets = other._offsets;
    _lconstraints = other._lconstraints;
    _container = other._container;
    return *this;
  }

  void clear(const Concept::LocalIndexSet auto& lspace) noexcept {
    resize(lspace);
  }

  void accumulate(const Concept::LeafTreeNode auto& node, auto dof, auto val) {
    data(node.path())[dof] += val;
  }

  [[nodiscard]] static constexpr Weight weight() noexcept { return 1; }

  friend void swap(LocalContainerBuffer& lhs, LocalContainerBuffer& rhs)
  {
    using std::swap;
    swap(lhs._buffer, rhs._buffer);
    swap(lhs._offsets, rhs._offsets);
    swap(lhs._container, rhs._container);
    swap(lhs._lconstraints, rhs._lconstraints);
  }

  void load(const Concept::LocalBasis auto& lspace, std::convertible_to<bool> auto is_correction) noexcept
  {
    this->clear(lspace);
    _lconstraints.bind(lspace.element());
    forEachLeafNode(lspace.tree(), [&](const auto& lspace_node, auto path) {
      const auto& lconstraints_node = PDELab::containerEntry(_lconstraints.tree(), path);
      auto data_ptr = data(lspace_node.path());
      for (std::size_t dof = 0; dof < lspace_node.size(); ++dof) {
        if (lconstraints_node.isConstrained(dof)) {
          if (not is_correction)
            data_ptr[dof] = lconstraints_node.translationValue(dof); // TODO handle units!
          for (auto [ci, weight] : lconstraints_node.linearCoefficients(dof))
            data_ptr[dof] += weight * containerEntry(_container.get(), lspace_node.index(dof));
        } else {
          data_ptr[dof] = localContainerEntry(_container.get(), lspace_node, dof);
        }
      }
    });
    _lconstraints.unbind();
  }

  void load(const Concept::LocalBasis auto& lspace) noexcept
  {
    load(lspace, std::false_type());
  }

  void fetch_add(Concept::LocalBasis auto& lspace, std::convertible_to<bool> auto is_correction) noexcept {
#if !(__cpp_lib_atomic_ref || __has_include(<boost/atomic/atomic_ref.hpp>))
    DUNE_THROW(NotImplemented, "To enable this feature std::atomic_ref or boost::atomic_ref is required");
#elif __cpp_lib_atomic_ref
    using std::atomic_ref;
#else
    using boost::atomic_ref;
#endif
    forEachEntryStore(lspace, is_correction,
      [](auto& lhs, auto rhs){lhs += rhs;},
      [](auto& lhs, auto rhs){atomic_ref(lhs).fetch_add(rhs, std::memory_order::relaxed);}
    );
  }

  void fetch_add(Concept::LocalBasis auto& lspace) noexcept {
    fetch_add(lspace, std::false_type());
  }

  void store(Concept::LocalBasis auto& lspace, std::convertible_to<bool> auto is_correction) noexcept {
#if !(__cpp_lib_atomic_ref || __has_include(<boost/atomic/atomic_ref.hpp>))
    DUNE_THROW(NotImplemented, "To enable this feature std::atomic_ref or boost::atomic_ref is required");
#elif __cpp_lib_atomic_ref
    using std::atomic_ref;
#else
    using boost::atomic_ref;
#endif
    forEachEntryStore(lspace, is_correction,
      [](auto& lhs, auto rhs){lhs = rhs;},
      [](auto& lhs, auto rhs){atomic_ref(lhs).store(rhs, std::memory_order::relaxed);}
    );
  }

  void store(Concept::LocalBasis auto& lspace) noexcept {
    store(lspace, std::false_type());
  }

  [[nodiscard]] const auto& operator()(const Concept::LeafTreeNode auto& node, auto dof) const noexcept {
    return data(node.path())[dof];
  }

  [[nodiscard]] auto& operator()(const Concept::LeafTreeNode auto& node, auto dof) noexcept {
    return data(node.path())[dof];
  }

  // get a pointer to the begining of const data for a given node
  template<Concept::MultiIndex Path>
  [[nodiscard]] auto data(Path path) const noexcept {
    using Node = TypeTree::ChildForTreePath<typename Basis::LocalView::Tree, Path>;
    using Field = typename LeafToField<Node>::Field;
    return std::launder(reinterpret_cast<Field const *>(_buffer.data() + _offsets[path]));
  }

  // get a pointer to the begining of data for a given node
  template<Concept::MultiIndex Path>
  [[nodiscard]] auto data(Path path) noexcept {
    using Node = TypeTree::ChildForTreePath<typename Basis::LocalView::Tree, Path>;
    using Field = typename LeafToField<Node>::Field;
    return std::launder(reinterpret_cast<Field*>(_buffer.data() + _offsets[path]));
  }

  [[nodiscard]] auto operator[](Concept::MultiIndex auto path) const noexcept {
    return data(path);
  }

  [[nodiscard]] auto operator[](Concept::MultiIndex auto path) noexcept {
    return data(path);
  }

private:

  // prepare buffer to hold maximum data size
  void reserve(const Concept::LocalIndexSet auto& lspace) {
    _offsets.resize(lspace.tree());
    std::size_t max_field_size = 0;
    std::size_t alignement_padding = 0;
    forEachLeafNode(lspace.tree(), [&]<class Node>(const Node& lspace_node) {
      using Field = typename LeafToField<Node>::Field;
      max_field_size = std::max<std::size_t>(max_field_size, sizeof(Field));
      alignement_padding += alignof(Field) - 1;
    });
    _buffer.resize(alignement_padding + max_field_size * lspace.maxSize());
  }

  // resizes and initializes data. Additionally, stores offsets to begin of data for each node
  void resize(const Concept::LocalIndexSet auto& lspace) noexcept {
    void* head = _buffer.data();
    std::size_t space = _buffer.size();
    forEachLeafNode(lspace.tree(), [&]<class Node>(const Node& lspace_node) {
      using Field = typename LeafToField<Node>::Field;
      std::size_t size = lspace_node.size();
      head = std::align(alignof(Field), size*sizeof(Field), head, space);
      // If assert fails, we did not reserve enough space for whatever reason.
      // Because the user passed the wrong local space, or because we reserved the wrong size.
      // In either case, it's a bug.
      assert(head != nullptr && "The reserved buffer cannot fit the local data");
      // Notice that next line initializes objects liftime with placement new, however,
      // the (non-UB) pointers are discarded. Thus, to avoid confusing
      // the compiler about the lifetime of the new objects we need to "launder"
      // their lifetime on every new access out of the old pointers, otherwise, we
      // would trigger UB. Also note that laundering extends to "reachable" objects
      // (i.e., bracket operator). See data() member function.
      auto tail = std::uninitialized_fill_n(reinterpret_cast<Field*>(head), size, Field{0});
      // We do not call destructor in each field. Instead, we end the liftime of the object
      // prematurely. This only well defined as long as the object is trivially destructible
      static_assert(std::is_trivially_destructible_v<Field>);
      // let's store an offset wrt the buffer where this node data begins
      _offsets[lspace_node.path()] = std::distance(_buffer.data(), static_cast<std::byte*>(head));
      // advance iterator to next position
      head = tail;
      space -= size*sizeof(Field);
    });
  }

  // stores and destroys data in buffer
  template<Concept::LocalBasis LocalBasis>
  void forEachEntryStore(LocalBasis& lspace, auto is_correction, auto fapply, auto fapply_atomic) noexcept
  {
    _lconstraints.bind(lspace.element());

    auto for_each_entry = [&](auto f) {
      return [&, f](const auto& lspace_node, auto path) {
        const auto& lconstraints_node = PDELab::containerEntry(_lconstraints.tree(), path);
        auto data_ptr = data(lspace_node.path());
        for (std::size_t dof = 0; dof != lspace_node.size(); ++dof) {
          if (lconstraints_node.isConstrained(dof)) {
            if (not is_correction)
              containerEntry(_container.get(), lspace_node.index(dof)) = lconstraints_node.translationValue(dof);
            const auto& linear_coeff = lconstraints_node.linearCoefficients(dof);
            if (not linear_coeff.empty() and lspace.partitionRegion() == EntitySetPartitioner::shared_region)
              DUNE_THROW(NotImplemented, "The concurrency model for non correction methods with generic affine constraints is missing. "
                                          "Consider using a local space without linear constraints or a single threaded assembly.");
            for (auto [ci, weight] : linear_coeff) {
              containerEntry(_container.get(), ci) += weight * data_ptr[dof];
            }
          } else {
            f(localContainerEntry(_container.get(), lspace_node, dof), data_ptr[dof]);
          }
        }
      };
    };

    if constexpr (Concept::BasicLockable<LocalBasis>) {
      auto scope_guard = [&](){
        if (lspace.partitionRegion() == EntitySetPartitioner::shared_region)
          return std::unique_lock{lspace};
        else
          return std::unique_lock{lspace, std::defer_lock};
      }();
      forEachLeafNode(lspace.tree(), for_each_entry(fapply));
    } else if (lspace.partitionRegion() == EntitySetPartitioner::shared_region) {
      forEachLeafNode(lspace.tree(), for_each_entry(fapply_atomic));
    } else {
      forEachLeafNode(lspace.tree(), for_each_entry(fapply));
    }

    _lconstraints.unbind();
  }

  std::vector<std::byte> _buffer;
  LocalTreeContainer _offsets;
  LocalConstraints _lconstraints;
  std::reference_wrapper<Container> _container;
};


template<Concept::LocalMutableContainer LocalContainer, class WeightType>
class WeightedLocalContainerView {
public:
  using Basis = typename LocalContainer::Basis;
  using Weight = decltype(typename LocalContainer::Weight{} * WeightType{});

  WeightedLocalContainerView(LocalContainer& lcontainer, WeightType weight)
    : _weight{weight}
    , _lcontainer{lcontainer}
  {}

  void accumulate(const Concept::LeafTreeNode auto& node, auto dof, auto val) noexcept {
    _lcontainer.accumulate(node, dof, weight() * val);
  }

  [[nodiscard]] const auto& operator()(const Concept::LeafTreeNode auto& node, auto dof) const noexcept {
    return _lcontainer(node, dof);
  }

  [[nodiscard]] auto& operator()(const Concept::LeafTreeNode auto& node, auto dof) noexcept {
    return _lcontainer(node, dof);
  }

  [[nodiscard]] Weight weight() const noexcept {
    return _lcontainer.weight() * _weight;
  }

  [[nodiscard]] auto data(Concept::MultiIndex auto path) const noexcept {
    return _lcontainer.data(path);
  }

  [[nodiscard]] auto data(Concept::MultiIndex auto path) noexcept {
    return _lcontainer.data(path);
  }

private:
  WeightType _weight;
  LocalContainer& _lcontainer;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_LOCAL_CONTAINER_HH
