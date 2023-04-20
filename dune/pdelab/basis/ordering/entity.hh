#ifndef DUNE_PDELAB_BASIS_ORDERING_ENTITY_HH
#define DUNE_PDELAB_BASIS_ORDERING_ENTITY_HH

#include <dune/pdelab/basis/ordering/entity_leaf.hh>
#include <dune/pdelab/basis/ordering/entity_composite.hh>

namespace Dune::PDELab::inline Experimental::Impl {


//! Returns a suitable Entity Ordering tree for a discrete-function-space tree
template<Concept::Impl::PreBasisTree PreBasis>
auto
makeEntityOrdering(const PreBasis& pre_basis)
{
  using MergingStrategy = typename PreBasis::Traits::MergingStrategy;

  static_assert(requires {typename MergingStrategy::EntitySet;},
    "Merging strategy for this node shall provide an entity set. "
    "Apply the following rules to avoid this error:"
    " * Leaf nodes usually need an entity grouping."
    " * If a node is grouped by entity, all of its children have to be grouped by entity."
  );

  if constexpr (Concept::LeafTreeNode<PreBasis>) {
    using EntityOrdering = LeafEntityOrdering<PreBasis>;
    return std::make_unique<EntityOrdering>(pre_basis);
  } else if constexpr (Concept::ArrayTreeNode<PreBasis>) {
    constexpr std::size_t degree = PreBasis::degree();
    using Child = std::decay_t<decltype(*makeEntityOrdering(pre_basis.child(0)))>;
    using EntityOrdering = ArrayEntityOrdering<MergingStrategy, Child, degree>;
    typename EntityOrdering::NodeStorage storage;
    for (std::size_t i = 0; i < degree; ++i)
      storage[i] = makeEntityOrdering(pre_basis.child(i));
    return std::make_unique<EntityOrdering>(std::move(storage), pre_basis.mergingStrategy());
  } else if constexpr (Concept::VectorTreeNode<PreBasis>) {
    std::size_t degree = pre_basis.degree();
    using Child = std::decay_t<decltype(*makeEntityOrdering(pre_basis.child(0)))>;
    using EntityOrdering = VectorEntityOrdering<MergingStrategy, Child>;
    typename EntityOrdering::NodeStorage storage(degree);
    for (std::size_t i = 0; i < degree; ++i)
      storage[i] = makeEntityOrdering(pre_basis.child(i));
    return std::make_unique<EntityOrdering>(std::move(storage), pre_basis.mergingStrategy());
  } else {
    static_assert(Concept::TupleTreeNode<PreBasis>);
    auto unfold_children = [&](auto... i) {
      using EntityOrdering = TupleEntityOrdering<
        MergingStrategy,
        std::decay_t<decltype(*makeEntityOrdering(pre_basis.child(i)))>...>;
      typename EntityOrdering::NodeStorage storage{ makeEntityOrdering(
        pre_basis.child(i))... };
      return std::make_unique<EntityOrdering>(std::move(storage), pre_basis.mergingStrategy());
    };
    auto indices = std::make_index_sequence<PreBasis::degree()>{};
    return unpackIntegerSequence(unfold_children, indices);
  }
}

} // namespace Dune::PDELab::inline Experimental::Impl

#endif // DUNE_PDELAB_BASIS_ORDERING_ENTITY_HH
