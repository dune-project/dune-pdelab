#ifndef DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH
#define DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH

//! maximum grid world dimension
#ifndef DUNE_PDELAB_MAX_WORLDDIM
#define DUNE_PDELAB_MAX_WORLDDIM 32
#endif

#include <dune/grid/concepts/gridview.hh>

#include <dune/pdelab/basis/ordering/lexicographic.hh>
#include <dune/pdelab/basis/ordering/entityset.hh>
#include <dune/pdelab/basis/ordering/entity.hh>
#include <dune/pdelab/basis/prebasis/concept.hh>

#include <bitset>

namespace Dune::PDELab::inline Experimental::inline Strategy {

template<bool ContainerBlocked>
struct DefaultStrategy
{
  using SizeType = std::size_t;
  using CodimFlag = std::bitset<DUNE_PDELAB_MAX_WORLDDIM>;
  static constexpr bool Blocked = ContainerBlocked;
};

template<bool ContainerBlocked>
void registerIndexMergingStrategy(DefaultStrategy<ContainerBlocked>);

template<bool Blocked>
struct Lexicographic : public DefaultStrategy<Blocked> {

  static auto makeOrdering(const Concept::Impl::PreBasisTree auto& pre_basis) {
    return Impl::makeLexicographicOrdering(pre_basis);
  }
};

using FlatLexicographic = Lexicographic<false>;
using BlockedLexicographic = Lexicographic<true>;

constexpr FlatLexicographic
flatLexicographic()
{
  return {};
}
constexpr BlockedLexicographic
blockedLexicographic()
{
  return {};
}


template<Dune::Concept::GridView ES, bool ContainerBlocked>
struct EntityGrouping : public DefaultStrategy<ContainerBlocked>
{
  using EntitySet = ES;

  EntityGrouping(EntitySet entity_set)
    : _entity_set{ std::move(entity_set) }
  {
  }

  EntitySet& entitySet() { return _entity_set; }
  const EntitySet& entitySet() const { return _entity_set; }

  static auto makeOrdering(const Concept::Impl::PreBasisTree auto& pre_basis) {
    return Impl::makeEntitySetOrdering(pre_basis);
  }

  static auto makeLocalOrdering(const Concept::Impl::PreBasisTree auto& pre_basis) {
    return Impl::makeEntityOrdering(pre_basis);
  }

private:
  EntitySet _entity_set;
};

template<Dune::Concept::GridView EntitySet>
constexpr auto
flatByEntity(const EntitySet& entity_Set)
{
  return EntityGrouping<EntitySet, false>{ entity_Set };
}

template<Dune::Concept::GridView EntitySet>
constexpr auto
blockedByEntity(const EntitySet& entity_Set)
{
  return EntityGrouping<EntitySet, true>{ entity_Set };
}

} // namespace Dune::PDELab::inline Experimental::inline Strategy

#endif // DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH
