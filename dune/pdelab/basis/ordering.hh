#ifndef DUNE_ASSEMBLER_SPACE_ORDERING_HH
#define DUNE_ASSEMBLER_SPACE_ORDERING_HH

#include <dune/assembler/space/concept.hh>
#include <dune/assembler/space/ordering/entity.hh>
#include <dune/assembler/space/ordering/entityset.hh>
#include <dune/assembler/space/ordering/lexicographic.hh>
#include <dune/assembler/space/merging_strategy.hh>

namespace Dune::Assembler::Impl {

  /**
   * @brief Make an ordering pointer for a given discrete function spaces tree
   * @note This function is inteded to be a customization point for spaces that
   * want a custom merging strategy and that produce an ordering. This is the
   * place to introduce ordering transformations over indices such as chunking,
   * mappings or simple change of multi-index types.
   *
   * @tparam Space  A discrete function space tree
   * @param space   The spaces for which to create an ordering
   * @return auto   A pointer to an ordering usable by DiscreteFunctionSpace
   */
  template <Concept::Impl::SpaceTree Space>
  auto makeOrdering(const Space &space)
  {
    return makeOrdering(space, space.mergingStrategy());
  }

} // namespace Dune::Assembler::Impl

#endif // DUNE_ASSEMBLER_SPACE_ORDERING_HH
