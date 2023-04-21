#ifndef DUNE_PDELAB_BASIS_ORDERING_HH
#define DUNE_PDELAB_BASIS_ORDERING_HH

#include <dune/pdelab/basis/prebasis/concept.hh>
#include <dune/pdelab/basis/ordering/entity.hh>
#include <dune/pdelab/basis/ordering/entityset.hh>
#include <dune/pdelab/basis/ordering/lexicographic.hh>
#include <dune/pdelab/basis/merging_strategy.hh>

namespace Dune::PDELab::Impl {

  /**
   * @brief Make an ordering pointer for a given discrete function spaces tree
   * @note This function is inteded to be a customization point for spaces that
   * want a custom merging strategy and that produce an ordering. This is the
   * place to introduce ordering transformations over indices such as chunking,
   * mappings or simple change of multi-index types.
   *
   * @tparam PreBasis   A discrete function space tree
   * @param pre_basis   The spaces for which to create an ordering
   * @return auto     A pointer to an ordering usable by Basis
   */
  template <Concept::Impl::PreBasisTree PreBasis>
  auto makeOrdering(const PreBasis &pre_basis)
  {
    return makeOrdering(pre_basis, pre_basis.mergingStrategy());
  }

} // namespace Dune::PDELab::Impl

#endif // DUNE_PDELAB_BASIS_ORDERING_HH
