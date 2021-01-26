#ifndef DUNE_PDELAB_COMMON_CONCEPTS_HH
#define DUNE_PDELAB_COMMON_CONCEPTS_HH


#include <utility>

#include <dune/common/concept.hh>
#include <dune/functions/functionspacebases/concepts.hh>

namespace Dune::PDELab {

namespace Concept {

using namespace Dune::Concept;

struct GridFunctionSpace
{
  template<class GFS>
  auto require(GFS&& gfs) -> decltype(
    requireType<typename GFS::NodeTag>(),
    requireType<typename GFS::Traits::Backend>(),
    requireConvertible<typename GFS::Traits::GridView>(gfs.gridView()),
    gfs.entitySet(),
    // only on the leafs
    // gfs.finiteElementMap(),
    // gfs.finiteElementMapStorage(),
    // gfs.constraints(),
    // gfs.constraintsStorage(),
    // requireConvertible<typename GFS::Traits::FiniteElementMap>(gfs.finiteElementMap()),
    // requireConvertible<typename GFS::Traits::ConstraintsType>(gfs.constraints()),
    requireConvertible<typename GFS::Ordering>(gfs.ordering())
    // requireConvertible<std::shared_ptr<typename GFS::Ordering>>(gfs.orderingStorage())
  );
};

struct BasisInfo // check for our thin wrapper around a dune-functions basis
{
  template<class B>
  auto require(B&& t) -> decltype(
    requireType<typename B::Basis>(),
    requireType<typename B::GridView>(),
    requireType<typename B::ContainerIndex>(),
    requireType<typename B::Traits::Backend>(),
    requireConcept<Dune::Functions::Concept::GlobalBasis<typename B::GridView>>(t.basis())
  );
};

} // namespace Dune::PDELab::Concept

/// Check if B models a PDELab GridFrunctionSpace
template<class B>
static constexpr bool isGridFunctionSpace()
{
  return models<Concept::GridFunctionSpace, B>();
}

/// Check if B is our thin wrapper around a dune-functions basis
template<class B>
static constexpr bool isBasisInfo()
{
  return models<Concept::BasisInfo, B>();
}

} // namespace Dune::PDELab


#endif // DUNE_PDELAB_COMMON_CONCEPTS_HH
