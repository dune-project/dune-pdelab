#ifndef DUNE_PDELAB_COMMON_CONCEPTS_HH
#define DUNE_PDELAB_COMMON_CONCEPTS_HH


#include <utility>

#include <dune/common/concept.hh>
#include <dune/functions/functionspacebases/concepts.hh>

namespace Dune::PDELab {

// /// Check if GFS models a PDELab GridFrunctionSpace
// template<class GFS>
// static constexpr bool isGridFunctionSpace()
// {
//   return models<Concept::GridFunctionSpace<Args&&...>, F>();
// }

namespace Concept {

using namespace Dune::Concept;

struct HasGridViewMethod
{
  template<class C>
  auto require(C&& c) -> decltype(
    c.gridView()
  );
};

struct HasFiniteElementMapMethod
{
  template<class C>
  auto require(C&& c) -> decltype(
    c.finiteElementMap()
  );
};

struct HasConstraintsMethod
{
  template<class C>
  auto require(C&& c) -> decltype(
    c.constraints()
  );
};

struct HasOrderingMethod
{
  template<class C>
  auto require(C&& c) -> decltype(
    c.ordering()
  );
};


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
    requireConvertible<typename GFS::Ordering>(gfs.ordering()),
    requireConvertible<std::shared_ptr<typename GFS::Ordering>>(gfs.orderingStorage())
  );
};

// derived from TypeTree::LeafNode
// type ConstraintsContainer
// methods
//    gridView
//    entitySet
//    finiteElementMap
//    finiteElementMapStorage
//    constraints
//    constraintsStorage
//    ordering
//    orderingStorage


} // namespace Dune::PDELab::Concept
} // namespace Dune::PDELab


#endif // DUNE_PDELAB_COMMON_CONCEPTS_HH
