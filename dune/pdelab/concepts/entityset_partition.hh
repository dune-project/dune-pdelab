#ifndef DUNE_PDELAB_CONCEPT_ENTITY_SET_PARTITION_HH
#define DUNE_PDELAB_CONCEPT_ENTITY_SET_PARTITION_HH

#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/gridview.hh>

#include <utility>
#include <ranges>
#include <concepts>

namespace Dune::PDELab::inline Experimental::Concept {

/**
 * @brief Models a partition over a set of entities (grid view)
 * @details The partition is defined as three nested containers that reflect the structure of the partition.
 * The first and outermost partition is a container of LabelSets. Each LabelSet is a container of
 * PatchSets, finally, a PatchSet is a container of grid entities
 * of an specific codimension where entity is contained only once in the whole partition.
 * Patches of a label set are indented to be traversed concurrently as long as the computations
 * relative to other patches respect the halo of the label set.
 *
 * The halo of the partition refers to the entities within a label set connected
 * to any patch set with at at-most N edges on the grid graph.
 * This halo helps to specify whether computations on a given entity are shared with another patch of the same label.
 * In particular, colored labels ensure that the there is no entities in the halo, meaning that computations may be lock-free,
 * whereas uncolored partitions will need some syncronization mechanism at the halo.
 * The halo distance is specified at construction of the partition and is a problem dependent requirement:
 * * Conforming grids with only volume integrals typically need a halo equal to 0.
 * * Hanging node grids typically need a halo of 1
 * * Discontinuos galerkin methods typically need a halo of 1.
 */
template<class ESP>
concept EntitySetPartition = requires(const ESP partition, const typename ESP::Entity& entity) {
  // holds specific kind of entities (codimension is fixed)
  requires Dune::Concept::Entity<typename ESP::Entity>;
  { partition.isHalo(entity) }   -> std::convertible_to<bool>;

  // partition holds the (unpartitioned) entity set
  requires Dune::Concept::GridView<typename ESP::EntitySet>;
  { partition.entitySet() }      -> std::convertible_to<typename ESP::EntitySet>;

  // range over labels in the partition
  requires std::ranges::range<ESP>;
  requires std::convertible_to<std::ranges::range_value_t<ESP>, const typename ESP::LabelSet&>;

  // range over patches in the label
  requires std::ranges::range<typename ESP::LabelSet>;
  requires std::convertible_to<std::ranges::range_value_t<typename ESP::LabelSet>, const typename ESP::PatchSet&>;

  // range over entites in the patch
  requires std::ranges::range<typename ESP::PatchSet>;
  requires std::convertible_to<std::ranges::range_value_t<typename ESP::PatchSet>, const typename ESP::Entity&>;
};

} // end namespace Dune::PDELab::Concept

#endif // DUNE_PDELAB_CONCEPT_ENTITY_SET_PARTITION_HH
