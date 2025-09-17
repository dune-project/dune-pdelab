#ifndef DUNE_PDELAB_COMMON_PARTITION_HALO_REGION_HH
#define DUNE_PDELAB_COMMON_PARTITION_HALO_REGION_HH

#include <dune/common/exceptions.hh>

#include <type_traits>
#include <limits>

namespace Dune::PDELab::EntitySetPartition {

//! Possible halo regions of an entity in a partition
enum class HaloRegion {Interior, Overlap};

//! Constant for interior halo region
inline static constexpr std::integral_constant<HaloRegion, HaloRegion::Interior> interior_halo_region = {};

//! Constant for overlap halo region
inline static constexpr std::integral_constant<HaloRegion, HaloRegion::Overlap> overlap_halo_region = {};

//! Constant for a halo distance with all halo regions being in the interior
inline static constexpr std::size_t all_interior_halo_region = std::numeric_limits<std::size_t>::max();

//! Constant for a halo distance with all halo regions being in the overlap
inline static constexpr std::size_t all_overlap_halo_region = std::numeric_limits<std::size_t>::max() - 1;

//! Exception for partition errors
struct PartitionError : public Dune::Exception {};

} //namespace Dune::PDELab::EntitySetPartition

#endif // DUNE_PDELAB_COMMON_PARTITION_HALO_REGION_HH
