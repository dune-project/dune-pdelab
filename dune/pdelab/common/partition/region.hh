#ifndef DUNE_PDELAB_COMMON_PARTITION_REGION_HH
#define DUNE_PDELAB_COMMON_PARTITION_REGION_HH

#include <type_traits>

namespace Dune::PDELab::inline Experimental::EntitySetPartitioner {

enum class Region {Unique, Shared};

inline static constexpr std::integral_constant<Region, Region::Unique> unique_region = {};
inline static constexpr std::integral_constant<Region, Region::Shared> shared_region = {};

} //namespace Dune::PDELab::inline Experimental::EntitySetPartitioner

#endif // DUNE_PDELAB_COMMON_PARTITION_REGION_HH
