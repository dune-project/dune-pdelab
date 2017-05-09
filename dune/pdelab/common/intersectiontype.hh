#ifndef DUNE_PDELAB_COMMON_INTERSECTIONTYPE_HH
#define DUNE_PDELAB_COMMON_INTERSECTIONTYPE_HH

#include <utility>
#include <tuple>

#include <dune/common/version.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune {
  namespace PDELab {

    //! Enum describing the type of an intersection.
    enum class IntersectionType
    {

      processor = 0, //!< processor boundary intersection (neighbor() == false && boundary() == false) or outside entity not in EntitySet
      skeleton = 1,  //!< skeleton intersection (neighbor() == true && boundary() == false)
      boundary = 2,  //!< domain boundary intersection (neighbor() == false && boundary() == true)
      periodic = 3   //!< periodic boundary intersection (neighbor() == true && boundary() == true)

    };

    //! Classifies the type of an intersection wrt to the passed EntitySet.
    /**
     * This function classifies the neighborship type of an intersection. It mirrors
     * the classification provided by the Intersection::boundary() and
     * Intersection::neighbor() methods, but also takes into account an EntitySet that
     * might not span all parallel PartitionTypes.
     *
     * As the classification often requires the function to obtain the outside entity, the function
     * returns a tuple containing both the IntersectionType and the outside entity to avoid obtaining
     * the outside entity twice. In case of a boundary intersection where an outside entity is not available,
     * it will return a default-constructed (invalid) entity.
     */
    template<typename EntitySet, typename Intersection>
    std::tuple<IntersectionType,typename EntitySet::Element> classifyIntersection(const EntitySet& entity_set, const Intersection& is)
    {
      auto type = static_cast<IntersectionType>(1* is.neighbor() + 2*is.boundary());
      if (type == IntersectionType::skeleton || type == IntersectionType::periodic)
#if DUNE_VERSION_NEWER_REV(DUNE_GRID,2,4,1)
        if (entity_set.partitions() == Partitions::all)
#else
        if (entity_set.partitions().partitionIterator() == Partitions::all.partitionIterator())
#endif
          return std::make_tuple(type,is.outside());
        else
          {
            auto outside_entity = is.outside();
            if (entity_set.partitions().contains(outside_entity.partitionType()))
              return std::make_tuple(type,outside_entity);
            else
              return std::make_tuple(IntersectionType::processor,std::move(outside_entity));
          }
      else
        return std::make_tuple(type,decltype(is.outside()){});
    }


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_INTERSECTIONTYPE_HH
