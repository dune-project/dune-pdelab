// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_PARTITIONINFOPROVIDER_HH
#define DUNE_PDELAB_COMMON_PARTITIONINFOPROVIDER_HH

#include <bitset>

#include <dune/grid/common/gridenums.hh>

namespace Dune {
  namespace PDELab {

    //! Mixin class for providing information about contained grid partitions.
    /**
     * This is a mixin class for orderings providing the common implementation
     * of the Dune::PartitionType query interface. As the number of partition types
     * is fixed, we can easily move the complete implementation into this mixin,
     * only requiring the ordering to update the contained information using the
     * protected API.
     */
    class PartitionInfoProvider
    {

    public:

      //! Returns whether this ordering contains entities with PartitionType partition.
      bool containsPartition(PartitionType partition) const
      {
        return _contained_partitions.test(static_cast<unsigned char>(partition));
      }

      //! Returns the internal representation of the set of contained entities.
      std::bitset<6> containedPartitions() const
      {
        return _contained_partitions;
      }

    protected:

      //! Empties the set of contained partitions.
      void clearPartitionSet()
      {
        _contained_partitions.reset();
      }

      //! Adds all partitions contained in r the set of contained partitions.
      void mergePartitionSet(const PartitionInfoProvider& r)
      {
        _contained_partitions |= r._contained_partitions;
      }

      //! Sets the set of contained partitions to the passed-in value.
      /**
       * \warning This is an internal interface that relies on the internal implementation
       *          of the partition set and may change without notice!
       */
      void setPartitionSet(const std::bitset<6>& partitions)
      {
        _contained_partitions = partitions;
      }

      //! Copies the set of contained partitions from r.
      void setPartitionSet(const PartitionInfoProvider& r)
      {
        _contained_partitions = r._contained_partitions;
      }

      //! Adds the partitions from all PartitionInfoProviders in the range [begin,end).
      /**
       * \note The passed-in iterators may yield both references and pointers to the
       *       PartitionInfoProviders in the range. This feature exists mostly to simplify
       *       implementation of the dynamic ordering base classes, which hold pointers
       *       to their children.
       */
      template<typename It>
      void mergePartitionSets(It begin, It end)
      {
        clearPartitionSet();
        for (; begin != end; ++begin)
          mergePartitionSet(reference(*begin));
      }

    private:

      // ********************************************************************************
      // The following two function are here to make mergePartitionSets() work with both
      // normal iterators and iterators of pointers. It would be nice to just use
      // boost::indirect_iterator, but alas...
      // ********************************************************************************

      static const PartitionInfoProvider& reference(const PartitionInfoProvider& provider)
      {
        return provider;
      }

      static const PartitionInfoProvider& reference(const PartitionInfoProvider* provider)
      {
        return *provider;
      }

      std::bitset<6> _contained_partitions;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_PARTITIONINFOPROVIDER_HH
