// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_SOLVERSTATISTICS_HH
#define DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_SOLVERSTATISTICS_HH

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>


namespace Dune{
  namespace PDELab{

    /** \file solverstatistics.hh
     * \brief Provides a class for collecting statistics on the number of
     * block-solves
     */

    /** \brief Statistics result structure */
    template <typename T>
    struct StatisticsResult {
      T minimum; // Global minimum
      T maximum; // Global maximum
      double avg; // Global average
      double stddev; // Standard deviation
      size_t size; // Total number of invocations
    };

    /** \class SolverStatistics
     * \brief Class for collecting statistics over several invocations
     *
     * Records data for every invocation and provides methods for calculating
     * min/max/avg/stddev over all invocations.
     * This can then be used to calculate statistics on the block-solves
     */
    template <typename T>
    class SolverStatistics {
    public:
      /** \brief Create new instance of class
       *
       * \param[in] comm_ Collective communication object
       */
      SolverStatistics(const Dune::CollectiveCommunication<MPI_Comm>& comm_)
        : data(), comm(comm_) {}

      /** \brief Add new data point
       *
       * \param[in] x Data point to add
       */
      void append(const T x) {
        data.push_back(x);
      }

      /** \brief clear out data
       */
      void clear() {
        data.clear();
      }

      /** \brief Total number of calls
       *
       * Calculates total number of invocations
       */
      const size_t size() const {
        size_t local_size = data.size();
        size_t global_size =  comm.sum(local_size);
        return global_size;
      }

      /** \brief Calculate global average
       *
       * Calculates global average over all processors and invocations
       */
      const double avg() const {
        double s_local = (double) std::accumulate(data.begin(),data.end(),0);
        size_t local_size = data.size();
        double global_size = (double) comm.sum(local_size);
        double s_global = comm.sum(s_local);
        return s_global / global_size;
      }

      /** \brief Calculate standard deviation
       *
       * Calculates standard deviation
       */
      const double stddev() const {
        using std::sqrt;
        double s_local = (double) std::accumulate(data.begin(),data.end(),0);
        double ss_local = (double) std::inner_product(data.begin(),data.end(),
                                                      data.begin(),0);
        size_t local_size = data.size();
        double global_size = (double) comm.sum(local_size);
        double s_global = comm.sum(s_local);
        double ss_global = comm.sum(ss_local);
        return sqrt(1./(global_size-1.)*(ss_global-s_global*s_global/global_size));
      }

      /** \brief Calculate global minimum
       *
       * Calculates global minimum over all processors and invocations
       */
      const T min() const {
        T min_local = *std::min_element(data.begin(),data.end());
        T min_global = comm.min(min_local);
        return min_global;
      }

      /** \brief Calculate global maximum
       *
       * Calculates global maximum over all processors and invocations
       */
      const T max() const {
        T max_local = *std::max_element(data.begin(),data.end());
        T max_global = comm.max(max_local);
        return max_global;
      }

      /** \brief Convert to statistics result
       */
      const StatisticsResult<T> result() const {
        StatisticsResult<T> result;
        result.minimum = min();
        result.maximum = max();
        result.avg = avg();
        result.stddev = stddev();
        result.size = size();
        return result;
      }

    private:
      // \brief local data
      std::vector<T> data;
      // \brief Collective communication object
      const Dune::CollectiveCommunication<MPI_Comm>& comm;
    };

    /** \brief Write statistics result to out stream */
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const StatisticsResult<T>& result) {
      os << "#calls = " << result.size;
      os << ", min = " << std::fixed << result.minimum;
      os << ", avg = " << std::fixed << result.avg;
      os << ", stddev = " << std::fixed << result.stddev;
      os << ", max = " << std::fixed << result.maximum;
      return os;
    }

  } // namespace PDELab
} // namespace Dune

#endif
