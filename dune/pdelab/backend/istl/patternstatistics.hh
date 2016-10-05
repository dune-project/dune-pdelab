// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_PATTERNSTATISTICS_HH
#define DUNE_PDELAB_BACKEND_ISTL_PATTERNSTATISTICS_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <iostream>

namespace Dune {
  namespace PDELab {
    namespace ISTL {

      //! Statistics about the pattern of a BCRSMatrix.
      template<typename T>
      class PatternStatistics
      {

      public:

        //! size_type of the associated BCRSMatrix.
        typedef T size_type;

#ifndef DOXYGEN

        PatternStatistics(size_type nnz,
                          size_type longest_row,
                          size_type overflow_count,
                          size_type estimate,
                          size_type rows)
          : _nnz(nnz)
          , _longest_row(longest_row)
          , _overflow_count(overflow_count)
          , _estimate(estimate)
          , _rows(rows)
        {}

#endif

        //! The total number of nonzero entries in the matrix.
        size_type nonZeros() const
        {
          return _nnz;
        }

        //! The maximum number of nonzero entries in any row of the matrix.
        size_type longestRow() const
        {
          return _longest_row;
        }

        //! The number of nonzero entries that had to be temporarily stored in the overflow area during pattern construction.
        size_type overflowCount() const
        {
          return _overflow_count;
        }

        //! The estimated number of nonzeros per row as provided by the user before pattern construction.
        size_type estimatedEntriesPerRow() const
        {
          return _estimate;
        }

        //! The number of matrix rows.
        size_type rows() const
        {
          return _rows;
        }

        //! The average number of nonzero entries per row, after matrix construction was completed.
        double averageEntriesPerRow() const
        {
          return static_cast<double>(_nnz) / _rows;
        }

        friend std::ostream& operator<<(std::ostream& os, const PatternStatistics& s)
        {
          std::cout << "==== Pattern statistics ====" << std::endl
                    << "matrix rows: " << s.rows() << std::endl
                    << "nonzero entries: " << s.nonZeros() << std::endl
                    << "maximum number of nonzeros per row: " << s.longestRow() << std::endl
                    << "user-provided estimate of nonzeros per row: " << s.estimatedEntriesPerRow() << std::endl
                    << "average nonzeros per row: " << s.averageEntriesPerRow() << std::endl
                    << "number of entries in overflow area during setup: " << s.overflowCount() << std::endl;
          return os;
        }

      private:

        size_type _nnz;
        size_type _longest_row;
        size_type _overflow_count;
        size_type _estimate;
        size_type _rows;

      };

    } // namespace ISTL
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_PATTERNSTATISTICS_HH
