// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_BCRSMATRIXBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTL_BCRSMATRIXBACKEND_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <dune/pdelab/backend/istl/bcrsmatrix.hh>
#include <dune/pdelab/backend/istl/bcrspattern.hh>
#include <dune/pdelab/backend/istl/patternstatistics.hh>

namespace Dune {
  namespace PDELab {
    namespace ISTL {

      // ********************************************************************************
      // infrastructure for deducing the pattern type from row and column orderings
      // ********************************************************************************

      namespace {

        template<typename M, typename RowOrdering, typename ColOrdering, bool pattern>
        struct _build_bcrs_pattern_type
        {
          // we use void as a marker for a nonexisting subpattern
          typedef void type;
        };

        // central TMP that is invoked recursively while traversing the ordering hierarchy
        // and builds up the possibly nested pattern.
        template<typename M, typename RowOrdering, typename ColOrdering>
        struct _build_bcrs_pattern_type<M,RowOrdering,ColOrdering,true>
        {

          // descend into blocks
          typedef typename _build_bcrs_pattern_type<
            typename M::block_type,
            RowOrdering,
            ColOrdering,
            requires_pattern<
              typename M::block_type
              >::value
            >::type BlockOrdering;

          // Leafs -> BCRSPattern, interior nodes -> NestedPattern
          typedef typename std::conditional<
            std::is_same<BlockOrdering,void>::value,
            BCRSPattern<
              RowOrdering,
              ColOrdering
              >,
            NestedPattern<
              RowOrdering,
              ColOrdering,
              BlockOrdering
              >
            >::type type;

        };

        // Wrapper TMP for constructing OrderingBase types from function spaces and for
        // shielding user from recursive implementation
        template<typename M, typename GFSV, typename GFSU, typename Tag>
        struct build_bcrs_pattern_type
        {

          typedef OrderingBase<
            typename GFSV::Ordering::Traits::DOFIndex,
            typename GFSV::Ordering::Traits::ContainerIndex
            > RowOrdering;

          typedef OrderingBase<
            typename GFSU::Ordering::Traits::DOFIndex,
            typename GFSU::Ordering::Traits::ContainerIndex
            > ColOrdering;

          typedef typename _build_bcrs_pattern_type<M,RowOrdering,ColOrdering,requires_pattern<M>::value>::type type;
        };

        // Specialization for forcibly flat backends
        template<typename M, typename GFSV, typename GFSU>
        struct build_bcrs_pattern_type<M,GFSV,GFSU,FlatContainerAllocationTag>
        {
          typedef BCRSPattern<typename GFSV::Ordering, typename GFSU::Ordering> type;
        };


        // leaf BCRSMatrix
        template<typename OrderingV, typename OrderingU, typename Pattern, typename Container, typename StatsVector>
        typename std::enable_if<
          std::is_same<typename Pattern::SubPattern,void>::value
          >::type
        allocate_bcrs_matrix(const OrderingV& ordering_v,
                             const OrderingU& ordering_u,
                             Pattern& p,
                             Container& c,
                             StatsVector& stats)
        {
          c.setSize(ordering_v.blockCount(),ordering_u.blockCount(),0);
          c.setBuildMode(Container::random);

          std::vector<typename Pattern::size_type> row_sizes(p.sizes());

          typename Pattern::size_type nnz = 0;
          typename Pattern::size_type longest_row = 0;

          for (typename Pattern::size_type i = 0; i < c.N(); ++i)
            {
              nnz += row_sizes[i];
              longest_row = std::max(longest_row,row_sizes[i]);
              c.setrowsize(i,row_sizes[i]);
            }
          c.endrowsizes();

          stats.push_back(typename StatsVector::value_type(nnz,longest_row,p.overflowCount(),p.entriesPerRow(),ordering_v.blockCount()));

          for (typename Pattern::size_type i = 0; i < c.N(); ++i)
            c.setIndices(i,p.begin(i),p.end(i));

          // free temporary index storage in pattern before allocating data array in matrix
          p.clear();
          // allocate data array
          c.endindices();
        }


        // ********************************************************************************
        // nested matrix allocation
        // In contrast to the older implementation, this code still uses BCRSMatrix for nested matrices,
        // but we do not attempt to keep the pattern of those interior matrices sparse and always allocate all
        // blocks. That greatly simplifies the code, and as those interior matrices really shouldn't be very
        // large, any performance impact is minimal.
        // The code also collects statistics about the pattern of the leaf BCRSMatrices and collates those stats
        // in a row-major ordering.
        // ********************************************************************************

        // interior BCRSMatrix
        template<typename OrderingV, typename OrderingU, typename Pattern, typename Container, typename StatsVector>
        typename std::enable_if<
          !std::is_same<typename Pattern::SubPattern,void>::value &&
           requires_pattern<Container>::value
        >::type
        allocate_bcrs_matrix(const OrderingV& ordering_v,
                             const OrderingU& ordering_u,
                             Pattern& p,
                             Container& c,
                             StatsVector& stats)
        {
          c.setSize(ordering_v.blockCount(),ordering_u.blockCount(),ordering_v.blockCount()*ordering_u.blockCount());
          c.setBuildMode(Container::random);

          for (std::size_t i = 0; i < c.N(); ++i)
            c.setrowsize(i,ordering_u.blockCount());
          c.endrowsizes();

          for (std::size_t i = 0; i < c.N(); ++i)
            for (std::size_t j = 0; j < c.M(); ++j)
              c.addindex(i,j);
          c.endindices();

          for (std::size_t i = 0; i < c.N(); ++i)
            for (std::size_t j = 0; j < c.M(); ++j)
              {
                allocate_bcrs_matrix(ordering_v.childOrdering(i),
                                     ordering_u.childOrdering(j),
                                     p.subPattern(i,j),
                                     c[i][j],
                                     stats);
              }
        }

      } // anonymous namespace



      //! Backend using (possibly nested) ISTL BCRSMatrices.
      /**
       * BCRSMatrixBackend is a matrix backend descriptor for ISTL matrices. It expects that
       * both the ansatz and the test function space use ISTL vectors and automatically deduces
       * the correct matrix type from those two vector backends.
       *
       * The backend uses an accelerated pattern construction scheme, which requires the average
       * number of non-zero entries per matrix row as a priori information. In constrast to the
       * older construction scheme, the improved version never requires more memory than the matrix
       * does after pattern construction and runs a lot faster, as long as it is provided with a
       * reasonable estimate for the number of non-zero entries per row.
       *
       */
      template<typename EntriesPerRow = std::size_t>
      struct BCRSMatrixBackend
      {

        //! The size type of the BCRSMatrix.
        typedef std::size_t size_type;

        //! The type of the object holding the statistics generated during pattern construction.
        typedef PatternStatistics<size_type> Statistics;

        //! The type of the pattern object passed to the GridOperator for pattern construction.
        template<typename Matrix, typename GFSV, typename GFSU>
        using Pattern = typename build_bcrs_pattern_type<
          typename Matrix::Container,
          GFSV,
          GFSU,
          typename GFSV::Ordering::ContainerAllocationTag
          >::type;

        template<typename VV, typename VU, typename E>
        struct MatrixHelper
        {
          typedef BCRSMatrix<
            typename VV::GridFunctionSpace,
            typename VU::GridFunctionSpace,
            typename build_matrix_type<
              E,
              typename VV::Container,
              typename VU::Container
              >::type,
            Statistics
            > type;
        };

        //! Builds the matrix pattern associated with grid_operator and initializes matrix with it.
        /**
         * \returns  a vector with statistics object for all leaf BCRSMatrices in row-major order.
         */
        template<typename GridOperator, typename Matrix>
        std::vector<Statistics> buildPattern(const GridOperator& grid_operator, Matrix& matrix) const
        {
          Pattern<
            Matrix,
            typename GridOperator::Traits::TestGridFunctionSpace,
            typename GridOperator::Traits::TrialGridFunctionSpace
            > pattern(grid_operator.testGridFunctionSpace().ordering(),grid_operator.trialGridFunctionSpace().ordering(),_entries_per_row);
          grid_operator.fill_pattern(pattern);
          std::vector<Statistics> stats;
          allocate_bcrs_matrix(grid_operator.testGridFunctionSpace().ordering(),
                               grid_operator.trialGridFunctionSpace().ordering(),
                               pattern,
                               Backend::native(matrix),
                               stats
                               );
          return stats;
        }

        //! Constructs a BCRSMatrixBackend.
        /**
         * TODO: Document and flesh out the way this should work for nested matrices (use a nested array as entries_per_row).
         *
         * \param entries_per_row  The average number of nonzero entries per row in matrices created with this backend.
         */
        BCRSMatrixBackend(const EntriesPerRow& entries_per_row)
          : _entries_per_row(entries_per_row)
        {}

      private:

        EntriesPerRow _entries_per_row;

      };

    } // namespace ISTL
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_BCRSMATRIXBACKEND_HH
