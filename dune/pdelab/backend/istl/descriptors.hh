// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH
#define DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH

#include <dune/pdelab/backend/istl/forwarddeclarations.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/utility.hh>
#include <cstddef>

namespace Dune {
  namespace PDELab {

    namespace ISTLParameters {

      enum Blocking
        {
          no_blocking,
          dynamic_blocking,
          static_blocking
        };
    }

    template<typename T>
    DUNE_CONSTEXPR bool deactivate_standard_blocking_for_ordering(const T&)
    {
      return false;
    }

    struct istl_vector_backend_tag {};

    template<ISTLParameters::Blocking blocking = ISTLParameters::no_blocking, std::size_t block_size_ = 1>
    struct ISTLVectorBackend
    {

      typedef istl_vector_backend_tag tag;

      static_assert((block_size_ > 0),"block size for FieldVector has to be positive");

      typedef std::size_t size_type;

      static const size_type blockSize = block_size_;

      struct Traits
      {
        static const ISTLParameters::Blocking block_type = blocking;
        static const size_type block_size = block_size_;

        static const bool blocked = blocking != ISTLParameters::no_blocking;

        static const size_type max_blocking_depth = blocked ? 1 : 0;
      };

      template<typename GFS>
      bool blocked(const GFS& gfs) const
      {
        if (deactivate_standard_blocking_for_ordering(gfs.orderingTag()))
          return false;
        // We have to make an exception for static blocking and block_size == 1:
        // In that case, the ISTL backends expect the redundant index information
        // at that level to be elided, and keeping it in here will break vector
        // and matrix accesses.
        // To work around that problem, we override the user and just turn off
        // blocking internally.
        return Traits::blocked && (blocking != ISTLParameters::static_blocking || !GFS::isLeaf || block_size_ > 1);
      }

    };

    //! Backend using ISTL matrices.
    /**
     * ISTLMatrixBackend is a matrix backend descriptor for ISTL matrices. It expects that
     * both the ansatz and the test function space use ISTL vectors and automatically deduces
     * the correct matrix type from those two vector backends.
     */
    struct ISTLMatrixBackend
    {

      typedef std::size_t size_type;

      // The default matrix construction method does not collect statistics, so provide a dummy type here.
      typedef int Statistics;

#if HAVE_TEMPLATE_ALIASES || DOXYGEN

      //! The type of the pattern object passed to the GridOperator for pattern construction.
      template<typename Matrix, typename GFSV, typename GFSU>
      using Pattern = typename istl::build_pattern_type<
        typename Matrix::Container,
        GFSV,
        GFSU,
        typename GFSV::Ordering::ContainerAllocationTag
        >::type;

#else // HAVE_TEMPLATE_ALIASES

      template<typename Matrix, typename GFSV, typename GFSU>
      struct Pattern
        : public istl::build_pattern_type<typename Matrix::Container,
                                          GFSV,
                                          GFSU,
                                          typename GFSV::Ordering::ContainerAllocationTag
                                          >::type
      {

        typedef OrderingBase<
          typename GFSV::Ordering::Traits::DOFIndex,
          typename GFSV::Ordering::Traits::ContainerIndex
          > RowOrdering;

        typedef OrderingBase<
          typename GFSU::Ordering::Traits::DOFIndex,
          typename GFSU::Ordering::Traits::ContainerIndex
          > ColOrdering;

        typedef typename istl::build_pattern_type<
          typename Matrix::Container,
          GFSV,
          GFSU,
          typename GFSV::Ordering::ContainerAllocationTag
          >::type BaseT;

        Pattern(const RowOrdering& row_ordering, const ColOrdering& col_ordering)
          : BaseT(row_ordering,col_ordering)
        {}

      };

#endif // HAVE_TEMPLATE_ALIASES

      template<typename VV, typename VU, typename E>
      struct MatrixHelper
      {
        typedef ISTLMatrixContainer<
          typename VV::GridFunctionSpace,
          typename VU::GridFunctionSpace,
          typename istl::build_matrix_type<
            E,
            typename VV::Container,
            typename VU::Container
            >::type,
          Statistics
          > type;
      };

      template<typename GridOperator, typename Matrix>
      std::vector<Statistics> buildPattern(const GridOperator& grid_operator, Matrix& matrix) const
      {
        Pattern<
          Matrix,
          typename GridOperator::Traits::TestGridFunctionSpace,
          typename GridOperator::Traits::TrialGridFunctionSpace
          > pattern(grid_operator.testGridFunctionSpace().ordering(),grid_operator.trialGridFunctionSpace().ordering());
        grid_operator.fill_pattern(pattern);
        allocate_matrix(grid_operator.testGridFunctionSpace().ordering(),
                        grid_operator.trialGridFunctionSpace().ordering(),
                        pattern,
                        Backend::native(matrix)
                        );
        return std::vector<Statistics>();
      }

    };

  } // namespace PDELab
} // namespace Dune



#endif // DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH
