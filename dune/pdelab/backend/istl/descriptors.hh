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

#ifndef DOXYGEN
      template<typename T>
      constexpr bool deactivate_standard_blocking_for_ordering(const T&)
      {
        return false;
      }
#endif

    namespace istl {

      //! The type of blocking employed at this node in the function space tree.
      enum class Blocking
      {
        //! No blocking at this level.
        none,
        //! Creates one macro block for each child space, each block is a BlockVector / BCRS matrix.
        bcrs,
        //! Create fixed size blocks that each group together a fixed number of DOFs from each child space.
        /**
         * Creates a block structure with fixed size child blocks that each group together a fixed number
         * of DOFs from each child space. Typically used for entity-wise blocking of DOFs across child spaces
         * for e.g. velocity in flow problems or concentrations in multi-component transport.
         *
         * \note This type of blocking cannot be nested due to limitations in ISTL.
         */
        fixed,
      };

      //! Tag describing an ISTL BlockVector backend.
      struct vector_backend_tag {};

      template<Blocking blocking = Blocking::none, std::size_t block_size_ = 1>
      struct VectorBackend
      {

        using tag = vector_backend_tag;

        static_assert((block_size_ > 0),"block size for FieldVector has to be positive");

        using size_type = std::size_t;

        static const size_type blockSize = block_size_;

        struct Traits
        {
          static const Blocking block_type = blocking;
          static const size_type block_size = block_size_;

          static const bool blocked = blocking != Blocking::none;

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
          return Traits::blocked && (blocking != Blocking::fixed || !GFS::isLeaf || block_size_ > 1);
        }

      };

    }

    namespace ISTLParameters {

      enum
      DUNE_DEPRECATED_MSG("ISTLParameters::blocking is deprecated and will be removed after PDELab 2.4. Use the new istl::VectorBackend and istl::Blocking instead. Note that the enum values of istl::Blocking are named differently!")
      Blocking
        {
          no_blocking,
          dynamic_blocking,
          static_blocking
        };
    }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

    template<ISTLParameters::Blocking blocking = ISTLParameters::no_blocking, std::size_t block_size = 1>
    using ISTLVectorBackend DUNE_DEPRECATED_MSG("ISTLVectorBackend is deprecated and will be removed after PDELab 2.4. Use istl::VectorBackend instead") = istl::VectorBackend<static_cast<istl::Blocking>(blocking),block_size>;

#pragma GCC diagnostic pop

    //! Backend using ISTL matrices.
    /**
     * ISTLMatrixBackend is a matrix backend descriptor for ISTL matrices. It expects that
     * both the ansatz and the test function space use ISTL vectors and automatically deduces
     * the correct matrix type from those two vector backends.
     */
    struct
    DUNE_DEPRECATED_MSG("ISTLMatrixBackend has been deprecated and will be removed after the release of PDELab 2.4. Use istl::BCRSMatrixBackend with the newer pattern construction method instead")
    ISTLMatrixBackend
    {

      typedef std::size_t size_type;

      // The default matrix construction method does not collect statistics, so provide a dummy type here.
      typedef int Statistics;

      //! The type of the pattern object passed to the GridOperator for pattern construction.
      template<typename Matrix, typename GFSV, typename GFSU>
      using Pattern = typename istl::build_pattern_type<
        typename Matrix::Container,
        GFSV,
        GFSU,
        typename GFSV::Ordering::ContainerAllocationTag
        >::type;

      template<typename VV, typename VU, typename E>
      struct MatrixHelper
      {
        typedef istl::BCRSMatrix<
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
