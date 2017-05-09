// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH
#define DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <dune/pdelab/backend/interface.hh>
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

    namespace ISTL {

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

#ifndef DOXYGEN

      namespace flat {

        // forward declaration for backend
        template<typename GFSV, typename GFSU>
        class Pattern;

      }

#endif // DOXYGEN

      template<typename Alloc>
      struct FlatVectorBackend
      {

        typedef Alloc Allocator;
        typedef Alloc allocator_type;

        typedef typename Alloc::size_type size_type;

        static const size_type blockSize = 1;

        struct Traits
        {
          static const size_type block_size = 1;
          static const bool blocked = false;
          static const size_type max_blocking_depth = 0;
        };

        template<typename GFS>
        bool blocked(const GFS& gfs) const
        {
          return false;
        }

      };

      template<typename Alloc>
      struct FlatMatrixBackend
      {

        typedef Alloc Allocator;
        typedef Alloc allocator_type;

        // The ELL matrix construction process does not collect statistics, so provide a dummy type here.
        typedef int Statistics;

        //! The type of the pattern object passed to the GridOperator for pattern construction.
        template<typename Matrix, typename GFSV, typename GFSU>
        using Pattern = flat::Pattern<typename GFSV::Ordering,typename GFSU::Ordering>;

        typedef typename Alloc::size_type size_type;

        template<typename VV, typename VU, typename E>
        struct MatrixHelper
        {
          typedef FlatELLMatrixContainer<
            typename VV::GridFunctionSpace,
            typename VU::GridFunctionSpace,
            Dune::ISTL::ELLMatrix<E,Alloc>
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

        std::size_t entriesPerRow() const
        {
          return _entries_per_row;
        }

        FlatMatrixBackend(std::size_t entries_per_row)
          : _entries_per_row(entries_per_row)
        {}

      private:

        std::size_t _entries_per_row;

      };



      template<typename Alloc>
      struct BlockVectorBackend
      {

        typedef Alloc Allocator;
        typedef Alloc allocator_type;

        typedef typename Alloc::size_type size_type;

        struct Traits
        {
          static const size_type max_blocking_depth = 1;
        };

        template<typename GFS>
        bool blocked(const GFS& gfs) const
        {
          return true;
        }

        size_type blockSize() const
        {
          return _block_size;
        }

        BlockVectorBackend(size_type block_size)
          : _block_size(block_size)
        {}

      private:

        const size_type _block_size;

      };

      template<typename Alloc>
      struct BELLMatrixBackend
      {

        typedef Alloc Allocator;
        typedef Alloc allocator_type;

        typedef typename Alloc::size_type size_type;

      // The BELL matrix construction process does not collect statistics, so provide a dummy type here.
      typedef int Statistics;

        //! The type of the pattern object passed to the GridOperator for pattern construction.
        template<typename Matrix, typename GFSV, typename GFSU>
        using Pattern = flat::Pattern<typename GFSV::Ordering,typename GFSU::Ordering>;

        template<typename VV, typename VU, typename E>
        struct MatrixHelper
        {
          typedef BELLMatrixContainer<
            typename VV::GridFunctionSpace,
            typename VU::GridFunctionSpace,
            Dune::ISTL::BELLMatrix<E,Alloc>
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

    }

  } // namespace PDELab
} // namespace Dune


#endif // DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH
