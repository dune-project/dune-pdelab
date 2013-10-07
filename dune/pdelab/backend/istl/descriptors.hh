// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH
#define DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH

#include <dune/common/static_assert.hh>
#include <dune/pdelab/backend/istl/forwarddeclarations.hh>
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

    struct istl_vector_backend_tag {};

    template<ISTLParameters::Blocking blocking = ISTLParameters::no_blocking, std::size_t block_size_ = 1>
    struct ISTLVectorBackend
    {

      typedef istl_vector_backend_tag tag;

      dune_static_assert((block_size_ > 0),"block size for FieldVector has to be positive");

      typedef std::size_t size_type;

      static const size_type blockSize = block_size_;

      struct Traits
      {
        static const ISTLParameters::Blocking block_type = blocking;
        static const size_type block_size = block_size_;
        static const bool blocked = blocking != ISTLParameters::no_blocking;
        static const size_type max_blocking_depth = blocked ? 1 : 0;
      };

      bool blocked() const
      {
        return Traits::blocked;
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

      template<typename VV, typename VU, typename E>
      struct MatrixHelper
      {
        typedef ISTLMatrixContainer<typename VV::GridFunctionSpace,typename VU::GridFunctionSpace,typename istl::build_matrix_type<E,typename VV::Container,typename VU::Container>::type > type;
      };
    };

    namespace istl {

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

        bool blocked() const
        {
          return false;
        }

      };

      template<typename Alloc>
      struct FlatMatrixBackend
      {

        typedef Alloc Allocator;
        typedef Alloc allocator_type;

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

        bool blocked() const
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

        template<typename VV, typename VU, typename E>
        struct MatrixHelper
        {
          typedef BELLMatrixContainer<
            typename VV::GridFunctionSpace,
            typename VU::GridFunctionSpace,
            Dune::ISTL::BELLMatrix<E,Alloc>
            > type;
        };
      };

    }

  } // namespace PDELab
} // namespace Dune



#endif // DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH
