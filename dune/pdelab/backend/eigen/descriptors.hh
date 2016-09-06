// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_EIGEN_DESCRIPTORS_HH
#define DUNE_PDELAB_BACKEND_EIGEN_DESCRIPTORS_HH

#include <vector>

#if HAVE_EIGEN

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    namespace Eigen {

      template<typename GFS, typename E>
      class VectorContainer;

      template<typename GFSV, typename GFSU, typename ET, int _Options>
      class MatrixContainer;

      template<typename M>
      struct MatrixPatternInserter;

    }

#endif // DOXYGEN

    namespace Eigen {

      /** Eigen backend for FunctionSpace
       *
       * works with EIGEN >= 3.1
       */
      struct VectorBackend
      {
        typedef std::size_t size_type;

        struct Traits
        {
          static const size_type max_blocking_depth = 0;
        };

        template<typename GFS>
        bool blocked(const GFS& gfs) const
        {
          return false;
        }

      };

      template<int _Options = ::Eigen::RowMajor>
      struct MatrixBackend
      {
        typedef std::size_t size_type;

        size_type avg_nz_per_row;

        DUNE_DEPRECATED_MSG("Please us the constructor taking the avg non-zeros")
        MatrixBackend() : avg_nz_per_row(0)
        {}

        MatrixBackend(size_type avg_nz_per_row_) : avg_nz_per_row(avg_nz_per_row_)
        {}

        //! The type of the pattern object passed to the GridOperator for pattern construction.
        template<typename Matrix, typename GFSV, typename GFSU>
        using Pattern = PDELab::Eigen::MatrixPatternInserter<typename Matrix::Container>;

        template<typename VV, typename VU, typename E>
        struct MatrixHelper
        {
          typedef PDELab::Eigen::MatrixContainer<typename VV::GridFunctionSpace,typename VU::GridFunctionSpace, E, _Options> type;
        };
      };

    } // namespace Eigen

    /** \brief For backward compatibility: access to Eigen::VectorBackend by its old name
     * \deprecated Use Eigen::VectorBackend instead!
     */
    using EigenVectorBackend DUNE_DEPRECATED_MSG("Use Eigen::VectorBackend instead of EigenVectorBackend")
        = Eigen::VectorBackend;

    /** \brief For backward compatibility: access to Eigen::MatrixBackend by its old name
     * \deprecated Use Eigen::MatrixBackend instead!
     */
    template<int _Options = ::Eigen::RowMajor>
    using EigenMatrixBackend DUNE_DEPRECATED_MSG("Use Eigen::MatrixBackend instead of EigenMatrixBackend")
        = Eigen::MatrixBackend<_Options>;

  } // namespace PDELab
} // namespace Dune

#elif defined HEADERCHECK
#warning Skipped header check due to missing Eigen.
#else
#error You need Eigen to use the Eigen backend
#endif // HAVE_EIGEN

#endif // DUNE_PDELAB_BACKEND_EIGEN_DESCRIPTORS_HH
