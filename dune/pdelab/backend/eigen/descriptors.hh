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

    namespace EIGEN {

      template<typename GFS, typename E>
      class VectorContainer;

      template<typename GFSV, typename GFSU, typename ET, int _Options>
      class MatrixContainer;

      template<typename M>
      struct MatrixPatternInserter;

    }

#endif // DOXYGEN

    /** Eigen backend for FunctionSpace
     *
     * works with EIGEN >= 3.1
     */
    struct EigenVectorBackend
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

    template<int _Options = Eigen::RowMajor>
    struct EigenMatrixBackend
    {

      typedef std::size_t size_type;

#if HAVE_TEMPLATE_ALIASES || DOXYGEN

      //! The type of the pattern object passed to the GridOperator for pattern construction.
      template<typename Matrix, typename GFSV, typename GFSU>
      using Pattern = EIGEN::MatrixPatternInserter<typename Matrix::Container>;

#else // HAVE_TEMPLATE_ALIASES

      template<typename Matrix, typename GFSV, typename GFSU>
      struct Pattern
        : public EIGEN::MatrixPatternInserter<typename Matrix::Container>
      {

        typedef EIGEN::MatrixPatternInserter<typename Matrix::Container> BaseT;

        Pattern(Matrix & m)
          : BaseT(m)
        {}

      };

#endif // HAVE_TEMPLATE_ALIASES

      template<typename VV, typename VU, typename E>
      struct MatrixHelper
      {
        typedef EIGEN::MatrixContainer<typename VV::GridFunctionSpace,typename VU::GridFunctionSpace, E, _Options> type;
      };
    };

  } // namespace PDELab
} // namespace Dune

#elif defined HEADERCHECK
#warning Skipped header check due to missing Eigen.
#else
#error You need Eigen to use the Eigen backend
#endif // HAVE_EIGEN

#endif // DUNE_PDELAB_BACKEND_EIGEN_DESCRIPTORS_HH
