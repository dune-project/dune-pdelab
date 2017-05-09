// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_SIMPLE_DESCRIPTORS_HH
#define DUNE_PDELAB_BACKEND_SIMPLE_DESCRIPTORS_HH

#include <vector>

#include <dune/pdelab/ordering/orderingbase.hh>

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    namespace Simple {

      template<typename GFS, typename C>
      class VectorContainer;

      template<typename GFSV, typename GFSU, typename C>
      class MatrixContainer;

      template<typename GFSV, typename GFSU, template<typename> class C, typename ET, typename I>
      class SparseMatrixContainer;

      class SparseMatrixPattern;

      template<typename E>
      using default_vector = std::vector<E>;

    }

#endif // DOXYGEN

    namespace Simple {

      template<template<typename> class Container = Simple::default_vector>
      struct VectorBackend
      {
        template<typename E>
        using vector_type = Container<E>;

        typedef typename vector_type<double>::size_type size_type;

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

      template<template<typename> class Container = Simple::default_vector>
      struct MatrixBackend
      {

        typedef std::size_t size_type;

        template<typename Matrix, typename GFSV, typename GFSU>
        struct Pattern
        {};

        template<typename VV, typename VU, typename E>
        struct MatrixHelper
        {
          typedef Simple::MatrixContainer<typename VV::GridFunctionSpace,typename VU::GridFunctionSpace,Container<E> > type;
        };
      };

      template<template<typename> class Container = Simple::default_vector, typename IndexType = std::size_t>
      struct SparseMatrixBackend
      {

        typedef IndexType size_type;

        //! The type of the pattern object passed to the GridOperator for pattern construction.
        template<typename Matrix, typename GFSV, typename GFSU>
        using Pattern = Simple::SparseMatrixPattern;

        template<typename VV, typename VU, typename E>
        struct MatrixHelper
        {
          typedef Simple::SparseMatrixContainer<typename VV::GridFunctionSpace,typename VU::GridFunctionSpace,Container, E, size_type> type;
        };
      };

    } // namespace Simple

    /** \brief For backward compatibility: access to Simple::VectorBackend by its old name
     * \deprecated Use Simple::VectorBackend instead!
     */
    template<template<typename> class Container = Simple::default_vector>
    using SimpleVectorBackend DUNE_DEPRECATED_MSG("Use Simple::VectorBackend instead of SimpleVectorBackend")
        = Simple::VectorBackend<Container>;

    /** \brief For backward compatibility: access to Simple::MatrixBackend by its old name
     * \deprecated Use Simple::MatrixBackend instead!
     */
    template<template<typename> class Container = Simple::default_vector>
    using SimpleMatrixBackend DUNE_DEPRECATED_MSG("Use Simple::MatrixBackend instead of SimpleMatrixBackend")
        = Simple::MatrixBackend<Container>;

    /** \brief For backward compatibility: access to Simple::SparseMatrixBackend by its old name
     * \deprecated Use Simple::SparseMatrixBackend instead!
     */
    template<template<typename> class Container = Simple::default_vector, typename IndexType = std::size_t>
    using SimpleSparseMatrixBackend DUNE_DEPRECATED_MSG("Use Simple::SparseMatrixBackend instead of SimpleSparseMatrixBackend")
        = Simple::SparseMatrixBackend<Container, IndexType>;

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_SIMPLE_DESCRIPTORS_HH
