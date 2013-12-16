// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_DENSE_DESCRIPTORS_HH
#define DUNE_PDELAB_BACKEND_DENSE_DESCRIPTORS_HH

#include <vector>

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    namespace dense {

      template<typename GFS, typename C>
      class VectorContainer;

      template<typename GFSV, typename GFSU, typename C>
      class MatrixContainer;

      template<typename E>
      using default_vector = std::vector<E>;

    }

#endif // DOXYGEN

    template<template<typename> class Container = dense::default_vector>
    struct DenseVectorBackend
    {
      template<typename E>
      using vector_type = Container<E>;

      typedef typename vector_type<double>::size_type size_type;

      struct Traits
      {
        static const size_type max_blocking_depth = 0;
      };

      bool blocked() const
      {
        return false;
      }

    };

    template<template<typename> class Container = dense::default_vector>
    struct DenseMatrixBackend
    {

      typedef std::size_t size_type;

      template<typename Matrix, typename GFSV, typename GFSU>
      struct Pattern
      {};

      template<typename VV, typename VU, typename E>
      struct MatrixHelper
      {
        typedef dense::MatrixContainer<typename VV::GridFunctionSpace,typename VU::GridFunctionSpace,Container<E> > type;
      };
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_DENSE_DESCRIPTORS_HH
