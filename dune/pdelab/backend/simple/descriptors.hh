// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_SIMPLE_DESCRIPTORS_HH
#define DUNE_PDELAB_BACKEND_SIMPLE_DESCRIPTORS_HH

#include <vector>

#include <dune/pdelab/ordering/orderingbase.hh>

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    namespace simple {

      template<typename GFS, typename C>
      class VectorContainer;

      template<typename GFSV, typename GFSU, typename C>
      class MatrixContainer;

      template<typename GFSV, typename GFSU, template<typename> class C, typename ET, typename I>
      class SparseMatrixContainer;

      template<typename _RowOrdering, typename _ColOrdering>
      class SparseMatrixPattern;

      template<typename E>
      using default_vector = std::vector<E>;

    }

#endif // DOXYGEN

    template<template<typename> class Container = simple::default_vector>
    struct SimpleVectorBackend
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

    template<template<typename> class Container = simple::default_vector>
    struct SimpleMatrixBackend
    {

      typedef std::size_t size_type;

      template<typename Matrix, typename GFSV, typename GFSU>
      struct Pattern
      {};

      template<typename VV, typename VU, typename E>
      struct MatrixHelper
      {
        typedef simple::MatrixContainer<typename VV::GridFunctionSpace,typename VU::GridFunctionSpace,Container<E> > type;
      };
    };

    template<template<typename> class Container = simple::default_vector, typename IndexType = std::size_t>
    struct SimpleSparseMatrixBackend
    {

      typedef IndexType size_type;

#if HAVE_TEMPLATE_ALIASES || DOXYGEN

      //! The type of the pattern object passed to the GridOperator for pattern construction.
      template<typename Matrix, typename GFSV, typename GFSU>
      using Pattern = simple::SparseMatrixPattern<
        OrderingBase<
          typename GFSV::Ordering::Traits::DOFIndex,
          typename GFSV::Ordering::Traits::ContainerIndex
          >,
        OrderingBase<
          typename GFSU::Ordering::Traits::DOFIndex,
          typename GFSU::Ordering::Traits::ContainerIndex> >;

#else // HAVE_TEMPLATE_ALIASES

      template<typename Matrix, typename GFSV, typename GFSU>
      struct Pattern
        : public simple::SparseMatrixPattern<
        OrderingBase<
          typename GFSV::Ordering::Traits::DOFIndex,
          typename GFSV::Ordering::Traits::ContainerIndex
          >,
        OrderingBase<
          typename GFSU::Ordering::Traits::DOFIndex,
          typename GFSU::Ordering::Traits::ContainerIndex> >
      {

        typedef SparseMatrixPattern<
          typedef OrderingBase<
            typename GFSV::Ordering::Traits::DOFIndex,
            typename GFSV::Ordering::Traits::ContainerIndex
            >,
          OrderingBase<
            typename GFSU::Ordering::Traits::DOFIndex,
            typename GFSU::Ordering::Traits::ContainerIndex> > BaseT;

        typedef BaseT::RowOrdering RowOrdering;
        typedef BaseT::ColOrdering ColOrdering;

        Pattern(const RowOrdering& row_ordering, const ColOrdering& col_ordering)
          : BaseT(row_ordering,col_ordering)
        {}

      };

#endif // HAVE_TEMPLATE_ALIASES

      template<typename VV, typename VU, typename E>
      struct MatrixHelper
      {
        typedef simple::SparseMatrixContainer<typename VV::GridFunctionSpace,typename VU::GridFunctionSpace,Container, E, size_type> type;
      };
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_SIMPLE_DESCRIPTORS_HH
