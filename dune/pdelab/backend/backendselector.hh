#ifndef DUNE_PDELAB_BACKENDSELECTOR_HH
#define DUNE_PDELAB_BACKENDSELECTOR_HH
namespace Dune {
  namespace PDELab {

    template<typename Backend, typename GridFunctionSpace, typename FieldType>
    struct BackendVectorSelectorHelper
    {};

    template<typename GridFunctionSpace, typename FieldType>
    struct BackendVectorSelector
    {
      typedef typename GridFunctionSpace::Traits::Backend Backend;
      typedef typename BackendVectorSelectorHelper<Backend, GridFunctionSpace, FieldType>::Type Type;
    };

    template<typename Backend, typename VU, typename VV, typename E>
    struct BackendMatrixSelector
    {
      typedef typename Backend::template MatrixHelper<VV,VU,E>::type Type;
    };

    namespace Backend
    {
      /**
       * \brief alias of the return type of BackendVectorSelector
       *
       * This alias can be used as a short hand for retrieving the vector type for
       * a grid function space and a given field type. The typedef
       * \code
       * typedef typename Dune::PDELab::BackendVectorSelector<GFS,FT>::Type Vec;
       * \endcode
       * simplifies to
       * \code
       * typedef Dune::PDELab::Backend::Vector<GFS,FT> Vec;
       * \endcode
       * or
       * \code
       * using Vec = Dune::PDELab::Backend::Vector<GFS,FT>;
       * \endcode
       *
       **/
      template<typename GridFunctionSpace, typename FieldType>
      using Vector = typename BackendVectorSelector<GridFunctionSpace, FieldType>::Type;

      /**
       * \brief alias of the return type of BackendMatrixSelector
       *
       * This alias can be used as a short hand for retrieving the matrix type for
       * given matrix backend, domain, range and field type.
       * \code
       * typedef typename Dune::PDELab::BackendMatrixSelector<Backend,VU,VV,E>::Type Mat;
       * \endcode
       * simplifies to
       * \code
       * typedef Dune::PDELab::Backend::Matrix<Backend,VU,VV,E> Mat;
       * \endcode
       * or
       * \code
       * using Mat = Dune::PDELab::Backend::Matrix<Backend,VU,VV,E>;
       * \endcode
       *
       **/
      template<typename Backend, typename VU, typename VV, typename E>
      using Matrix = typename BackendMatrixSelector<Backend, VU, VV, E>::Type;
    }
  }
}

#endif
