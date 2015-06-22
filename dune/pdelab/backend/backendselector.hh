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
    }
  }
}

#endif
