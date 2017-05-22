#ifndef DUNE_PDELAB_BACKEND_INTERFACE_HH
#define DUNE_PDELAB_BACKEND_INTERFACE_HH

#include <type_traits>
#include <utility>

namespace Dune {
  namespace PDELab {
    namespace Backend {

#ifndef DOXYGEN

      namespace impl {

        // This class needs to be specialized by each backend and return the correct
        // vector wrapper type via the nested type named "type"
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

        // marker mixin type used by the generic implementation below to decide whether a type
        // is a PDELab wrapper around a native object - for internal use only!
        struct WrapperBase
        {};

        // mixin class for PDELab wrappers around vectors and matrices.
        // All backend wrappers should inherit from this class and set the
        // template parameter to the type of the native container of the backend that they are
        // wrapping.
        // Moreover, they have to provide methods
        //
        // NativeContainer& native()
        // const NativeContainer& native() const
        //
        // that provide access to the wrapped data structure. These two methods should be private;
        // in that case, the Wrapper<NativeContainer> mixin must be a friend of the wrapper implementation
        // so that it can access the native() methods.
        template<typename NativeContainer>
        struct Wrapper
          : public WrapperBase
        {

          using native_type = NativeContainer;

          template<typename U>
          static auto access_native(U&& u) -> decltype(u.native())
          {
            // make sure the wrapper actually returns the right type of object
            static_assert(
              std::is_same<
                typename std::decay<
                  decltype(u.native())
                  >::type,
                NativeContainer
                >::value,
              "u.native() must return a cv-qualified xvalue of type T"
              );

            return u.native();
          }

        };

      } // namespace impl

    } // namespace Backend


    namespace Backend{

#endif // DOXYGEN


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
      using Vector = typename impl::BackendVectorSelector<GridFunctionSpace, FieldType>::Type;

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
      using Matrix = typename impl::BackendMatrixSelector<Backend, VU, VV, E>::Type;


      namespace {

        // helper TMP to deduce the native type of a wrapper. It works by checking whether
        // T inherits from WrapperBase. In that case, it returns the nested typedef native_type,
        // otherwise T is returned unchanged.
        template<typename T>
        struct native_type
        {

          // We need to defer the (possible) extraction of the nested type until we are sure
          // that the nested type actually exists. That'w what these lazy_... structs are for:
          // a std::conditional picks the correct version, and the actual evaluation only happens
          // after the evaluation of the std::conditional.

          struct lazy_identity
          {
            template<typename U>
            struct evaluate
            {
              using type = U;
            };

          };

          struct lazy_native_type
          {

            template<typename U>
            struct evaluate
            {
              using type = typename U::native_type;
            };

          };

          using type = typename std::conditional<
            std::is_base_of<impl::WrapperBase,T>::value,
            lazy_native_type,
            lazy_identity
            >::type::template evaluate<T>::type;
        };

      }

      //! Alias of the native container type associated with T or T itself if it is not a backend wrapper.
      template<typename T>
      using Native = typename native_type<T>::type;

#ifdef DOXYEN

      //! Returns the native container embedded in t or t itself if it is not a container.
      template<typename T>
      Native<T> native(T&& t);

#else // DOXYGEN

      // version for mutable reference to wrapper
      template<typename T>
      typename std::enable_if<
        std::is_base_of<impl::WrapperBase,T>::value,
        Native<T>&
        >::type
      native(T& t)
      {
        return impl::Wrapper<Native<T>>::access_native(t);
      }

      // version for const reference to wrapper
      template<typename T>
      typename std::enable_if<
        std::is_base_of<impl::WrapperBase,T>::value,
        const Native<T>&
        >::type
      native(const T& t)
      {
        return impl::Wrapper<Native<T>>::access_native(t);
      }

      // version for non-wrapper class. Important: Don't drop the std::decay<>!
      template<typename T>
      typename std::enable_if<
        !std::is_base_of<impl::WrapperBase,typename std::decay<T>::type>::value,
        decltype(std::forward<T>(std::declval<T&&>()))
        >::type
      native(T&& t)
      {
        return std::forward<T>(t);
      }

#endif // DOXYGEN

    } // namespace Backend
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_INTERFACE_HH
