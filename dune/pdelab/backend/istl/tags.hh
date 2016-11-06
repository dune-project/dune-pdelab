// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_TAGS_HH
#define DUNE_PDELAB_BACKEND_ISTL_TAGS_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <dune/common/documentation.hh>
#include <dune/pdelab/backend/istl/forwarddeclarations.hh>
#include <cstddef>

namespace Dune {

  namespace PDELab {

    namespace ISTL {

      // ********************************************************************************
      // tag definitions
      // ********************************************************************************

      namespace tags {

        //! Tag describing a BlockVector.
        struct block_vector
        {
          typedef block_vector base_tag;
        };

        //! Tag describing a DynamicVector.
        struct dynamic_vector
        {
          typedef dynamic_vector base_tag;
        };

        //! Tag describing an arbitrary FieldVector.
        /**
         * The tagging mechanism will never tag a FieldVector directly with this tag,
         * but will instead return a more specific tag which describes the blocking structure
         * and which inherits from this tag.
         * If you dispatch using a function, you can ignore this detail, as the derived tag will
         * automatically be upcast to this tag, but if you use the tags in a struct, you can
         * always use the member type base_tag to access the common base tag.
         */
        struct field_vector
        {
          //! Base tag for this tag category.
          typedef field_vector base_tag;
        };

        //! Tag describing a field vector with block size 1.
        struct field_vector_1
          : public field_vector
        {};

        //! Tag describing a field vector with block size > 1.
        struct field_vector_n
          : public field_vector
        {};

        //! Tag describing a BCRSMatrix.
        struct bcrs_matrix
        {
          typedef bcrs_matrix base_tag;
        };

        //! Tag describing a DynamicMatrix.
        struct dynamic_matrix
        {
          typedef dynamic_matrix base_tag;
        };

        //! Tag describing an arbitrary FieldMatrix.
        /**
         * The tagging mechanism will never tag a FieldMatrix directly with this tag,
         * but will instead return a more specific tag which describes the blocking structure
         * and which inherits from this tag.
         * If you dispatch using a function, you can ignore this detail, as the derived tag will
         * automatically be upcast to this tag, but if you use the tags in a struct, you can
         * always use the member type base_tag to access the common base tag.
         */
        struct field_matrix
        {
          //! Base tag for this tag category.
          typedef field_matrix base_tag;
        };

        //! Tag describing a FieldMatrix with row block size 1 and arbitrary column block size.
        struct field_matrix_1_any
        {};

        //! Tag describing a FieldMatrix with row block size > 1 and arbitrary column block size.
        struct field_matrix_n_any
        {};

        //! Tag describing a FieldMatrix with arbitrary row block size and column block size 1.
        struct field_matrix_any_1
        {};

        //! Tag describing a FieldMatrix with arbitrary row block size and column block size > 1.
        struct field_matrix_any_m
        {};

        //! Tag describing a FieldMatrix with row block size 1 and column block size 1.
        struct field_matrix_1_1
          : public field_matrix
          , public field_matrix_1_any
          , public field_matrix_any_1
        {};

        //! Tag describing a FieldMatrix with row block size > 1 and column block size 1.
        struct field_matrix_n_1
          : public field_matrix
          , public field_matrix_n_any
          , public field_matrix_any_1
        {};

        //! Tag describing a FieldMatrix with row block size 1 and column block size > 1.
        struct field_matrix_1_m
          : public field_matrix
          , public field_matrix_1_any
          , public field_matrix_any_m
        {};

        //! Tag describing a FieldMatrix with row block size > 1 and column block size > 1.
        struct field_matrix_n_m
          : public field_matrix
          , public field_matrix_n_any
          , public field_matrix_any_m
        {};


        // ********************************************************************************
        // Tag extraction
        // ********************************************************************************

#ifdef DOXYGEN

        //! Extracts the container tag from T.
        /**
         * \tparam T  The container type to tag.
         */
        template<typename T>
        struct container
        {
          //! The container tag associated with T.
          typedef ImplementationDefined type;
        };

#else // DOXYGEN

        // There is no standard implementation.
        template<typename T>
        struct container;


        template<typename Block, typename Alloc>
        struct container<Dune::BlockVector<Block,Alloc> >
        {
          typedef block_vector type;
        };


        // DynamicVector grew allocator support some time after the 2.3 release,
        // so we have to adjust the forward declaration accordingly

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,4)

        template<typename F, typename Allocator>
        struct container<DynamicVector<F,Allocator> >
        {
          typedef dynamic_vector type;
        };

#else

        template<typename F>
        struct container<DynamicVector<F> >
        {
          typedef dynamic_vector type;
        };

#endif

        template<typename F, int n>
        struct container<FieldVector<F,n> >
        {
          typedef field_vector_n type;
        };

        template<typename F>
        struct container<FieldVector<F,1> >
        {
          typedef field_vector_1 type;
        };


        template<typename Block, typename Alloc>
        struct container<Dune::BCRSMatrix<Block,Alloc> >
        {
          typedef bcrs_matrix type;
        };

        template<typename F>
        struct container<DynamicMatrix<F> >
        {
          typedef dynamic_matrix type;
        };

        template<typename F, int n, int m>
        struct container<FieldMatrix<F,n,m> >
        {
          typedef field_matrix_n_m type;
        };

        template<typename F, int n>
        struct container<FieldMatrix<F,n,1> >
        {
          typedef field_matrix_n_1 type;
        };

        template<typename F, int m>
        struct container<FieldMatrix<F,1,m> >
        {
          typedef field_matrix_1_m type;
        };

        template<typename F>
        struct container<FieldMatrix<F,1,1> >
        {
          typedef field_matrix_1_1 type;
        };

#endif // DOXYGEN

      } // namespace tags

      //! Gets instance of container tag associated with T
      /**
       * Returns an instance of the container tag for T. This function
       * is convenient when doing function-based tag dispatch, as it
       * saves on a lot of typing.
       *
       * \tparam T  The container for which to return a tag.
       * \returns   An instance of the associated container tag.
       */
      template<typename T>
      typename tags::container<T>::type container_tag(const T&)
      {
        return typename tags::container<T>::type();
      }

    } // namespace ISTL

  } // namespace PDELab
} // namespace Dune



#endif // DUNE_PDELAB_BACKEND_ISTL_TAGS_HH
