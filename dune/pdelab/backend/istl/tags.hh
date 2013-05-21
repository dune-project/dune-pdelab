// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_TAGS_HH
#define DUNE_PDELAB_BACKEND_ISTL_TAGS_HH

#include <dune/common/static_assert.hh>
#include <dune/common/documentation.hh>
#include <dune/pdelab/backend/istl/forwarddeclarations.hh>
#include <cstddef>

namespace Dune {

  namespace PDELab {

    namespace istl {

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
        struct container<BlockVector<Block,Alloc> >
        {
          typedef block_vector type;
        };

        template<typename F>
        struct container<DynamicVector<F> >
        {
          typedef dynamic_vector type;
        };

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
        struct container<BCRSMatrix<Block,Alloc> >
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

#ifndef DOXYGEN

      namespace impl {

        template<typename T, std::size_t depth, typename Tag>
        struct nesting_depth;

        template<typename T, std::size_t depth>
        struct nesting_depth<T,depth,tags::block_vector>
          : public nesting_depth<typename T::block_type,depth+1,typename tags::container<typename T::block_type>::type::base_tag>
        {};

        template<typename T, std::size_t depth>
        struct nesting_depth<T,depth,tags::dynamic_vector>
          : public integral_constant<std::size_t,depth+1>
        {};

        template<typename T, std::size_t depth>
        struct nesting_depth<T,depth,tags::field_vector>
          : public integral_constant<std::size_t,depth+1>
        {};

        template<typename T, std::size_t depth>
        struct nesting_depth<T,depth,tags::bcrs_matrix>
          : public nesting_depth<typename T::block_type,depth+1,typename tags::container<typename T::block_type>::type::base_tag>
        {};

        template<typename T, std::size_t depth>
        struct nesting_depth<T,depth,tags::dynamic_matrix>
          : public integral_constant<std::size_t,depth+1>
        {};

        template<typename T, std::size_t depth>
        struct nesting_depth<T,depth,tags::field_matrix>
          : public integral_constant<std::size_t,depth+1>
        {};

      }

#endif // DOXYGEN

      //! TMP for figuring out the depth up to which ISTL containers are nested.
      /**
       * This TMP calculates the nesting depth of ISTL containers. A FieldVector or
       * FieldMatrix has a depth of 1.
       */
      template<typename T>
      struct nesting_depth
        : public impl::nesting_depth<T,0,typename tags::container<T>::type::base_tag>
      {};

    } // namespace istl

  } // namespace PDELab
} // namespace Dune



#endif // DUNE_PDELAB_BACKEND_ISTL_TAGS_HH
