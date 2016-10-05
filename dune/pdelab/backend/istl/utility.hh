// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_UTILITY_HH
#define DUNE_PDELAB_BACKEND_ISTL_UTILITY_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <dune/common/typetraits.hh>
#include <dune/common/deprecated.hh>

#include <dune/pdelab/backend/istl/tags.hh>

namespace Dune {

  namespace PDELab {

    namespace ISTL {

#ifndef DOXYGEN

      // ********************************************************************************
      // Helpers for the nesting_depth TMP
      // ********************************************************************************

      namespace impl {

        template<typename T, std::size_t depth, typename Tag>
        struct nesting_depth;

        template<typename T, std::size_t depth>
        struct nesting_depth<T,depth,tags::block_vector>
          : public nesting_depth<typename T::block_type,depth+1,typename tags::container<typename T::block_type>::type::base_tag>
        {};

        template<typename T, std::size_t depth>
        struct nesting_depth<T,depth,tags::dynamic_vector>
          : public std::integral_constant<std::size_t,depth+1>
        {};

        template<typename T, std::size_t depth>
        struct nesting_depth<T,depth,tags::field_vector>
          : public std::integral_constant<std::size_t,depth+1>
        {};

        template<typename T, std::size_t depth>
        struct nesting_depth<T,depth,tags::bcrs_matrix>
          : public nesting_depth<typename T::block_type,depth+1,typename tags::container<typename T::block_type>::type::base_tag>
        {};

        template<typename T, std::size_t depth>
        struct nesting_depth<T,depth,tags::dynamic_matrix>
          : public std::integral_constant<std::size_t,depth+1>
        {};

        template<typename T, std::size_t depth>
        struct nesting_depth<T,depth,tags::field_matrix>
          : public std::integral_constant<std::size_t,depth+1>
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

    } // namespace ISTL
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_UTILITY_HH
