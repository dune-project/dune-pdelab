// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_UTILITY_HH
#define DUNE_PDELAB_BACKEND_ISTL_UTILITY_HH

#include <dune/common/typetraits.hh>
#include <dune/common/deprecated.hh>

#include <dune/pdelab/backend/istl/tags.hh>

namespace Dune {

  namespace PDELab {

    namespace istl {

      // ********************************************************************************
      // Helper functions for uniform access to ISTL containers
      //
      // The following suite of raw() functions should be used in places where an
      // algorithm might work on either the bare ISTL container or the PDELab
      // wrapper and has to access the bare container.
      // ********************************************************************************

      //! Returns the raw ISTL object associated with v, or v itself it is already an ISTL object.
      template<typename V>
      DUNE_DEPRECATED_MSG("raw() is deprecated and will be removed after PDELab 2.4. Use Backend::native() instead")
      V& raw(V& v)
      {
        return v;
      }

      //! Returns the raw ISTL object associated with v, or v itself it is already an ISTL object.
      template<typename V>
      DUNE_DEPRECATED_MSG("raw() is deprecated and will be removed after PDELab 2.4. Use Backend::native() instead")
      const V& raw(const V& v)
      {
        return v;
      }

      //! Returns the raw ISTL type associated with C, or C itself it is already an ISTL type.
      template<typename C>
      DUNE_DEPRECATED_MSG("raw_type<> is deprecated and will be removed after PDELab 2.4. Use Backend::Native<> instead")
      struct raw_type
      {
        typedef C type;
      };

#ifndef DOXYGEN

      template<typename GFS, typename C>
      typename ISTLBlockVectorContainer<GFS,C>::Container&
      raw(ISTLBlockVectorContainer<GFS,C>& v)
      {
        return v.base();
      }

      template<typename GFS, typename C>
      const typename ISTLBlockVectorContainer<GFS,C>::Container&
      raw(const ISTLBlockVectorContainer<GFS,C>& v)
      {
        return v.base();
      }

      template<typename GFSU, typename GFSV, typename C, typename Stats>
      typename ISTLMatrixContainer<GFSU,GFSV,C,Stats>::Container&
      raw(ISTLMatrixContainer<GFSU,GFSV,C,Stats>& m)
      {
        return m.base();
      }

      template<typename GFSU, typename GFSV, typename C, typename Stats>
      const typename ISTLMatrixContainer<GFSU,GFSV,C,Stats>::Container&
      raw(const ISTLMatrixContainer<GFSU,GFSV,C,Stats>& m)
      {
        return m.base();
      }

      template<typename GFS, typename C>
      struct raw_type<ISTLBlockVectorContainer<GFS,C> >
      {
        typedef C type;
      };

      template<typename GFSU, typename GFSV, typename C, typename Stats>
      struct raw_type<ISTLMatrixContainer<GFSU,GFSV,C,Stats> >
      {
        typedef C type;
      };


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

#endif // DUNE_PDELAB_BACKEND_ISTL_UTILITY_HH
