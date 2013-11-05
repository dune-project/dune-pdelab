// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_UTILITY_HH
#define DUNE_PDELAB_BACKEND_ISTL_UTILITY_HH

#include <dune/common/typetraits.hh>
#include <dune/common/parametertree.hh>

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
      V& raw(V& v)
      {
        return v;
      }

      //! Returns the raw ISTL object associated with v, or v itself it is already an ISTL object.
      template<typename V>
      const V& raw(const V& v)
      {
        return v;
      }

      //! Returns the raw ISTL type associated with C, or C itself it is already an ISTL type.
      template<typename C>
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

      template<typename GFS, typename C>
      typename FlatVectorContainer<GFS,C>::Container&
      raw(FlatVectorContainer<GFS,C>& v)
      {
        return v.base();
      }

      template<typename GFS, typename C>
      const typename FlatVectorContainer<GFS,C>::Container&
      raw(const FlatVectorContainer<GFS,C>& v)
      {
        return v.base();
      }

      template<typename GFS, typename C>
      typename BlockVectorContainer<GFS,C>::Container&
      raw(BlockVectorContainer<GFS,C>& v)
      {
        return v.base();
      }

      template<typename GFS, typename C>
      const typename BlockVectorContainer<GFS,C>::Container&
      raw(const BlockVectorContainer<GFS,C>& v)
      {
        return v.base();
      }

      template<typename GFSU, typename GFSV, typename C>
      typename ISTLMatrixContainer<GFSU,GFSV,C>::Container&
      raw(ISTLMatrixContainer<GFSU,GFSV,C>& m)
      {
        return m.base();
      }

      template<typename GFSU, typename GFSV, typename C>
      const typename ISTLMatrixContainer<GFSU,GFSV,C>::Container&
      raw(const ISTLMatrixContainer<GFSU,GFSV,C>& m)
      {
        return m.base();
      }

      template<typename GFSU, typename GFSV, typename C>
      typename FlatELLMatrixContainer<GFSU,GFSV,C>::Container&
      raw(FlatELLMatrixContainer<GFSU,GFSV,C>& m)
      {
        return m.base();
      }

      template<typename GFSU, typename GFSV, typename C>
      const typename FlatELLMatrixContainer<GFSU,GFSV,C>::Container&
      raw(const FlatELLMatrixContainer<GFSU,GFSV,C>& m)
      {
        return m.base();
      }

      template<typename GFSU, typename GFSV, typename C>
      typename BELLMatrixContainer<GFSU,GFSV,C>::Container&
      raw(BELLMatrixContainer<GFSU,GFSV,C>& m)
      {
        return m.base();
      }

      template<typename GFSU, typename GFSV, typename C>
      const typename BELLMatrixContainer<GFSU,GFSV,C>::Container&
      raw(const BELLMatrixContainer<GFSU,GFSV,C>& m)
      {
        return m.base();
      }

      template<typename GFS, typename C>
      struct raw_type<ISTLBlockVectorContainer<GFS,C> >
      {
        typedef C type;
      };

      template<typename GFS, typename C>
      struct raw_type<FlatVectorContainer<GFS,C> >
      {
        typedef C type;
      };

      template<typename GFS, typename C>
      struct raw_type<BlockVectorContainer<GFS,C> >
      {
        typedef C type;
      };

      template<typename GFSU, typename GFSV, typename C>
      struct raw_type<ISTLMatrixContainer<GFSU,GFSV,C> >
      {
        typedef C type;
      };

      template<typename GFSU, typename GFSV, typename C>
      struct raw_type<FlatELLMatrixContainer<GFSU,GFSV,C> >
      {
        typedef C type;
      };

      template<typename GFSU, typename GFSV, typename C>
      struct raw_type<BELLMatrixContainer<GFSU,GFSV,C> >
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


      class MatrixParameters
      {

      public:

        typedef std::size_t size_type;

        explicit MatrixParameters(size_type entries_per_row, double overflow_fraction = 0.05)
          : _entries_per_row(entries_per_row)
          , _overflow_fraction(overflow_fraction)
        {}

        explicit MatrixParameters(const ParameterTree& parameters)
          : _entries_per_row(parameters.get<size_type>("entries_per_row"))
          , _overflow_fraction(parameters.get<double>("overflow_fraction",0.05))
        {}

        void setEntriesPerRow(size_type entries_per_row)
        {
          _entries_per_row = entries_per_row;
        }

        void setOverflowFraction(double overflow_fraction)
        {
          _overflow_fraction = overflow_fraction;
        }

        size_type entriesPerRow() const
        {
          return _entries_per_row;
        }

        double overflowFraction() const
        {
          return _overflow_fraction;
        }

      private:

        size_type _entries_per_row;
        double _overflow_fraction;

      };

    } // namespace istl
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_UTILITY_HH
