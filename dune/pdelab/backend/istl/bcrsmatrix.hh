// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_BCRSMATRIX_HH
#define DUNE_PDELAB_BACKEND_ISTL_BCRSMATRIX_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <dune/common/typetraits.hh>
#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/common/aliasedmatrixview.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>

namespace Dune {
  namespace PDELab {

    namespace ISTL {

      template<typename GFSV, typename GFSU, typename C, typename Stats>
      class BCRSMatrix
        : public Backend::impl::Wrapper<C>
      {

        friend Backend::impl::Wrapper<C>;

      public:

        typedef typename C::field_type ElementType;
        typedef ElementType E;
        typedef C Container;
        typedef typename C::field_type field_type;
        typedef typename C::block_type block_type;
        typedef typename C::size_type size_type;

        typedef GFSU TrialGridFunctionSpace;
        typedef GFSV TestGridFunctionSpace;

        typedef typename GFSV::Ordering::Traits::ContainerIndex RowIndex;
        typedef typename GFSU::Ordering::Traits::ContainerIndex ColIndex;

        typedef typename ISTL::build_pattern_type<C,GFSV,GFSU,typename GFSV::Ordering::ContainerAllocationTag>::type Pattern;

        typedef Stats PatternStatistics;

        using value_type = E;

#ifndef DOXYGEN

        // some trickery to avoid exposing average users to the fact that there might
        // be multiple statistics objects
        typedef typename std::conditional<
          (C::blocklevel > 2),
          std::vector<PatternStatistics>,
          PatternStatistics
          >::type StatisticsReturnType;

#endif // DOXYGEN

        template<typename RowCache, typename ColCache>
        using LocalView = UncachedMatrixView<BCRSMatrix,RowCache,ColCache>;

        template<typename RowCache, typename ColCache>
        using ConstLocalView = ConstUncachedMatrixView<const BCRSMatrix,RowCache,ColCache>;

        template<typename RowCache, typename ColCache>
        using AliasedLocalView = AliasedMatrixView<BCRSMatrix,RowCache,ColCache>;

        template<typename RowCache, typename ColCache>
        using ConstAliasedLocalView = ConstAliasedMatrixView<const BCRSMatrix,RowCache,ColCache>;

        template<typename GO>
        explicit BCRSMatrix (const GO& go)
          : _container(std::make_shared<Container>())
        {
          _stats = go.matrixBackend().buildPattern(go,*this);
        }

        /** \brief Construct matrix container using an externally given matrix as storage
         *
         * \tparam GO GridOperator type used to assemble into the matrix
         *
         * \param go GridOperator object used to assemble into the matrix
         * \param container ISTL matrix type that stores the actual data
         *
         * This BCRSMatrix constructor will reassemble the matrix occupation pattern.
         */
        template<typename GO>
        BCRSMatrix (const GO& go, Container& container)
          : _container(Dune::stackobject_to_shared_ptr(container))
        {
          _stats = go.matrixBackend().buildPattern(go,*this);
        }

        template<typename GO>
        BCRSMatrix (const GO& go, const E& e)
          : _container(std::make_shared<Container>())
        {
          _stats = go.matrixBackend().buildPattern(go,*this);
          (*_container) = e;
        }

        //! Creates an BCRSMatrix without allocating an underlying ISTL matrix.
        explicit BCRSMatrix (Backend::unattached_container = Backend::unattached_container())
        {}

        //! Creates an BCRSMatrix with an empty underlying ISTL matrix.
        explicit BCRSMatrix (Backend::attached_container)
          : _container(std::make_shared<Container>())
        {}

        BCRSMatrix(const BCRSMatrix& rhs)
          : _container(std::make_shared<Container>(*(rhs._container)))
        {}

        BCRSMatrix& operator=(const BCRSMatrix& rhs)
        {
          if (this == &rhs)
            return *this;
          _stats.clear();
          if (attached())
            {
              (*_container) = (*(rhs._container));
            }
          else
            {
              _container = std::make_shared<Container>(*(rhs._container));
            }
          return *this;
        }

        //! Returns pattern statistics for all contained BCRSMatrix objects.
        const StatisticsReturnType& patternStatistics() const
        {
          return patternStatistics(std::integral_constant<bool,(C::blocklevel > 2)>());
        }

#ifndef DOXYGEN

      private:

        const PatternStatistics& patternStatistics(std::false_type multiple) const
        {
          if (_stats.empty())
            DUNE_THROW(InvalidStateException,"no pattern statistics available");
          return _stats[0];
        }

        const std::vector<PatternStatistics>& patternStatistics(std::true_type multiple) const
        {
          if (_stats.empty())
            DUNE_THROW(InvalidStateException,"no pattern statistics available");
          return _stats;
        }

      public:

#endif

        void detach()
        {
          _container.reset();
          _stats.clear();
        }

        void attach(std::shared_ptr<Container> container)
        {
          _container = container;
        }

        bool attached() const
        {
          return bool(_container);
        }

        const std::shared_ptr<Container>& storage() const
        {
          return _container;
        }

        size_type N() const
        {
          return _container->N();
        }

        size_type M() const
        {
          return _container->M();
        }

        BCRSMatrix& operator= (const E& e)
        {
          (*_container) = e;
          return *this;
        }

        BCRSMatrix& operator*= (const E& e)
        {
          (*_container) *= e;
          return *this;
        }

        E& operator()(const RowIndex& ri, const ColIndex& ci)
        {
          return ISTL::access_matrix_element(ISTL::container_tag(*_container),*_container,ri,ci,ri.size()-1,ci.size()-1);
        }

        const E& operator()(const RowIndex& ri, const ColIndex& ci) const
        {
          return ISTL::access_matrix_element(ISTL::container_tag(*_container),*_container,ri,ci,ri.size()-1,ci.size()-1);
        }

      private:

        const Container& native() const
        {
          return *_container;
        }

        Container& native()
        {
          return *_container;
        }

      public:

        template<typename RowCache, typename ColCache>
        value_type* data(const RowCache& row_cache, const ColCache& col_cache)
        {
          return &((*this)(row_cache.containerIndex(0),col_cache.containerIndex(0)));
        }

        template<typename RowCache, typename ColCache>
        const value_type* data(const RowCache& row_cache, const ColCache& col_cache) const
        {
          return &((*this)(row_cache.containerIndex(0),col_cache.containerIndex(0)));
        }

        void flush()
        {}

        void finalize()
        {}

        void clear_row(const RowIndex& ri, const E& diagonal_entry)
        {
          ISTL::clear_matrix_row(ISTL::container_tag(*_container),*_container,ri,ri.size()-1);
          ISTL::write_matrix_element_if_exists(diagonal_entry,ISTL::container_tag(*_container),*_container,ri,ri,ri.size()-1,ri.size()-1);
        }

        void clear_row_block(const RowIndex& ri, const E& diagonal_entry)
        {
          ISTL::clear_matrix_row_block(ISTL::container_tag(*_container),*_container,ri,ri.size()-1);
          ISTL::write_matrix_element_if_exists_to_block(diagonal_entry,ISTL::container_tag(*_container),*_container,ri,ri,ri.size()-1,ri.size()-1);
        }

      private:

        std::shared_ptr<Container> _container;
        std::vector<PatternStatistics> _stats;

      };

    } // namespace ISTL

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_BCRSMATRIX_HH
