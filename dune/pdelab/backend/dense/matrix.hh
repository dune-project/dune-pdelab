// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_DENSE_MATRIX_HH
#define DUNE_PDELAB_BACKEND_DENSE_MATRIX_HH

#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>

#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/pdelab/backend/tags.hh>
#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/backend/dense/descriptors.hh>

namespace Dune {
  namespace PDELab {
    namespace dense {

      template<typename GFSV, typename GFSU, typename C>
      class MatrixContainer
      {

      public:

        typedef C Container;
        typedef typename Container::value_type ElementType;
        typedef ElementType E;

        typedef ElementType field_type;
        typedef typename Container::size_type size_type;

        typedef typename GFSV::Ordering::Traits::ContainerIndex RowIndex;
        typedef typename GFSU::Ordering::Traits::ContainerIndex ColIndex;

        // We have to provide some kind of type for the pattern to make the GridOperator happy
        typedef int Pattern;

        template<typename RowCache, typename ColCache>
        class LocalView
        {

          dune_static_assert((is_same<typename RowCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFSV>::value),
                             "The RowCache passed to LocalView must belong to the underlying GFSV");

          dune_static_assert((is_same<typename ColCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFSU>::value),
                             "The ColCache passed to LocalView must belong to the underlying GFSU");

        public:

          typedef typename Container::value_type E;
          typedef E ElementType;

          typedef RowCache RowIndexCache;
          typedef ColCache ColIndexCache;

          typedef typename RowCache::LocalFunctionSpace LFSV;
          typedef typename ColCache::LocalFunctionSpace LFSU;

          typedef typename LFSV::Traits::DOFIndex RowDOFIndex;
          typedef typename LFSV::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex RowContainerIndex;

          typedef typename LFSU::Traits::DOFIndex ColDOFIndex;
          typedef typename LFSU::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex ColContainerIndex;

          LocalView()
            : _container(nullptr)
            , _row_cache(nullptr)
            , _col_cache(nullptr)
          {}

          LocalView(MatrixContainer& container)
            : _container(&container)
            , _row_cache(nullptr)
            , _col_cache(nullptr)
          {}

          const RowIndexCache& rowIndexCache() const
          {
            assert(_row_cache);
            return *_row_cache;
          }

          const ColIndexCache& colIndexCache() const
          {
            assert(_col_cache);
            return *_col_cache;
          }

          void attach(MatrixContainer& container)
          {
            _container = &container;
          }

          void detach()
          {
            _container = nullptr;
          }

          void bind(const RowCache& row_cache, const ColCache& col_cache)
          {
            _row_cache = &row_cache;
            _col_cache = &col_cache;
          }

          void unbind()
          {
          }

          size_type N() const
          {
            return _row_cache->size();
          }

          size_type M() const
          {
            return _col_cache->size();
          }

          template<typename LC>
          void read(LC& local_container) const
          {
            for (size_type i = 0; i < N(); ++i)
              for (size_type j = 0; j < M(); ++j)
                local_container.getEntry(i,j) = (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
          }

          template<typename LC>
          void write(const LC& local_container) const
          {
            for (size_type i = 0; i < N(); ++i)
              for (size_type j = 0; j < M(); ++j)
                (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j)) = local_container.getEntry(i,j);
          }

          template<typename LC>
          void add(const LC& local_container) const
          {
            for (size_type i = 0; i < N(); ++i)
              for (size_type j = 0; j < M(); ++j)
                (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j)) += local_container.getEntry(i,j);
          }

          void commit()
          {
          }


          ElementType& operator()(size_type i, size_type j)
          {
            return (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
          }

          const ElementType& operator()(size_type i, size_type j) const
          {
            return (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
          }

          ElementType& operator()(const RowDOFIndex& i, const ColDOFIndex& j)
          {
            return (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
          }

          const ElementType& operator()(const RowDOFIndex& i, const ColDOFIndex& j) const
          {
            return (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
          }

          ElementType& operator()(const RowContainerIndex& i, const ColContainerIndex& j)
          {
            return (*_container)(i,j);
          }

          const ElementType& operator()(const RowContainerIndex& i, const ColContainerIndex& j) const
          {
            return (*_container)(i,j);
          }

          ElementType& operator()(const RowContainerIndex& i, size_type j)
          {
            return (*_container)(i,_col_cache->containerIndex(j));
          }

          const ElementType& operator()(const RowContainerIndex& i, size_type j) const
          {
            return (*_container)(i,_col_cache->containerIndex(j));
          }

          ElementType& operator()(size_type i, const ColContainerIndex& j)
          {
            return (*_container)(_row_cache->containerIndex(i),j);
          }

          const ElementType& operator()(size_type i, const ColContainerIndex& j) const
          {
            return (*_container)(_row_cache->containerIndex(i),j);
          }


          void add(size_type i, size_type j, const ElementType& v)
          {
            (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j)) += v;
          }

          void add(const RowDOFIndex& i, const ColDOFIndex& j, const ElementType& v)
          {
            (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j)) += v;
          }

          void add(const RowContainerIndex& i, const ColContainerIndex& j, const ElementType& v)
          {
            (*_container)(i,j) += v;
          }

          void add(const RowContainerIndex& i, size_type j, const ElementType& v)
          {
            (*_container)(i,_col_cache->containerIndex(j)) += v;
          }

          void add(size_type i, const ColContainerIndex& j, const ElementType& v)
          {
            (*_container)(_row_cache->containerIndex(i),j) += v;
          }

          MatrixContainer& global_container()
          {
            return *_container;
          }

          const MatrixContainer& global_container() const
          {
            return *_container;
          }

        private:

          MatrixContainer* _container;
          const RowCache* _row_cache;
          const ColCache* _col_cache;

        };

        template<typename GO>
        MatrixContainer(const GO& go)
          : _rows(go.testGridFunctionSpace().size())
          , _cols(go.trialGridFunctionSpace().size())
          , _container(make_shared<Container>(_rows*_cols,E(0)))
        {}

        template<typename GO>
        MatrixContainer(const GO& go, const E& e)
          : _rows(go.testGridFunctionSpace().size())
          , _cols(go.trialGridFunctionSpace().size())
          , _container(make_shared<Container>(_rows*_cols,e))
        {}

        //! Creates an ISTLMatrixContainer without allocating an underlying ISTL matrix.
        explicit MatrixContainer(tags::unattached_container = tags::unattached_container())
          : _rows(0)
          , _cols(0)
        {}

        //! Creates an ISTLMatrixContainer with an empty underlying ISTL matrix.
        explicit MatrixContainer(tags::attached_container)
          : _rows(0)
          , _cols(0)
          , _container(make_shared<Container>())
        {}

        MatrixContainer(const MatrixContainer& rhs)
          : _rows(rhs._rows)
          , _cols(rhs._cols)
          , _container(make_shared<Container>(*(rhs._container)))
        {}

        MatrixContainer& operator=(const MatrixContainer& rhs)
        {
          if (this == &rhs)
            return *this;
          if (_rows == 0 && _cols == 0)
            {
              _rows = rhs._rows;
              _cols = rhs._cols;
            }
          if (attached())
            {
              (*_container) = (*(rhs._container));
            }
          else
            {
              _container = make_shared<Container>(*(rhs._container));
            }
          return *this;
        }

        void detach()
        {
          _container.reset();
        }

        void attach(shared_ptr<Container> container)
        {
          _container = container;
        }

        bool attached() const
        {
          return bool(_container);
        }

        const shared_ptr<Container>& storage() const
        {
          return _container;
        }

        size_type N() const
        {
          return _rows;
        }

        size_type M() const
        {
          return _cols;
        }

        MatrixContainer& operator=(const E& e)
        {
          std::fill(_container->begin(),_container->end(),e);
          return *this;
        }

        MatrixContainer& operator*=(const E& e)
        {
          using namespace std::placeholders;
          std::transform(_container->begin(),_container->end(),_container->begin(),std::bind(std::multiplies<E>(),e,_1));
          return *this;
        }

        template<typename V>
        void mv(const V& x, V& y) const
        {
          auto rowit = _container->begin();
          for (auto& v : y)
            {
              v = std::inner_product(rowit,rowit + _cols,x.begin(),E(0));
              rowit += _cols;
            }
        }

        template<typename V>
        void usmv(const E alpha, const V& x, V& y) const
        {
          auto rowit = _container->begin();
          for (auto& v : y)
            {
              v += alpha * std::inner_product(rowit,rowit + _cols,x.begin(),E(0));
              rowit += _cols;
            }
        }

        E& operator()(const RowIndex& ri, const ColIndex& ci)
        {
          return (*_container)[ri[0]*_cols + ci[0]];
        }

        const E& operator()(const RowIndex& ri, const ColIndex& ci) const
        {
          return (*_container)[ri[0]*_cols + ci[0]];
        }

        const Container& base() const
        {
          return *_container;
        }

        Container& base()
        {
          return *_container;
        }

        void flush()
        {}

        void finalize()
        {}

        void clear_row(const RowIndex& ri, const E& diagonal_entry)
        {
          std::fill(_container->begin() + ri[0]*_cols,_container->begin() + (ri[0]+1)*_cols,E(0));
          (*this)(ri,ri) = diagonal_entry;
        }

      private:

        std::size_t _rows;
        std::size_t _cols;
        shared_ptr<Container> _container;

      };

    } // namespace dense

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_DENSE_MATRIX_HH
