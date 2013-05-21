// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH

#include <dune/pdelab/backend/tags.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>

namespace Dune {
  namespace PDELab {

    template<typename GFSV, typename GFSU, typename C>
    class ISTLMatrixContainer
    {

    public:

      typedef typename C::field_type ElementType;
      typedef ElementType E;
      typedef C Container;
      typedef C BaseT;
      typedef typename C::field_type field_type;
      typedef typename C::block_type block_type;
      typedef typename C::size_type size_type;

      typedef typename GFSV::Ordering::Traits::ContainerIndex RowIndex;
      typedef typename GFSU::Ordering::Traits::ContainerIndex ColIndex;

      typedef typename istl::build_pattern_type<C,GFSV,GFSU,typename GFSV::Ordering::ContainerAllocationTag>::type Pattern;

      template<typename RowCache, typename ColCache>
      class LocalView
      {

        dune_static_assert((is_same<typename RowCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFSV>::value),
                           "The RowCache passed to LocalView must belong to the underlying GFSV");

        dune_static_assert((is_same<typename ColCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFSU>::value),
                           "The ColCache passed to LocalView must belong to the underlying GFSU");

      public:

        typedef typename C::field_type E;
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

        LocalView(ISTLMatrixContainer& container)
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

        void attach(ISTLMatrixContainer& container)
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

        ISTLMatrixContainer& global_container()
        {
          return *_container;
        }

        const ISTLMatrixContainer& global_container() const
        {
          return *_container;
        }

      private:

        ISTLMatrixContainer* _container;
        const RowCache* _row_cache;
        const ColCache* _col_cache;

      };

      template<typename GO>
      ISTLMatrixContainer (const GO& go)
        : _container(make_shared<Container>())
      {
        Pattern pattern(go.testGridFunctionSpace().ordering(),go.trialGridFunctionSpace().ordering());
        go.fill_pattern(pattern);
        allocate_matrix(go.testGridFunctionSpace().ordering(),
                        go.trialGridFunctionSpace().ordering(),
                        pattern,
                        *_container);
      }

      template<typename GO>
      ISTLMatrixContainer (const GO& go, const E& e)
        : _container(make_shared<Container>())
      {
        Pattern pattern(go.testGridFunctionSpace().ordering(),go.trialGridFunctionSpace().ordering());
        go.fill_pattern(pattern);
        allocate_matrix(go.testGridFunctionSpace().ordering(),
                        go.trialGridFunctionSpace().ordering(),
                        pattern,
                        *_container);
        _container = e;
      }


      //! Creates an ISTLMatrixContainer without allocating an underlying ISTL matrix.
      explicit ISTLMatrixContainer (tags::unattached_container = tags::unattached_container())
      {}

      //! Creates an ISTLMatrixContainer with an empty underlying ISTL matrix.
      explicit ISTLMatrixContainer (tags::attached_container)
      {}

      ISTLMatrixContainer(const ISTLMatrixContainer& rhs)
        : _container(make_shared<Container>(*(rhs._container)))
      {}

      ISTLMatrixContainer& operator=(const ISTLMatrixContainer& rhs)
      {
        if (this == &rhs)
          return *this;
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
        return _container->N();
      }

      size_type M() const
      {
        return _container->M();
      }

      ISTLMatrixContainer& operator= (const E& e)
      {
        (*_container) = e;
        return *this;
      }

      ISTLMatrixContainer& operator*= (const E& e)
      {
        (*_container) *= e;
        return *this;
      }

      E& operator()(const RowIndex& ri, const ColIndex& ci)
      {
        return istl::access_matrix_element(istl::container_tag(*_container),*_container,ri,ci,ri.size()-1,ci.size()-1);
      }

      const E& operator()(const RowIndex& ri, const ColIndex& ci) const
      {
        return istl::access_matrix_element(istl::container_tag(*_container),*_container,ri,ci,ri.size()-1,ci.size()-1);
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
        istl::clear_matrix_row(istl::container_tag(*_container),*_container,ri,ri.size()-1);
        (*this)(ri,ri) = diagonal_entry;
      }

    private:

      shared_ptr<Container> _container;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH
