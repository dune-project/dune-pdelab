// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_COMMON_THREADEDMATRIXVIEW_HH
#define DUNE_PDELAB_BACKEND_COMMON_THREADEDMATRIXVIEW_HH

#include <cstddef>
#include <mutex>
#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
  namespace PDELab {

    template<typename RowIndex, typename ColIndex, typename Value>
    struct MatrixBufferEntry {
      MatrixBufferEntry(const RowIndex &row_, const ColIndex &column_,
                        const Value &value_) :
        row(row_), column(column_), value(value_)
      { }

      RowIndex row;
      ColIndex column;
      Value value;
    };

    template<typename M, typename Mutex>
    class ThreadedMatrixBuffer
    {
    public:
      typedef M Container;

      typedef typename Container::ElementType value_type;
      typedef typename Container::RowIndex row_index_type;
      typedef typename Container::ColIndex col_index_type;

      class reference;

      ThreadedMatrixBuffer(M& container, Mutex &mutex)
        : _container(container)
        , _mutex(mutex)
        , _autocommit_after(0)
      { }

      ThreadedMatrixBuffer(const ThreadedMatrixBuffer &other)
        : _container(other._container)
        , _mutex(other._mutex)
        , _autocommit_after(other._autocommit_after)
      {
        if(_autocommit_after == 0)
          _add_buffer.reserve(other._add_buffer.capacity());
        else
          _add_buffer.reserve(_autocommit_after);
      }

      void set_autocommit(std::size_t autocommit)
      {
        _autocommit_after = autocommit;
        _add_buffer.reserve(_autocommit_after);
      }

      void commit()
      {
        if(_add_buffer.empty())
          return;
        {
          std::lock_guard<Mutex> guard(_mutex);
          for(const auto &entry : _add_buffer)
            container()(entry.row, entry.column) += entry.value;
        }
        _add_buffer.clear();
      }

      const Container& container() const
      {
        return _container;
      }

      Container& container()
      {
        return _container;
      }

      void add(const row_index_type &row, const col_index_type &column,
               const value_type &value)
      {
        if(_autocommit_after > 0 && _add_buffer.size() >= _autocommit_after)
          commit();
        _add_buffer.emplace_back(row, column, value);
      }

      reference operator()(const row_index_type& row,
                           const col_index_type &column)
      {
        return reference(*this, row, column);
      }

    private:
      M &_container;
      std::vector<MatrixBufferEntry<row_index_type, col_index_type,
                                    value_type> > _add_buffer;
      Mutex &_mutex;
      std::size_t _autocommit_after;
    };

    template<typename M, typename Mutex>
    class ThreadedMatrixBuffer<M, Mutex>::reference
    {
      ThreadedMatrixBuffer &buffer_;
      row_index_type row_;
      col_index_type column_;

    public:
      reference(ThreadedMatrixBuffer &buffer, const row_index_type &row,
                const col_index_type &column) :
        buffer_(buffer), row_(row), column_(column)
      { }

      reference &operator+=(const value_type &value)
      {
        buffer_.add(row_, column_, value);
        return *this;
      }
      reference &operator-=(const value_type &value)
      {
        buffer_.add(row_, column_, -value);
        return *this;
      }
    };

    //! A threadsave view for PDELab matrix containers.
    /**
     * The only action supported is add() and friends.  Values to be added are
     * collected in an internal std::vector and commited when the autocommit
     * limit is reached, or commit is called explicitly.
     *
     * view[index] += value (and view[index] -= value) is supported with the
     * help of a proxy object.
     */
    template<typename M_, typename RowCache, typename ColCache,
             typename Mutex = std::mutex>
    class ThreadedMatrixView
    {

    public:
      typedef typename remove_const<M_>::type Container;
      typedef ThreadedMatrixBuffer<Container, Mutex> Buffer;
      typedef typename Buffer::reference Proxy;

      static_assert(
        (is_same<
           typename RowCache::LocalFunctionSpace::Traits::GridFunctionSpace,
           typename Container::TestGridFunctionSpace
           >::value),
        "The RowCache passed to ThreadedMatrixView must belong to the "
        "underlying GFSV"
        );

      static_assert(
        (is_same<
           typename ColCache::LocalFunctionSpace::Traits::GridFunctionSpace,
           typename Container::TrialGridFunctionSpace
           >::value),
        "The ColCache passed to ThreadedMatrixView must belong to the "
        "underlying GFSU"
        );

    public:

      typedef typename Container::field_type E;
      typedef typename Container::size_type size_type;

      typedef E ElementType;

      typedef RowCache RowIndexCache;
      typedef ColCache ColIndexCache;

      typedef typename RowCache::LocalFunctionSpace LFSV;
      typedef typename ColCache::LocalFunctionSpace LFSU;

      typedef typename LFSV::Traits::DOFIndex RowDOFIndex;
      typedef typename LFSV::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex RowContainerIndex;

      typedef typename LFSU::Traits::DOFIndex ColDOFIndex;
      typedef typename LFSU::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex ColContainerIndex;

      ThreadedMatrixView()
        : _row_cache(nullptr)
        , _col_cache(nullptr)
      {}

      ThreadedMatrixView(const shared_ptr<Buffer>& buffer)
        : _row_cache(nullptr)
        , _col_cache(nullptr)
        , _buffer(buffer)
      {}

      void commit()
      {
        _buffer->commit();
      }

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

      void attach(const shared_ptr<Buffer>& buffer)
      {
        _buffer = buffer;
      }

      void detach()
      {
        _buffer.reset();
      }

      void bind(const RowCache& row_cache, const ColCache& col_cache)
      {
        _row_cache = &row_cache;
        _col_cache = &col_cache;
      }

      void unbind()
      {}

      size_type N() const
      {
        return rowIndexCache().size();
      }

      size_type M() const
      {
        return colIndexCache().size();
      }

      // operations on local_container
      template<typename LC>
      void read(LC& local_container) const
      {
        static_assert(AlwaysFalse<LC>::value, "read not supported for "
                      "ThreadedMatrixView");
      }

      template<typename LC>
      void write(const LC& local_container)
      {
        static_assert(AlwaysFalse<LC>::value, "write not supported for "
                      "ThreadedMatrixView");
      }

      template<typename LC>
      void add(const LC& local_container)
      {
        for (size_type i = 0; i < N(); ++i)
          for (size_type j = 0; j < M(); ++j)
            _buffer->add(rowIndexCache().containerIndex(i),
                         colIndexCache().containerIndex(j),
                         local_container.getEntry(i,j));
      }

      // subscription
      template<typename RowIndexType, typename ColIndexType>
      const ElementType &operator()(const RowIndexType &row,
                                    const ColIndexType &column) const
      {
        static_assert(AlwaysFalse<RowIndexType>::value, "element access "
                      "(read) not supported for ThreadedMatrixView");
      }

      Proxy operator()(const RowContainerIndex& i, const ColContainerIndex& j)
      {
        return (*_buffer)(i, j);
      }

      template<typename RowIndexType>
      Proxy operator()(const RowIndexType& i, const ColDOFIndex& j)
      {
        return (*this)(i, colIndexCache().containerIndex(j));
      }

      Proxy operator()(const RowDOFIndex& i, const ColContainerIndex& j)
      {
        return (*this)(rowIndexCache().containerIndex(i), j);
      }

      template<typename RowIndexType>
      Proxy operator()(const RowIndexType& i, size_type j)
      {
        return (*this)(i, colIndexCache().containerIndex(j));
      }

      Proxy operator()(size_type i, const ColContainerIndex& j)
      {
        return (*this)(rowIndexCache().containerIndex(i), j);
      }

      template<typename RowIndexType, typename ColIndexType>
      void add(const RowIndexType& i, const ColIndexType& j,
               const ElementType& v)
      {
        (*this)(i, j) += v;
      }

      const Container& container() const
      {
        return _buffer->container();
      }

      Container& container()
      {
        return _buffer->container();
      }

      const Container& global_container() const DUNE_DEPRECATED_MSG("global_container() is deprecated, use container() instead.")
      {
        return container();
      }

      Container& global_container() DUNE_DEPRECATED_MSG("global_container() is deprecated, use container() instead.")
      {
        return container();
      }

      const shared_ptr<Buffer> &buffer()
      {
        return _buffer;
      }

    protected:
      const RowCache* _row_cache;
      const ColCache* _col_cache;
      shared_ptr<Buffer> _buffer;
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_COMMON_THREADEDMATRIXVIEW_HH
