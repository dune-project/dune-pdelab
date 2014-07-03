// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_COMMON_THREADEDVECTORVIEW_HH
#define DUNE_PDELAB_BACKEND_COMMON_THREADEDVECTORVIEW_HH

#include <algorithm>
#include <cstddef>
#include <mutex>
#include <utility>
#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/shared_ptr.hh>

namespace Dune {
  namespace PDELab {

    template<typename V, typename Mutex>
    class ThreadedVectorViewBuffer
    {
    public:
      typedef V Container;

      typedef typename Container::ElementType value_type;
      typedef typename Container::ContainerIndex index_type;
      class reference;

      ThreadedVectorViewBuffer(V& container, Mutex &mutex)
        : _container(container)
        , _mutex(mutex)
        , _autocommit_after(0)
      { }

      ThreadedVectorViewBuffer(const ThreadedVectorViewBuffer &other)
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
          for(const auto &action : _add_buffer)
            container()[action.first] += action.second;
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

      void add(const index_type &index, const value_type &value)
      {
        if(_autocommit_after > 0 && _add_buffer.size() >= _autocommit_after)
          commit();
        _add_buffer.emplace_back(index, value);
      }

      reference operator[](index_type i)
      {
        return reference(*this, i);
      }

    private:
      V &_container;
      std::vector<std::pair<index_type, value_type> > _add_buffer;
      Mutex &_mutex;
      std::size_t _autocommit_after;
    };

    template<typename V, typename Mutex>
    class ThreadedVectorViewBuffer<V, Mutex>::reference
    {
      ThreadedVectorViewBuffer &buffer_;
      index_type index_;

    public:
      reference(ThreadedVectorViewBuffer &buffer, const index_type &index) :
        buffer_(buffer), index_(index)
      { }

      reference &operator+=(const value_type &value)
      {
        buffer_.add(index_, value);
      }
      reference &operator-=(const value_type &value)
      {
        buffer_.add(index_, -value);
      }
    };

    //! A threadsave view for PDELab vector containers.
    /**
     * The only action supported is add() and friends.  Values to be added are
     * collected in an internal std::vector and commited when the autocommit
     * limit is reached, or commit is called explicitly.
     *
     * view[index] += value (and view[index] -= value) is supported with the
     * help of a proxy object.
     */
    template<typename V, typename LFSC, typename Mutex = std::mutex>
    class ThreadedVectorView
    {
    public:
      typedef ThreadedVectorViewBuffer<V, Mutex> Buffer;

      typedef typename Buffer::reference Proxy;

      typedef V Container;
      typedef LFSC LFSCache;

      typedef typename Container::ElementType ElementType;
      typedef typename Container::size_type size_type;
      typedef typename LFSCache::DOFIndex DOFIndex;
      typedef typename LFSCache::ContainerIndex ContainerIndex;

      ThreadedVectorView()
        : _lfs_cache(nullptr)
        , _buffer(nullptr)
      { }

      ThreadedVectorView(const shared_ptr<Buffer> &buffer)
        : _lfs_cache(nullptr)
        , _buffer(buffer)
      { }

      void attach(const shared_ptr<Buffer> &buffer)
      {
        _buffer = buffer;
      }

      void detach()
      {
        _buffer.reset();
      }

      void bind(const LFSCache& lfs_cache)
      {
        _lfs_cache = &lfs_cache;
      }

      void unbind()
      {
      }

      size_type size() const
      {
        return cache().size();
      }

      void commit()
      {
        _buffer->commit();
      }

      // operations on local_container
      template<typename LC>
      void read(LC& local_container) const
      {
        static_assert(AlwaysFalse<LC>::value, "read not supported for "
                      "ThreadedVectorView");
      }

      template<typename LC>
      void write(const LC& local_container)
      {
        static_assert(AlwaysFalse<LC>::value, "write not supported for "
                      "ThreadedVectorView");
      }

      template<typename LC>
      void add(const LC& local_container)
      {
        for (size_type i = 0; i < size(); ++i)
          _buffer->add(cache().containerIndex(i),
                       accessBaseContainer(local_container)[i]);
      }

      // operations on child_lfs and local_container
      template<typename ChildLFS, typename LC>
      void read(const ChildLFS& child_lfs, LC& local_container) const
      {
        static_assert(AlwaysFalse<LC>::value, "read not supported for "
                      "ThreadedVectorView");
      }

      template<typename ChildLFS, typename LC>
      void write(const ChildLFS& child_lfs, const LC& local_container)
      {
        static_assert(AlwaysFalse<LC>::value, "write not supported for "
                      "ThreadedVectorView");
      }

      template<typename ChildLFS, typename LC>
      void add(const ChildLFS& child_lfs, const LC& local_container)
      {
        for (size_type i = 0; i < child_lfs.size(); ++i)
        {
          const size_type local_index = child_lfs.localIndex(i);
          _buffer->add(cache().containerIndex(local_index),
                       accessBaseContainer(local_container)[local_index]);
        }
      }

      // sub-operations on child_lfs and local_container
      template<typename ChildLFS, typename LC>
      void read_sub_container(const ChildLFS& child_lfs, LC& local_container) const
      {
        static_assert(AlwaysFalse<LC>::value, "read_sub_container not "
                      "supported for ThreadedVectorView");
      }

      template<typename ChildLFS, typename LC>
      void write_sub_container(const ChildLFS& child_lfs, const LC& local_container)
      {
        static_assert(AlwaysFalse<LC>::value, "write_sub_container not "
                      "supported for ThreadedVectorView");
      }

      template<typename ChildLFS, typename LC>
      void add_sub_container(const ChildLFS& child_lfs, const LC& local_container)
      {
        for (size_type i = 0; i < child_lfs.size(); ++i)
        {
          const size_type local_index = child_lfs.localIndex(i);
          _buffer->add(cache().containerIndex(local_index),
                       accessBaseContainer(local_container)[i]);
        }
      }

      // subscription
      template<typename IndexType>
      const ElementType& operator[](const IndexType& i) const
      {
        static_assert(AlwaysFalse<IndexType>::value, "element access "
                      "(read) not supported for ThreadedVectorView");
      }

      Proxy operator[](size_type i)
      {
        return (*_buffer)[cache().containerIndex(i)];
      }


      Proxy operator[](const DOFIndex& di)
      {
        return (*_buffer)[cache().containerIndex(di)];
      }


      Proxy operator[](const ContainerIndex& ci)
      {
        return (*_buffer)[ci];
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
        return _buffer->container();
      }

      Container& global_container() DUNE_DEPRECATED_MSG("global_container() is deprecated, use container() instead.")
      {
        return _buffer->container();
      }

      const LFSCache& cache() const
      {
        return *_lfs_cache;
      }

      const shared_ptr<Buffer> &buffer()
      {
        return _buffer;
      }

    protected:

      const LFSCache* _lfs_cache;
      shared_ptr<Buffer> _buffer;
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_COMMON_THREADEDVECTORVIEW_HH
