#ifndef DUNE_PDELAB_COMMON_CONCURRENCY_SHARED_STASH_HH
#define DUNE_PDELAB_COMMON_CONCURRENCY_SHARED_STASH_HH

#include <functional>
#include <queue>
#include <memory>

#if HAVE_TBB
#include <oneapi/tbb/concurrent_queue.h>
#endif


namespace Dune::PDELab::inline Experimental {


// behaves like a unique pointer, however on destruction the object is not discarded but stored for later reuse.
// when this object is copied, it searches for a new stash, or if there is none, creates one anew.
// this is useful for multi-threaded code when each object must hold its own state but the factory stage is very expensive.
template<class T>
struct SharedStash
{
  struct Data {

    std::function<std::unique_ptr<T>()> _factory;
#ifdef HAVE_TBB
    oneapi::tbb::concurrent_queue<std::unique_ptr<T>> _stash_queue;
#else
    std::queue<std::unique_ptr<T>> _stash_queue;
    std::mutex _mutex;
#endif
  };

  SharedStash(std::function<std::unique_ptr<T>()> factory)
  {
    _data = std::make_shared<Data>();
    _data->_factory = std::move(factory);
    _stash = _data->_factory();
  }

  SharedStash(const SharedStash& other)
    : _data{other._data}
  {
    _stash = obtain_stash();
  }

  SharedStash(SharedStash&&) = default;

  ~SharedStash() { release_stash(); }

  T const * operator->() const {
    assert(_stash);
    return _stash.get();
  }

  T* operator->() {
    assert(_stash);
    return _stash.get();
  }

  explicit operator bool() const noexcept { return bool(_stash); }

private:

  std::unique_ptr<T> obtain_stash() {
    std::unique_ptr<T> ptr;
#if HAVE_TBB
    if (not _data->_stash_queue.try_pop(ptr))
      ptr = _data->_factory();
#else
    auto guard = std::unique_lock{_data->_mutex};
    if (_data->_stash_queue.empty()) {
      ptr = _data->_factory();
    } else {
      ptr = std::move(_data->_stash_queue.front());
      _data->_stash_queue.pop();
    }
#endif
    return ptr;
  }

  void release_stash() noexcept {
    if(_stash) {
#if !HAVE_TBB
      auto guard = std::unique_lock{_data->_mutex};
#endif
      _data->_stash_queue.push(std::move(_stash));
    }
  }

  std::shared_ptr<Data> _data;
  std::unique_ptr<T> _stash;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_CONCURRENCY_SHARED_STASH_HH
