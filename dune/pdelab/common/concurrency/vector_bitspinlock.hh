#ifndef DUNE_PDELAB_COMMON_CONCURRENCY_VECTOR_BITSPINLOCK_HH
#define DUNE_PDELAB_COMMON_CONCURRENCY_VECTOR_BITSPINLOCK_HH

#include <vector>
#include <atomic>
#include <climits>
#include <cassert>

namespace Dune::PDELab::inline Experimental {

/**
 * @brief Vector of bit spin locks
 * Each bit in the vector is a different lock.
 *
 * @warning Notice that locks don't respect the cache line sizes, so
 * heterogeneous thread access of neighboring locks will result in false sharing.
 */
class VectorBitSpinLock
{
  using T = unsigned char;
  static constexpr std::size_t T_BIT =  CHAR_BIT*sizeof(T);
  static constexpr T T_NULL =  T{} ^ T{};
  static_assert(std::atomic<T>::is_always_lock_free);

public:
  /**
   * @brief Handle of an bit spin lock
   * @details Manupulates the acquisition and release of the lock and fullfils
   * the standard Lockable requirement
   */
  struct LockHandle {
    LockHandle(std::size_t i, std::atomic<T>& chunk)
      : _lock{static_cast<T>(1 << i%T_BIT)}
      , _chunk{chunk}
    {}

    LockHandle(const LockHandle&) = delete;

    /**
     * @brief Acquire ownership of the lock
     * If locked, this lock will spin until unlocked by another thread
     */
     void lock() {
      if (try_lock()) [[likely]]
        return;
      do {
        wait();
      } while (not try_lock());
    }

    //! Try to acquire ownership of the lock. Return true in success
     [[nodiscard]] bool try_lock() noexcept {
      return not (_chunk.fetch_or(_lock, std::memory_order_acquire) & _lock);
    }

    //! Relesase ownership of the lock
     void unlock() noexcept {
      _chunk.fetch_and(~_lock, std::memory_order_release);
    }

    //! Wait until lock could be acquired
     void wait() const noexcept {
      do {
        pause();
      } while (locked());
    }

    [[nodiscard]] bool operator==(const LockHandle& other) const {
      return (&_chunk == &other._chunk) and (_lock == other._lock);
    }

  private:
    [[nodiscard]] bool locked() const noexcept {
      return _chunk.load(std::memory_order_relaxed) & _lock;
    }

    inline void pause() const noexcept {
#if defined(__x86_64__) || defined(__i386__)
      __asm__ __volatile__ ("pause");
#endif
    }

    const T _lock;
    std::atomic<T>& _chunk;
  };

  //! Constructs a vector with a given size
  VectorBitSpinLock(std::size_t size = 0) {
    resize(size);
  }

  //! Get the size of the vector
  std::size_t size() const
  {
    return _size;
  }

  //! Access the i-th lock
  LockHandle operator[](std::size_t i) {
    assert(i < _size);
    assert(i/T_BIT < _data.size());
    return {i, _data[i/T_BIT]};
  }

  //! Sets new size of the lock vector. All previous locks are dumped
  void resize(std::size_t new_size) {
    _size = new_size;
    std::size_t chunks = (_size+T_BIT-1)/T_BIT;
    _data = std::vector<std::atomic<T>>(chunks);

    for (auto& v : _data)
      v = T_NULL;
  }

private:
  std::vector<std::atomic<T>> _data;
  std::size_t _size;
};

} //namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_CONCURRENCY_VECTOR_BITSPINLOCK_HH
