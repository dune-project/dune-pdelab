#ifndef DUNE_PDELAB_COMMON_LOCKMANAGER_HH
#define DUNE_PDELAB_COMMON_LOCKMANAGER_HH

#include <cstddef>

#include <dune/pdelab/common/elementmapper.hh>

namespace Dune{
  namespace PDELab{

    //! LockManager that uses one global lock.
    template<class Mutex>
    class GlobalLockManager
    {
      Mutex mutex_;

    public:
      typedef Mutex value_type;

      template<class Entity>
      Mutex & operator[](const Entity &e)
      {
        return mutex_;
      }
    };

    //! Hold a vector of mutexes
    /**
     * Mutexes usually can only be default-constructed, but cannot be copied
     * or moved, neither by construction nor by assignment.  This vector can
     * hold a variable number of such mutexes.  The vector can be resized, but
     * that reconstructs all mutexes, so references to old mutexes are no
     * longer valid.
     */
    template<class Mutex>
    class MutexVector
    {
      Mutex* data_;
      std::size_t size_;

    public:
      typedef Mutex value_type;
      typedef Mutex* iterator;

      MutexVector(std::size_t size = 0) :
        data_(0), size_(0)
      {
        resize(size);
      }
      ~MutexVector()
      {
        clear();
      }

      void clear() {
        if(data_)
          delete[] data_;
        data_ = 0;
        size_ = 0;
      }

      void resize(std::size_t newsize) {
        clear();
        if(newsize) {
          data_ = new Mutex[newsize];
          size_ = newsize;
        }
      }

      Mutex &operator[](std::size_t i)
      {
        return data_[i];
      }

      std::size_t size() const
      {
        return size_;
      }

      iterator begin()
      {
        return data_;
      }
      iterator end()
      {
        return data_ + size_;
      }

    };

    //! LockManager that uses one lock per Element.
    template<class GridView, class Mutex>
    class PerElementLockManager :
      public MutexVector<Mutex>
    {
      ElementMapper<GridView> emapper;

    public:
      PerElementLockManager(const GridView &gv) :
        emapper(gv)
      {
        this->resize(emapper.size());
      }

      Mutex & operator[](const typename GridView::template Codim<0>::Entity &e)
      {
        return MutexVector<Mutex>::operator[](emapper.map(e));
      }
    };

  }
}
#endif // DUNE_PDELAB_COMMON_LOCKMANAGER_HH
