#ifndef DUNE_PDELAB_COMMON_LOCKMANAGER_HH
#define DUNE_PDELAB_COMMON_LOCKMANAGER_HH

#include <cstddef>

#include <dune/common/documentation.hh>

#include <dune/pdelab/common/elementmapper.hh>

namespace Dune{
  namespace PDELab{

    //! Interface for lock managers
    class LockManagerInterface {
    public:
      //! type of the mutex as returned by operator[]()
      typedef ImplementationDefined value_type;

      //! obtain mutex for the given entity
      /**
       * Multiple entities may share a mutex.  When a thread must lock
       * multiple entities, and it is not known that the mutexes are
       * recursive, the thread must ensure not to try and lock the same mutex
       * twice.  The thread can find out whether two mutexes are the same by
       * comparing their adresses.
       *
       * \note The LockManager is not required to implement this function as a
       *       template functions.  It just must be possible to call this
       *       function with any entity the LockManager supports as an
       *       argument.  An alternate valid implementation would be via
       *       overloading.
       */
      template<class Entity>
      value_type &operator[](const Entity &);
    };

    //! LockManager that uses one global lock.
    /**
     * \implements LockManagerInterface
     */
    template<class Mutex>
    class GlobalLockManager
    {
      Mutex mutex_;

    public:
      //! type of the mutex
      typedef Mutex value_type;

      //! obtain mutex for the given entity
      /**
       * \note This lock manager will always return the same mutex.
       */
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
      //! type of the mutexes
      typedef Mutex value_type;
      //! iterator type
      typedef Mutex* iterator;

      //! construct
      /**
       * \param size Number of mutexes held.
       */
      MutexVector(std::size_t size = 0) :
        data_(0), size_(0)
      {
        resize(size);
      }
      //! destroy
      /**
       * This will call clear()
       */
      ~MutexVector()
      {
        clear();
      }

      //! free the storage for all mutexes
      /**
       * Precodition: none of the mutexes are locked.
       *
       * Invalidates any iterators.
       */
      void clear() {
        if(data_)
          delete[] data_;
        data_ = 0;
        size_ = 0;
      }

      //! resize the vector
      /**
       * This frees the old vector and reallocates a new vector of the correct
       * size.
       *
       * Precodition: none of the mutexes are locked.
       *
       * Postcondition: none of the mutexes is locked.
       *
       * Invalidates any iterators.
       */
      void resize(std::size_t newsize) {
        clear();
        if(newsize) {
          data_ = new Mutex[newsize];
          size_ = newsize;
        }
      }

      //! return a reference to the i'th mutex
      Mutex &operator[](std::size_t i)
      {
        return data_[i];
      }

      //! number of mutexes currently in the vector
      std::size_t size() const
      {
        return size_;
      }

      //! iterator to the first mutex
      iterator begin()
      {
        return data_;
      }
      //! passed-the-end iterator
      iterator end()
      {
        return data_ + size_;
      }

    };

    //! LockManager that uses one lock per Element.
    /**
     * \implements LockManagerInterface
     */
    template<class GridView, class Mutex>
    class PerElementLockManager :
      public MutexVector<Mutex>
    {
      ElementMapper<GridView> emapper;

    public:
      //! construct
      /**
       * The gridview is needed to construct a mapper of the elements.
       */
      PerElementLockManager(const GridView &gv) :
        emapper(gv)
      {
        this->resize(emapper.size());
      }

      //! obtain mutex for the given entity
      /**
       * \note This lock manager will return different mutexes for different
       *       elements.  However, a user of this LockManager must still take
       *       care not to lock a mutex twice e.g. when assembling diagonal
       *       entries of a jacobian.  See the hints in the documantation of
       *       LockManagerInterface::operator[]().
       */
      Mutex & operator[](const typename GridView::template Codim<0>::Entity &e)
      {
        return MutexVector<Mutex>::operator[](emapper.map(e));
      }
    };

  }
}
#endif // DUNE_PDELAB_COMMON_LOCKMANAGER_HH
