#ifndef DUNE_PDELAB_COMMON_LOCKMANAGER_HH
#define DUNE_PDELAB_COMMON_LOCKMANAGER_HH

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

  }
}
#endif // DUNE_PDELAB_COMMON_LOCKMANAGER_HH
