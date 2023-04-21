#ifndef DUNE_PDELAB_CONCEPT_LOCKABLE_HH
#define DUNE_PDELAB_CONCEPT_LOCKABLE_HH

#include <concepts>


namespace Dune::PDELab::inline Experimental::Concept {

  template<class T>
  concept BasicLockable = requires(T m)
  {
    m.lock();
    m.unlock();
  };

  template<class T>
  concept Lockable = BasicLockable<T> && requires(T m)
  {
    { m.try_lock() } -> std::convertible_to<bool>;
  };

} // namespace Dune::PDELab::inline Experimental::Concept

#endif // DUNE_PDELAB_CONCEPT_LOCKABLE_HH
