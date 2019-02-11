#ifndef DUNE_PDELAB_COMMON_DESTRUCTIBLESINGLETONHOLDER_HH
#define DUNE_PDELAB_COMMON_DESTRUCTIBLESINGLETONHOLDER_HH

#include <memory>
#include <dune/pdelab/common/exceptions.hh>

namespace Dune::PDELab {


  template<typename Singleton, typename Factory>
  class DestructibleSingletonHolder
  {

    std::unique_ptr<Singleton> _data;
    Factory _factory;

  public:

    template<typename... T>
    void create(T&... args)
    {
      if (_data)
        DUNE_THROW(Exception,"Singleton already created");
      _data = _factory(std::forward<T>(args)...);
    }

    operator bool() const
    {
      return bool(_data);
    }

    Singleton& get()
    {
      return *_data;
    }

    void destroy()
    {
      if (not _data)
        DUNE_THROW(Exception,"Singleton not present");
      _data.reset();
    }

    DestructibleSingletonHolder(Factory factory)
      : _factory(std::move(factory))
    {}

  };

  template<typename Singleton, typename Factory>
  DestructibleSingletonHolder<Singleton,Factory>& destructibleSingleton(Factory factory)
  {
    static DestructibleSingletonHolder<Singleton,Factory> holder(std::move(factory));
    return holder;
  }

} // Dune::PDELab


#endif // DUNE_PDELAB_COMMON_DESTRUCTIBLESINGLETONHOLDER_HH
