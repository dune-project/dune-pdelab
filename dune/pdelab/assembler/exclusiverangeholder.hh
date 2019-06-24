// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_EXCLUSIVERANGEHOLDER_HH
#define DUNE_PDELAB_ASSEMBLER_EXCLUSIVERANGEHOLDER_HH

#include <cassert>
#include <utility>

#include <dune/pdelab/common/checks.hh>
#include <dune/pdelab/common/typetraits.hh>
#include <dune/pdelab/common/exceptions.hh>


namespace Dune::PDELab::Experimental {

  template<typename Implementation, typename Container>
  class ExclusiveRangeHolder
  {

  public:

    class Proxy {

    public:

      using value_type = const typename Container::value_type;
      using const_reference = typename Container::const_reference;
      using reference = const_reference;
      using const_iterator = typename Container::const_iterator;
      using iterator   = const_iterator;
      using size_type      = typename Container::size_type;

      size_type size() const
      {
        assert(_range_holder);
        return _range_holder->accessContainer().size();
      }

      iterator begin() const
      {
        assert(_range_holder);
        return _range_holder->accessContainer().begin();
      }

      iterator end() const
      {
        assert(_range_holder);
        return _range_holder->accessContainer().end();
      }

      const_reference operator[](size_type i) const
      {
        assert(_range_holder);
        return _range_holder->accessContainer()[i];
      }

      Proxy() noexcept
      : _range_holder(nullptr)
      {}

      Proxy(ExclusiveRangeHolder& range_holder)
        : _range_holder(&range_holder)
      {
#if DUNE_PDELAB_DEBUG_RANGE_PROXY
        ++_range_holder->_proxy_count;
#endif
      }

#if DUNE_PDELAB_DEBUG_RANGE_PROXY
      ~Proxy() noexcept
      {
        if (_range_holder)
          --_range_holder->_proxy_count;
      }

      Proxy(const Proxy& proxy) noexcept
        : _range_holder(proxy._range_holder)
      {
        if (_range_holder)
          ++_range_holder->_proxy_count;
      }

      friend void swap(Proxy& a, Proxy& b) noexcept
      {
        using std::swap;
        swap(a._range_holder,b._range_holder);
      }

      Proxy(Proxy&& proxy) noexcept
        : Proxy()
      {
        swap(*this,proxy);
      }

      Proxy& operator=(const Proxy& proxy) noexcept
      {
        Proxy tmp(proxy);
        swap(*this,tmp);
        return *this;
      }

      Proxy& operator=(Proxy&& proxy) noexcept
      {
        swap(*this,proxy);
        return *this;
      }
#endif

    private:

      const ExclusiveRangeHolder* _range_holder;

    };

    Proxy proxy()
    {
      return Proxy{*this};
    }

    Container& container()
    {
#if DUNE_PDELAB_DEBUG_RANGE_PROXY
      if (_proxy_count > 0)
        DUNE_THROW(Dune::Exception,"Attempt to access underlying container while there are still " << _proxy_count << " proxies around.");
#endif
      return accessContainer();
    }

    void assertWriteAccess() const
    {
#if DUNE_PDELAB_DEBUG_RANGE_PROXY
      if (_proxy_count > 0)
        DUNE_THROW(Dune::Exception,"Attempt to access underlying container while there are still " << _proxy_count << " proxies around.");
#endif
    }

  private:

    Container& accessContainer()
    {
      return static_cast<Implementation*>(this)->accessContainer();
    }

    const Container& accessContainer() const
    {
      return static_cast<const Implementation*>(this)->accessContainer();
    }

#if DUNE_PDELAB_DEBUG_RANGE_PROXY
    mutable std::size_t _proxy_count;
#endif

  };

} // namespace Dune::PDELab::Experimental

#endif // DUNE_PDELAB_ASSEMBLER_EXCLUSIVERANGEHOLDER_HH
