// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_FLATVECTORBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTL_FLATVECTORBACKEND_HH

#include <dune/common/fvector.hh>
#include <dune/istl/vector/host.hh>
#include <dune/typetree/typetree.hh>

#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/common/uncachedvectorview.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/backend/istl/vectorhelpers.hh>
#include <dune/pdelab/backend/istl/vectoriterator.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/backend/interface.hh>

namespace Dune {
  namespace PDELab {
    namespace istl {

    template<typename GFS, typename C>
    class FlatVectorContainer
      : public Backend::impl::Wrapper<C>
    {

      friend Backend::impl::Wrapper<C>;

    public:
      typedef typename C::value_type ElementType;
      typedef ElementType E;
      typedef C Container;
      typedef GFS GridFunctionSpace;
      typedef Container BaseT;
      typedef typename Container::value_type field_type;
      typedef typename Container::value_type value_type;
      typedef typename Container::size_type size_type;

      typedef typename GFS::Ordering::Traits::ContainerIndex ContainerIndex;

      typedef typename C::iterator iterator;
      typedef typename C::const_iterator const_iterator;

      template<typename LFSCache>
      using LocalView = UncachedVectorView<FlatVectorContainer,LFSCache>;

      template<typename LFSCache>
      using ConstLocalView = ConstUncachedVectorView<const FlatVectorContainer,LFSCache>;


      FlatVectorContainer(const FlatVectorContainer& rhs)
        : _gfs(rhs._gfs)
        , _container(std::make_shared<Container>(Backend::native(rhs)))
      {}

      FlatVectorContainer (const GFS& gfs, Dune::PDELab::Backend::attached_container = Dune::PDELab::Backend::attached_container())
        : _gfs(gfs)
        , _container(std::make_shared<Container>(gfs.ordering().blockCount()))
      {}

      //! Creates an FlatVectorContainer without allocating an underlying ISTL vector.
      FlatVectorContainer(const GFS& gfs, Dune::PDELab::Backend::unattached_container)
        : _gfs(gfs)
      {}

      FlatVectorContainer(const GFS& gfs, const E& e)
        : _gfs(gfs)
        , _container(std::make_shared<Container>(gfs.ordering().blockCount()))
      {
        (*_container)=e;
      }

      void detach()
      {
        _container.reset();
      }

      void attach(std::shared_ptr<Container> container)
      {
        _container = container;
      }

      bool attached() const
      {
        return bool(_container);
      }

      const std::shared_ptr<Container>& storage() const
      {
        return _container;
      }

      size_type N() const
      {
        return _container->N();
      }

      FlatVectorContainer& operator=(const FlatVectorContainer& r)
      {
        if (this == &r)
          return *this;
        if (attached())
          {
            (*_container) = Backend::native(r);
          }
        else
          {
            _container = std::make_shared<Container>(Backend::native(r));
          }
        return *this;
      }

      FlatVectorContainer& operator=(const E& e)
      {
        (*_container) = e;
        return *this;
      }

      FlatVectorContainer& operator*=(const E& e)
      {
        (*_container) *= e;
        return *this;
      }


      FlatVectorContainer& operator+=(const E& e)
      {
        (*_container) += e;
        return *this;
      }

      FlatVectorContainer& operator+=(const FlatVectorContainer& e)
      {
        (*_container) += Backend::native(e);
        return *this;
      }

      FlatVectorContainer& operator-=(const FlatVectorContainer& e)
      {
        (*_container) -= Backend::native(e);
        return *this;
      }

      E& block(std::size_t i)
      {
        return (*_container)[i];
      }

      const E& block(std::size_t i) const
      {
        return (*_container)[i];
      }

      E& operator[](const ContainerIndex& ci)
      {
        assert(ci.size() == 1);
        return block(ci[0]);
      }

      const E& operator[](const ContainerIndex& ci) const
      {
        assert(ci.size() == 1);
        return block(ci[0]);
      }

      typename Dune::template FieldTraits<E>::real_type two_norm() const
      {
        return _container->two_norm();
      }

      typename Dune::template FieldTraits<E>::real_type one_norm() const
      {
        return _container->one_norm();
      }

      typename Dune::template FieldTraits<E>::real_type infinity_norm() const
      {
        return _container->infinity_norm();
      }

      E operator*(const FlatVectorContainer& y) const
      {
        return (*_container) * Backend::native(y);
      }

      E dot(const FlatVectorContainer& y) const
      {
        return _container->dot(Backend::native(y));
      }

      FlatVectorContainer& axpy(const E& a, const FlatVectorContainer& y)
      {
        _container->axpy(a, Backend::native(y));
        return *this;
      }

      iterator begin()
      {
        return _container->begin();
      }

      const_iterator begin() const
      {
        return const_cast<const Container&>(*_container).begin();
      }

      iterator end()
      {
        return _container->end();
      }

      const_iterator end() const
      {
        return const_cast<const Container&>(*_container).end();
      }

      size_t flatsize() const
      {
        return _container->size();
      }

      const GFS& gridFunctionSpace() const
      {
        return _gfs;
      }

    private:

      Container& native()
      {
        return *_container;
      }

      const Container& native() const
      {
        return *_container;
      }

      const GFS& _gfs;
      std::shared_ptr<Container> _container;

    };


    } // namespace istl


#ifndef DOXYGEN

    namespace Backend {
      namespace impl {

        template<typename Allocator, typename GFS, typename E>
        struct BackendVectorSelectorHelper<istl::FlatVectorBackend<Allocator>, GFS, E>
        {

          typedef istl::FlatVectorContainer<
            GFS,
            Dune::ISTL::Vector<
              E,
              Allocator
              >
            > type;

          typedef type Type;

        };

      }
    }

#endif // DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_FLATVECTORBACKEND_HH
