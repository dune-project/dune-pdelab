// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_BLOCKVECTORBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTL_BLOCKVECTORBACKEND_HH

#include <dune/common/fvector.hh>
#include <dune/istl/blockvector/host.hh>
#include <dune/typetree/typetree.hh>

#include <dune/pdelab/backend/tags.hh>
#include <dune/pdelab/backend/common/uncachedvectorview.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/backend/istl/vectorhelpers.hh>
#include <dune/pdelab/backend/istl/vectoriterator.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/backend/backendselector.hh>

namespace Dune {
  namespace PDELab {
    namespace istl {

    template<typename GFS, typename C>
    class BlockVectorContainer
    {

    public:
      typedef typename C::value_type ElementType;
      typedef ElementType E;
      typedef E value_type;
      typedef C Container;
      typedef GFS GridFunctionSpace;
      typedef Container BaseT;
      typedef typename Container::value_type field_type;
      typedef typename Container::size_type size_type;

      typedef typename GFS::Ordering::Traits::ContainerIndex ContainerIndex;

      typedef typename C::iterator iterator;
      typedef typename C::const_iterator const_iterator;


#if HAVE_TEMPLATE_ALIASES

      template<typename LFSCache>
      using LocalView = UncachedVectorView<BlockVectorContainer,LFSCache>;

      template<typename LFSCache>
      using ConstLocalView = ConstUncachedVectorView<const BlockVectorContainer,LFSCache>;

#else

      template<typename LFSCache>
      struct LocalView
        : public UncachedVectorView<BlockVectorContainer,LFSCache>
      {

      LocalView()
      {}

      LocalView(BlockVectorContainer& vc)
        : UncachedVectorView<BlockVectorContainer,LFSCache>(vc)
      {}

    };

      template<typename LFSCache>
      struct ConstLocalView
        : public ConstUncachedVectorView<const BlockVectorContainer,LFSCache>
      {

      ConstLocalView()
      {}

      ConstLocalView(const BlockVectorContainer& vc)
        : BlockVectorView<const BlockVectorContainer,LFSCache>(vc)
      {}

    };

#endif // HAVE_TEMPLATE_ALIASES


      BlockVectorContainer(const BlockVectorContainer& rhs)
        : _gfs(rhs._gfs)
        , _container(make_shared<Container>(raw(rhs)))
      {}

      BlockVectorContainer (const GFS& gfs, Dune::PDELab::tags::attached_container = Dune::PDELab::tags::attached_container())
        : _gfs(gfs)
        , _container(make_shared<Container>(gfs.ordering().blockCount(),gfs.backend().blockSize()))
      {}

      //! Creates an BlockVectorContainer without allocating an underlying ISTL vector.
      BlockVectorContainer(const GFS& gfs, Dune::PDELab::tags::unattached_container)
        : _gfs(gfs)
      {}

      BlockVectorContainer(const GFS& gfs, const E& e)
        : _gfs(gfs)
        , _container(make_shared<Container>(gfs.ordering().blockCount(),gfs.backend().blockSize()))
      {
        (*_container)=e;
      }

      void detach()
      {
        _container.reset();
      }

      void attach(shared_ptr<Container> container)
      {
        _container = container;
      }

      bool attached() const
      {
        return bool(_container);
      }

      const shared_ptr<Container>& storage() const
      {
        return _container;
      }

      size_type N() const
      {
        return _container->N();
      }

      BlockVectorContainer& operator=(const BlockVectorContainer& r)
      {
        if (this == &r)
          return *this;
        if (attached())
          {
            (*_container) = raw(r);
          }
        else
          {
            _container = make_shared<Container>(raw(r));
          }
        return *this;
      }

      BlockVectorContainer& operator=(const E& e)
      {
        (*_container) = e;
        return *this;
      }

      BlockVectorContainer& operator*=(const E& e)
      {
        (*_container) *= e;
        return *this;
      }


      BlockVectorContainer& operator+=(const E& e)
      {
        (*_container) += e;
        return *this;
      }

      BlockVectorContainer& operator+=(const BlockVectorContainer& e)
      {
        (*_container) += raw(e);
        return *this;
      }

      BlockVectorContainer& operator-=(const BlockVectorContainer& e)
      {
        (*_container) -= raw(e);
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
        assert(ci.size() == 2);
        return (*_container)(ci[1],ci[0]);
      }

      const E& operator[](const ContainerIndex& ci) const
      {
        assert(ci.size() == 2);
        return (*_container)(ci[1],ci[0]);
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

      E operator*(const BlockVectorContainer& y) const
      {
        return (*_container) * raw(y);
      }

      E dot(const BlockVectorContainer& y) const
      {
        return _container->dot(raw(y));
      }

      BlockVectorContainer& axpy(const E& a, const BlockVectorContainer& y)
      {
        _container->axpy(a, raw(y));
        return *this;
      }

      // for debugging and AMG access
      Container& base ()
      {
        return *_container;
      }

      const Container& base () const
      {
        return *_container;
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
      const GFS& _gfs;
      shared_ptr<Container> _container;
    };


    } // namespace istl


#ifndef DOXYGEN

    template<typename Allocator, typename GFS, typename E>
    struct BackendVectorSelectorHelper<istl::BlockVectorBackend<Allocator>, GFS, E>
    {

      typedef istl::BlockVectorContainer<
        GFS,
        Dune::ISTL::BlockVector<
          E,
          Allocator
          >
        > type;

      typedef type Type;

    };

#endif // DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_BLOCKVECTORBACKEND_HH
