// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_SIMPLE_VECTOR_HH
#define DUNE_PDELAB_BACKEND_SIMPLE_VECTOR_HH

#include <algorithm>
#include <functional>
#include <memory>
#include <numeric>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/common/uncachedvectorview.hh>
#include <dune/pdelab/backend/simple/descriptors.hh>

namespace Dune {
  namespace PDELab {

    namespace Simple {

      namespace {

        // For some reason std::bind cannot directly deduce the correct version
        // of Dune::fvmeta::abs2, so we package it into a functor to help it along.
        template<typename K>
        struct abs2
        {
          typename FieldTraits<K>::real_type operator()(const K& k) const
          {
            return Dune::fvmeta::abs2(k);
          }
        };

      }

      template<typename GFS, typename C>
      class VectorContainer
        : public Backend::impl::Wrapper<C>
      {

        friend Backend::impl::Wrapper<C>;

      public:
        typedef C Container;
        typedef typename Container::value_type ElementType;
        typedef ElementType E;

        // for ISTL solver compatibility
        typedef ElementType field_type;

        typedef GFS GridFunctionSpace;
        typedef typename Container::size_type size_type;

        typedef typename GFS::Ordering::Traits::ContainerIndex ContainerIndex;

        typedef typename Container::iterator iterator;
        typedef typename Container::const_iterator const_iterator;

        template<typename LFSCache>
        using LocalView = UncachedVectorView<VectorContainer,LFSCache>;

        template<typename LFSCache>
        using ConstLocalView = ConstUncachedVectorView<const VectorContainer,LFSCache>;


        VectorContainer(const VectorContainer& rhs)
          : _gfs(rhs._gfs)
          , _container(std::make_shared<Container>(*(rhs._container)))
        {}

        VectorContainer (const GFS& gfs, Backend::attached_container = Backend::attached_container())
          : _gfs(gfs)
          , _container(std::make_shared<Container>(gfs.ordering().blockCount()))
        {}

        //! Creates a VectorContainer without allocating storage.
        VectorContainer(const GFS& gfs, Backend::unattached_container)
          : _gfs(gfs)
        {}

        /** \brief Constructs an VectorContainer for an explicitly given vector object
         *
         * \param gfs GridFunctionSpace that determines the size and the blocking of the vector
         * \param container The actual container class
         */
        VectorContainer (const GFS& gfs, Container& container)
          : _gfs(gfs)
          , _container(stackobject_to_shared_ptr(container))
        {
          _container->resize(gfs.ordering().blockCount());
        }

        VectorContainer (const GFS& gfs, const E& e)
          : _gfs(gfs)
          , _container(std::make_shared<Container>(gfs.ordering().blockCount(),e))
        {}

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
          return _container->size();
        }

        VectorContainer& operator=(const VectorContainer& r)
        {
          if (this == &r)
            return *this;
          if (attached())
            {
              (*_container) = (*r._container);
            }
          else
            {
              _container = std::make_shared<Container>(*(r._container));
            }
          return *this;
        }

        VectorContainer& operator=(const E& e)
        {
          std::fill(_container->begin(),_container->end(),e);
          return *this;
        }

        VectorContainer& operator*=(const E& e)
        {
          std::transform(_container->begin(),_container->end(),_container->begin(),
                         std::bind(std::multiplies<E>(),e,std::placeholders::_1));
          return *this;
        }


        VectorContainer& operator+=(const E& e)
        {
          std::transform(_container->begin(),_container->end(),_container->begin(),
                         std::bind(std::plus<E>(),e,std::placeholders::_1));
          return *this;
        }

        VectorContainer& operator+=(const VectorContainer& y)
        {
          std::transform(_container->begin(),_container->end(),y._container->begin(),
                         _container->begin(),std::plus<E>());
          return *this;
        }

        VectorContainer& operator-= (const VectorContainer& y)
        {
          std::transform(_container->begin(),_container->end(),y._container->begin(),
                         _container->begin(),std::minus<E>());
          return *this;
        }

        E& operator[](const ContainerIndex& ci)
        {
          return (*_container)[ci[0]];
        }

        const E& operator[](const ContainerIndex& ci) const
        {
          return (*_container)[ci[0]];
        }

        typename Dune::template FieldTraits<E>::real_type two_norm() const
        {
          using namespace std::placeholders;
          typedef typename Dune::template FieldTraits<E>::real_type Real;
          return std::sqrt(std::accumulate(_container->begin(),_container->end(),Real(0),std::bind(std::plus<Real>(),_1,std::bind(abs2<E>(),_2))));
        }

        typename Dune::template FieldTraits<E>::real_type one_norm() const
        {
          using namespace std::placeholders;
          typedef typename Dune::template FieldTraits<E>::real_type Real;
          return std::accumulate(_container->begin(),_container->end(),Real(0),std::bind(std::plus<Real>(),_1,std::bind(std::abs<E>,_2)));
        }

        typename Dune::template FieldTraits<E>::real_type infinity_norm() const
        {
          if (_container->size() == 0)
            return 0;
          using namespace std::placeholders;
          typedef typename Dune::template FieldTraits<E>::real_type Real;
          return *std::max_element(_container->begin(),_container->end(),std::bind(std::less<Real>(),std::bind(std::abs<E>,_1),std::bind(std::abs<E>,_2)));
        }

        E operator*(const VectorContainer& y) const
        {
          return std::inner_product(_container->begin(),_container->end(),y._container->begin(),E(0));
        }

        E dot(const VectorContainer& y) const
        {
          return std::inner_product(_container->begin(),_container->end(),y._container->begin(),E(0),std::plus<E>(),Dune::dot<E,E>);
        }

        VectorContainer& axpy(const E& a, const VectorContainer& y)
        {
          using namespace std::placeholders;
          std::transform(_container->begin(),_container->end(),y._container->begin(),
                         _container->begin(),std::bind(std::plus<E>(),_1,std::bind(std::multiplies<E>(),a,_2)));
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

      private:

        Container& native ()
        {
          return *_container;
        }

        const Container& native () const
        {
          return *_container;
        }

      public:

        iterator begin()
        {
          return _container->begin();
        }

        const_iterator begin() const
        {
          return _container->begin();
        }

        iterator end()
        {
          return _container->end();
        }

        const_iterator end() const
        {
          return _container->end();
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
        std::shared_ptr<Container> _container;
      };


    } // namespace Simple


#ifndef DOXYGEN

    template<typename GFS, typename E>
    struct SimpleVectorSelectorHelper
    {

      using vector_type = typename GFS::Traits::Backend::template vector_type<E>;

      using Type = Simple::VectorContainer<GFS,vector_type>;

    };

    namespace Backend {
      namespace impl {

        template<template<typename> class Container, typename GFS, typename E>
        struct BackendVectorSelectorHelper<Simple::VectorBackend<Container>, GFS, E>
          : public SimpleVectorSelectorHelper<GFS,E>
        {};

      } // namespace impl
    } // namespace Backend

#endif // DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_SIMPLE_VECTOR_HH
