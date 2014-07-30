#ifndef DUNE_PDELAB_BACKEND_EIGEN_VECTOR_HH
#define DUNE_PDELAB_BACKEND_EIGEN_VECTOR_HH

#if HAVE_EIGEN

#include <dune/common/shared_ptr.hh>

#include <dune/istl/bvector.hh>

#include <dune/pdelab/backend/tags.hh>
#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/backend/common/uncachedvectorview.hh>
#include "descriptors.hh"
#include <Eigen/Dense>

namespace Dune {
  namespace PDELab {

    namespace EIGEN {

      template<typename GFS, typename ET>
      class VectorContainer
      {
      public:
        typedef Eigen::Matrix<ET, Eigen::Dynamic, 1> Container;
        typedef ET ElementType;
        typedef ET E;

        // for ISTL solver compatibility
        typedef ElementType field_type;

        typedef GFS GridFunctionSpace;
          typedef std::size_t size_type;

        typedef typename GFS::Ordering::Traits::ContainerIndex ContainerIndex;

        typedef ElementType* iterator;
        typedef const ElementType* const_iterator;
        // #warning iterators does not work well with Eigen
        // typedef typename Container::iterator iterator;
        // typedef typename Container::const_iterator const_iterator;

#if HAVE_TEMPLATE_ALIASES

        template<typename LFSCache>
        using LocalView = UncachedVectorView<VectorContainer,LFSCache>;

        template<typename LFSCache>
        using ConstLocalView = ConstUncachedVectorView<const VectorContainer,LFSCache>;

#else

        template<typename LFSCache>
        struct LocalView
          : public UncachedVectorView<VectorContainer,LFSCache>
        {

          LocalView()
          {}

          LocalView(VectorContainer& vc)
            : UncachedVectorView<VectorContainer,LFSCache>(vc)
          {}

        };

        template<typename LFSCache>
        struct ConstLocalView
          : public ConstUncachedVectorView<const VectorContainer,LFSCache>
        {

          ConstLocalView()
          {}

          ConstLocalView(const VectorContainer& vc)
            : ConstUncachedVectorView<const VectorContainer,LFSCache>(vc)
          {}

        };

#endif // HAVE_TEMPLATE_ALIASES

        VectorContainer(const VectorContainer& rhs)
          : _gfs(rhs._gfs)
          , _container(make_shared<Container>(*(rhs._container)))
        {}

        VectorContainer (const GFS& gfs, tags::attached_container = tags::attached_container())
          : _gfs(gfs)
          , _container(make_shared<Container>(gfs.ordering().blockCount()))
        {}

        //! Creates a VectorContainer without allocating storage.
        VectorContainer(const GFS& gfs, tags::unattached_container)
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
          , _container(make_shared<Container>(Container::Constant(gfs.ordering().blockCount(),e)))
        {}

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
              _container = make_shared<Container>(*(r._container));
            }
          return *this;
        }

        VectorContainer& operator=(const E& e)
        {
          (*_container) = Container::Constant(N(),e);
          return *this;
        }

        VectorContainer& operator*=(const E& e)
        {
          (*_container) *= e;
          return *this;
        }


        VectorContainer& operator+=(const E& e)
        {
          (*_container) += Container::Constant(N(),e);
          return *this;
        }

        VectorContainer& operator+=(const VectorContainer& y)
        {
          (*_container) += (*y._container);
          return *this;
        }

        VectorContainer& operator-= (const VectorContainer& y)
        {
          (*_container) -= (*y._container);
          return *this;
        }

        E& operator[](const ContainerIndex& ci)
        {
          return (*_container)(ci[0]);
        }

        const E& operator[](const ContainerIndex& ci) const
        {
          return (*_container)(ci[0]);
        }

        typename Dune::template FieldTraits<E>::real_type two_norm() const
        {
          return _container->norm();
        }

        typename Dune::template FieldTraits<E>::real_type one_norm() const
        {
          return _container->template lpNorm<1>();
        }

        typename Dune::template FieldTraits<E>::real_type infinity_norm() const
        {
          return _container->template lpNorm<Eigen::Infinity>();
        }

        //! (*this)^T y
        E operator*(const VectorContainer& y) const
        {
          return _container->transpose() * (*y._container);
        }

        //! (*this)^+ y
        E dot(const VectorContainer& y) const
        {
          return _container->dot(*y._container);
        }

        //! vector space axpy operation ( *this += a y )
        VectorContainer& axpy(const E& a, const VectorContainer& y)
        {
          (*_container) += a * (*y._container);
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
          return _container->data();
        }

        const_iterator begin() const
        {
          return _container->data();
        }

        iterator end()
        {
          return _container->data() + N();
        }

        const_iterator end() const
        {
          return _container->data() + N();
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

    } // end namespace EIGEN


#ifndef DOXYGEN

    template<typename GFS, typename E>
    struct EigenVectorSelectorHelper
    {
      using Type = EIGEN::VectorContainer<GFS, E>;
    };

    template<typename GFS, typename E>
    struct BackendVectorSelectorHelper<EigenVectorBackend, GFS, E>
      : public EigenVectorSelectorHelper<GFS,E>
    {};

#endif // DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif

#endif // DUNE_PDELAB_BACKEND_EIGEN_VECTOR_HH

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
