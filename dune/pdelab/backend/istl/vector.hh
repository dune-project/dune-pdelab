// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_VECTOR_HH
#define DUNE_PDELAB_BACKEND_ISTL_VECTOR_HH

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/typetree/typetree.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/common/uncachedvectorview.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/backend/istl/vectorhelpers.hh>
#include <dune/pdelab/backend/istl/vectoriterator.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

namespace Dune {
  namespace PDELab {
    namespace istl {

      template<typename GFS, typename C>
      class BlockVector
        : public Backend::impl::Wrapper<C>
      {

        friend Backend::impl::Wrapper<C>;

      public:
        typedef typename C::field_type ElementType;
        typedef ElementType E;
        typedef C Container;
        typedef GFS GridFunctionSpace;
        typedef Container BaseT;
        typedef typename Container::field_type field_type;
        typedef typename Container::block_type block_type;
        typedef typename Container::size_type size_type;

        typedef typename GFS::Ordering::Traits::ContainerIndex ContainerIndex;

        typedef istl::vector_iterator<C> iterator;
        typedef istl::vector_iterator<const C> const_iterator;


#if HAVE_TEMPLATE_ALIASES

        template<typename LFSCache>
        using LocalView = UncachedVectorView<BlockVector,LFSCache>;

        template<typename LFSCache>
        using ConstLocalView = ConstUncachedVectorView<const BlockVector,LFSCache>;

#else

        template<typename LFSCache>
        struct LocalView
          : public UncachedVectorView<BlockVector,LFSCache>
        {

        LocalView()
        {}

        LocalView(BlockVector& vc)
          : UncachedVectorView<BlockVector,LFSCache>(vc)
        {}

      };

        template<typename LFSCache>
        struct ConstLocalView
          : public ConstUncachedVectorView<const BlockVector,LFSCache>
        {

        ConstLocalView()
        {}

        ConstLocalView(const BlockVector& vc)
          : ConstUncachedVectorView<const BlockVector,LFSCache>(vc)
        {}

      };

#endif // HAVE_TEMPLATE_ALIASES


        BlockVector(const BlockVector& rhs)
          : _gfs(rhs._gfs)
          , _container(std::make_shared<Container>(_gfs.ordering().blockCount()))
        {
          istl::dispatch_vector_allocation(_gfs.ordering(),*_container,typename GFS::Ordering::ContainerAllocationTag());
          (*_container) = rhs.base();
        }

        BlockVector (const GFS& gfs, Backend::attached_container = Backend::attached_container())
          : _gfs(gfs)
          , _container(std::make_shared<Container>(gfs.ordering().blockCount()))
        {
          istl::dispatch_vector_allocation(gfs.ordering(),*_container,typename GFS::Ordering::ContainerAllocationTag());
        }

        //! Creates an BlockVector without allocating an underlying ISTL vector.
        BlockVector(const GFS& gfs, Backend::unattached_container)
          : _gfs(gfs)
        {}

        /** \brief Constructs an BlockVector for an explicitly given vector object
         *
         * \param gfs GridFunctionSpace that determines the size and the blocking of the vector
         * \param container The actual ISTL container class
         */
        BlockVector (const GFS& gfs, Container& container)
          : _gfs(gfs)
          , _container(stackobject_to_shared_ptr(container))
        {
          _container->resize(gfs.ordering().blockCount());
          istl::dispatch_vector_allocation(gfs.ordering(),*_container,typename GFS::Ordering::ContainerAllocationTag());
        }

        BlockVector (const GFS& gfs, const E& e)
          : _gfs(gfs)
          , _container(std::make_shared<Container>(gfs.ordering().blockCount()))
        {
          istl::dispatch_vector_allocation(gfs.ordering(),*_container,typename GFS::Ordering::ContainerAllocationTag());
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

        BlockVector& operator= (const BlockVector& r)
        {
          if (this == &r)
            return *this;
          if (attached())
            {
              (*_container) = r.base();
            }
          else
            {
              _container = std::make_shared<Container>(r.base());
            }
          return *this;
        }

        BlockVector& operator= (const E& e)
        {
          (*_container)=e;
          return *this;
        }

        BlockVector& operator*= (const E& e)
        {
          (*_container)*=e;
          return *this;
        }


        BlockVector& operator+= (const E& e)
        {
          (*_container)+=e;
          return *this;
        }

        BlockVector& operator+= (const BlockVector& e)
        {
          (*_container)+= e.base();
          return *this;
        }

        BlockVector& operator-= (const BlockVector& e)
        {
          (*_container)-= e.base();
          return *this;
        }

        block_type& block(std::size_t i)
        {
          return (*_container)[i];
        }

        const block_type& block(std::size_t i) const
        {
          return (*_container)[i];
        }

        E& operator[](const ContainerIndex& ci)
        {
          return istl::access_vector_element(istl::container_tag(*_container),*_container,ci,ci.size()-1);
        }

        const E& operator[](const ContainerIndex& ci) const
        {
          return istl::access_vector_element(istl::container_tag(*_container),*_container,ci,ci.size()-1);
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

        E operator*(const BlockVector& y) const
        {
          return (*_container)*y.base();
        }

        E dot(const BlockVector& y) const
        {
          return _container->dot(y.base());
        }

        BlockVector& axpy(const E& a, const BlockVector& y)
        {
          _container->axpy(a, y.base());
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

        // for debugging and AMG access
        Container& native ()
        {
          return *_container;
        }

        const Container& native () const
        {
          return *_container;
        }

      public:

        operator Container&()
        {
          return *_container;
        }

        operator const Container&() const
        {
          return *_container;
        }

        iterator begin()
        {
          return iterator(*_container,false);
        }


        const_iterator begin() const
        {
          return const_iterator(*_container,false);
        }

        iterator end()
        {
          return iterator(*_container,true);
        }


        const_iterator end() const
        {
          return const_iterator(*_container,true);
        }

        size_t flatsize() const
        {
          return _container->dim();
        }

        const GFS& gridFunctionSpace() const
        {
          return _gfs;
        }

      private:
        const GFS& _gfs;
        std::shared_ptr<Container> _container;
      };

#ifndef DOXYGEN

      // helper struct invoking the GFS tree -> ISTL vector reduction
      template<typename GFS, typename E>
      struct BlockVectorSelectorHelper
      {

        typedef typename TypeTree::AccumulateType<
          GFS,
          istl::vector_creation_policy<E>
          >::type vector_descriptor;

        typedef BlockVector<GFS,typename vector_descriptor::vector_type> Type;

      };

#endif // DOXYGEN

    // can't have the closing of the namespace inside the #ifndef DOXYGEN block
    } // namespace istl

#ifndef DOXYGEN

    namespace Backend {
      namespace impl {

        template<Dune::PDELab::istl::Blocking blocking, std::size_t block_size, typename GFS, typename E>
        struct BackendVectorSelectorHelper<istl::VectorBackend<blocking,block_size>, GFS, E>
          : public istl::BlockVectorSelectorHelper<GFS,E>
        {};

      } // namespace impl
    } // namespace Backend

#endif // DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_VECTOR_HH
