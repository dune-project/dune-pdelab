// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTLVECTORBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTLVECTORBACKEND_HH

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/backend/tags.hh>
#include <dune/pdelab/backend/istl/vectorhelpers.hh>
#include <dune/pdelab/backend/istl/vectoriterator.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/backend/backendselector.hh>

namespace Dune {
  namespace PDELab {

    template<typename GFS, typename C>
    class ISTLBlockVectorContainer
    {

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

      template<typename LFSCache>
      struct LocalView
      {

        //dune_static_assert((is_same<typename LFSCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFS>::value),
        //                   "The LocalFunctionSpace passed to LocalView must belong to the underlying GridFunctionSpace.");

        typedef E ElementType;
        //typedef typename LFSCache::LocalFunctionSpace LFS;
        typedef typename LFSCache::DOFIndex DOFIndex;
        typedef typename LFSCache::ContainerIndex ContainerIndex;

        LocalView()
          : _container(nullptr)
          , _lfs_cache(nullptr)
        {}

        LocalView(ISTLBlockVectorContainer& container)
          : _container(&container)
          , _lfs_cache(nullptr)
        {}

        void attach(ISTLBlockVectorContainer& container)
        {
          _container = &container;
        }

        void detach()
        {
          _container = nullptr;
        }

        void bind(const LFSCache& lfs_cache)
        {
          _lfs_cache = &lfs_cache;
        }

        void unbind()
        {
        }

        size_type size() const
        {
          return _lfs_cache->size();
        }

        template<typename LC>
        void read(LC& local_container) const
        {
          for (size_type i = 0; i < size(); ++i)
            {
              accessBaseContainer(local_container)[i] = (*_container)[_lfs_cache->containerIndex(i)];
            }
        }

        template<typename LC>
        void write(const LC& local_container)
        {
          for (size_type i = 0; i < size(); ++i)
            {
              (*_container)[_lfs_cache->containerIndex(i)] = accessBaseContainer(local_container)[i];
            }
        }

        template<typename LC>
        void add(const LC& local_container)
        {
          for (size_type i = 0; i < size(); ++i)
            {
              (*_container)[_lfs_cache->containerIndex(i)] += accessBaseContainer(local_container)[i];
            }
        }

        template<typename ChildLFS, typename LC>
        void read(const ChildLFS& child_lfs, LC& local_container) const
        {
          for (size_type i = 0; i < child_lfs.size(); ++i)
            {
              const size_type local_index = child_lfs.localIndex(i);
              accessBaseContainer(local_container)[local_index] = (*_container)[_lfs_cache->containerIndex(local_index)];
            }
        }

        template<typename ChildLFS, typename LC>
        void write(const ChildLFS& child_lfs, const LC& local_container)
        {
          for (size_type i = 0; i < child_lfs.size(); ++i)
            {
              const size_type local_index = child_lfs.localIndex(i);
              (*_container)[_lfs_cache->containerIndex(local_index)] = accessBaseContainer(local_container)[local_index];
            }
        }

        template<typename ChildLFS, typename LC>
        void add(const ChildLFS& child_lfs, const LC& local_container)
        {
          for (size_type i = 0; i < child_lfs.size(); ++i)
            {
              const size_type local_index = child_lfs.localIndex(i);
              (*_container)[_lfs_cache->containerIndex(local_index)] += accessBaseContainer(local_container)[local_index];
            }
        }


        template<typename ChildLFS, typename LC>
        void read_sub_container(const ChildLFS& child_lfs, LC& local_container) const
        {
          for (size_type i = 0; i < child_lfs.size(); ++i)
            {
              const size_type local_index = child_lfs.localIndex(i);
              accessBaseContainer(local_container)[i] = (*_container)[_lfs_cache->containerIndex(local_index)];
            }
        }

        template<typename ChildLFS, typename LC>
        void write_sub_container(const ChildLFS& child_lfs, const LC& local_container)
        {
          for (size_type i = 0; i < child_lfs.size(); ++i)
            {
              const size_type local_index = child_lfs.localIndex(i);
              (*_container)[_lfs_cache->containerIndex(local_index)] = accessBaseContainer(local_container)[i];
            }
        }

        template<typename ChildLFS, typename LC>
        void add_sub_container(const ChildLFS& child_lfs, const LC& local_container)
        {
          for (size_type i = 0; i < child_lfs.size(); ++i)
            {
              const size_type local_index = child_lfs.localIndex(i);
              (*_container)[_lfs_cache->containerIndex(local_index)] += accessBaseContainer(local_container)[i];
            }
        }

        void commit()
        {
        }


        ElementType& operator[](size_type i)
        {
          return (*_container)[_lfs_cache->containerIndex(i)];
        }

        const ElementType& operator[](size_type i) const
        {
          return (*_container)[_lfs_cache->containerIndex(i)];
        }

        ElementType& operator[](const DOFIndex& di)
        {
          return (*_container)[_lfs_cache->containerIndex(di)];
        }

        const ElementType& operator[](const DOFIndex& di) const
        {
          return (*_container)[_lfs_cache->containerIndex(di)];
        }

        ElementType& operator[](const ContainerIndex& ci)
        {
          return (*_container)[ci];
        }

        const ElementType& operator[](const ContainerIndex& ci) const
        {
          return (*_container)[ci];
        }

        ISTLBlockVectorContainer& global_container()
        {
          return *_container;
        }

        const ISTLBlockVectorContainer& global_container() const
        {
          return *_container;
        }

        const LFSCache& cache() const
        {
          return *_lfs_cache;
        }

      private:

        ISTLBlockVectorContainer* _container;
        const LFSCache* _lfs_cache;

      };


      template<typename LFSCache>
      struct ConstLocalView
      {

        // dune_static_assert((is_same<typename LFSCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFS>::value),
        //                    "The LocalFunctionSpace passed to LocalView must belong to the underlying GridFunctionSpace.");

        typedef E ElementType;
        typedef typename LFSCache::LocalFunctionSpace LFS;
        typedef LFS LocalFunctionSpace;
        typedef typename LFS::Traits::DOFIndex DOFIndex;
        typedef typename LFS::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex ContainerIndex;

        ConstLocalView()
          : _container(nullptr)
          , _lfs_cache(nullptr)
        {}

        ConstLocalView(const ISTLBlockVectorContainer& container)
          : _container(&container)
          , _lfs_cache(nullptr)
        {}

        void attach(const ISTLBlockVectorContainer& container)
        {
          _container = &container;
        }

        void detach()
        {
          _container = nullptr;
        }

        void bind(const LFSCache& lfs_cache)
        {
          _lfs_cache = &lfs_cache;
        }

        void unbind()
        {
        }

        size_type size() const
        {
          return _lfs_cache->size();
        }

        template<typename LC>
        void read(LC& local_container) const
        {
          for (size_type i = 0; i < size(); ++i)
            {
              accessBaseContainer(local_container)[i] = (*_container)[_lfs_cache->containerIndex(i)];
            }
        }

        template<typename ChildLFS, typename LC>
        void read(const ChildLFS& child_lfs, LC& local_container) const
        {
          for (size_type i = 0; i < child_lfs.size(); ++i)
            {
              const size_type local_index = child_lfs.localIndex(i);
              accessBaseContainer(local_container)[local_index] = (*_container)[_lfs_cache->containerIndex(local_index)];
            }
        }

        template<typename ChildLFS, typename LC>
        void read_sub_container(const ChildLFS& child_lfs, LC& local_container) const
        {
          for (size_type i = 0; i < child_lfs.size(); ++i)
            {
              const size_type local_index = child_lfs.localIndex(i);
              accessBaseContainer(local_container)[i] = (*_container)[_lfs_cache->containerIndex(local_index)];
            }
        }


        const ElementType& operator[](size_type i) const
        {
          return (*_container)[_lfs_cache->containerIndex(i)];
        }

        const ElementType& operator[](const DOFIndex& di) const
        {
          return (*_container)[_lfs_cache->containerIndex(di)];
        }

        const ElementType& operator[](const ContainerIndex& ci) const
        {
          return (*_container)[ci];
        }

        const ISTLBlockVectorContainer& global_container() const
        {
          return *_container;
        }


      private:

        const ISTLBlockVectorContainer* _container;
        const LFSCache* _lfs_cache;

      };


      ISTLBlockVectorContainer(const ISTLBlockVectorContainer& rhs)
        : _gfs(rhs._gfs)
        , _container(make_shared<Container>(_gfs.ordering().blockCount()))
      {
        istl::dispatch_vector_allocation(_gfs.ordering(),*_container,typename GFS::Ordering::ContainerAllocationTag());
        (*_container) = rhs.base();
      }

      ISTLBlockVectorContainer (const GFS& gfs, tags::attached_container = tags::attached_container())
        : _gfs(gfs)
        , _container(make_shared<Container>(gfs.ordering().blockCount()))
      {
        istl::dispatch_vector_allocation(gfs.ordering(),*_container,typename GFS::Ordering::ContainerAllocationTag());
      }

      //! Creates an ISTLBlockVectorContainer without allocating an underlying ISTL vector.
      ISTLBlockVectorContainer(const GFS& gfs, tags::unattached_container)
        : _gfs(gfs)
      {}

      ISTLBlockVectorContainer (const GFS& gfs, const E& e)
        : _gfs(gfs)
        , _container(make_shared<Container>(gfs.ordering().blockCount()))
      {
        istl::dispatch_vector_allocation(gfs.ordering(),*_container,typename GFS::Ordering::ContainerAllocationTag());
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

      ISTLBlockVectorContainer& operator= (const ISTLBlockVectorContainer& r)
      {
        if (this == &r)
          return *this;
        if (attached())
          {
            (*_container) = r.base();
          }
        else
          {
            _container = make_shared<Container>(r.base());
          }
        return *this;
      }

      ISTLBlockVectorContainer& operator= (const E& e)
      {
        (*_container)=e;
        return *this;
      }

      ISTLBlockVectorContainer& operator*= (const E& e)
      {
        (*_container)*=e;
        return *this;
      }


      ISTLBlockVectorContainer& operator+= (const E& e)
      {
        (*_container)+=e;
        return *this;
      }

      ISTLBlockVectorContainer& operator+= (const ISTLBlockVectorContainer& e)
      {
        (*_container)+= e.base();
        return *this;
      }

      ISTLBlockVectorContainer& operator-= (const ISTLBlockVectorContainer& e)
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

      E operator*(const ISTLBlockVectorContainer& y) const
      {
        return (*_container)*y.base();
      }

      E dot(const ISTLBlockVectorContainer& y) const
      {
        return _container->dot(y.base());
      }

      ISTLBlockVectorContainer& axpy(const E& a, const ISTLBlockVectorContainer& y)
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
      shared_ptr<Container> _container;
    };





#ifndef DOXYGEN

    // helper struct invoking the GFS tree -> ISTL vector reduction
    template<typename GFS, typename E>
    struct ISTLVectorSelectorHelper
    {

      typedef typename TypeTree::AccumulateType<
        GFS,
        istl::vector_creation_policy<E>
        >::type vector_descriptor;

      typedef ISTLBlockVectorContainer<GFS,typename vector_descriptor::vector_type> Type;

    };

    template<ISTLParameters::Blocking blocking, std::size_t block_size, typename GFS, typename E>
    struct BackendVectorSelectorHelper<ISTLVectorBackend<blocking,block_size>, GFS, E>
      : public ISTLVectorSelectorHelper<GFS,E>
    {};

#endif // DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTLVECTORBACKEND_HH
