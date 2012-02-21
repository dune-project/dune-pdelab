// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTLVECTORBACKEND_HH
#define DUNE_ISTLVECTORBACKEND_HH

#include <vector>
#include <stack>
#include <unordered_map>

#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>

#include <dune/istl/bvector.hh>

#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/gridfunctionspace/lfscontainerindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>


#include "backendselector.hh"
#include "istlmatrixbackend.hh"


namespace Dune {
  namespace PDELab {

    template<std::size_t block_size = 1>
    struct ISTLFieldVectorBackend
    {

      dune_static_assert((block_size > 0),"block size for FieldVector has to be positive");


      typedef std::size_t size_type;

      static const size_type blockSize = block_size;

      struct Traits
      {
        static const size_type max_blocking_depth = block_size > 1 ? 1 : 0;
      };

      bool blocked() const
      {
        return block_size > 1;
      }
    };

    template<bool block = false>
    struct ISTLVectorBackend
    {

      typedef std::size_t size_type;

      static const bool is_blocked = block;

      struct Traits
      {
        static const size_type max_blocking_depth = block ? 1 : 0;
      };


      bool blocked() const
      {
        return block;
      }

    };

    template<typename E,typename Node, typename Tag>
    struct extract_istl_vector_helper
    {
      // export node type, as the actual calculation takes place in the reduction functor
      typedef Node type;
    };

    struct LeafGridFunctionSpaceTag;

    template<typename E, typename Node>
    struct extract_istl_vector_helper<E,Node,LeafGridFunctionSpaceTag>
    {
      typedef BlockVector<FieldVector<E,Node::Traits::Backend::blockSize> > type;
    };

    template<typename E>
    struct extract_istl_vector
    {

      template<typename Node, typename TreePath>
      struct doVisit
      {
        static const bool value = true;
      };

      template<typename Node, typename TreePath>
      struct visit
      {
        typedef typename extract_istl_vector_helper<E,Node,typename Node::ImplementationTag>::type type;
      };

    };

    struct combine_istl_vector_siblings
    {

      template<typename V1, typename V2>
      struct reduce
      {
        dune_static_assert((is_same<V1,V2>::value) || // identical child types
                           (is_same<V1,void>::value), // special case to catch start type
                           "This backend only supports sub-vectors of identical type.");

        typedef V2 type;
      };

    };

    struct combine_istl_vector_parent
    {

      template<typename ChildContainer, typename Node>
      struct reduce
      {
        typedef typename SelectType<
          Node::Traits::Backend::is_blocked,
          BlockVector<ChildContainer>,
          ChildContainer
          >::Type type;
      };

    };


    template<typename CI, typename Block>
    typename enable_if<Block::blocklevel != 1,typename Block::field_type&>::type
    access_istl_vector_element(Block& b, const CI& ci, std::size_t i)
    {
      return access_istl_vector_element(b[ci[i]],ci,i-1);
    }

    template<typename CI, typename Block>
    typename enable_if<Block::blocklevel == 1 && Block::dimension == 1,typename Block::field_type&>::type
    access_istl_vector_element(Block& b, const CI& ci, std::size_t i)
    {
      assert(i == 0);
      return b[0];
    }

    template<typename CI, typename Block>
    typename enable_if<Block::blocklevel == 1 && Block::dimension != 1,typename Block::field_type&>::type
    access_istl_vector_element(Block& b, const CI& ci, std::size_t i)
    {
      assert(i == 0);
      return b[ci[i]];
    }


    template<typename CI, typename Block>
    typename enable_if<Block::blocklevel != 1,const typename Block::field_type&>::type
    access_istl_vector_element(const Block& b, const CI& ci, std::size_t i)
    {
      return access_istl_vector_element(b[ci[i]],ci,i-1);
    }

    template<typename CI, typename Block>
    typename enable_if<Block::blocklevel == 1 && Block::dimension == 1,const typename Block::field_type&>::type
    access_istl_vector_element(const Block& b, const CI& ci, std::size_t i)
    {
      assert(i == 0);
      return b[0];
    }

    template<typename CI, typename Block>
    typename enable_if<Block::blocklevel == 1 && Block::dimension != 1,const typename Block::field_type&>::type
    access_istl_vector_element(const Block& b, const CI& ci, std::size_t i)
    {
      assert(i == 0);
      return b[ci[i]];
    }


    template<typename Vector>
    void resize_istl_vector(Vector& v, std::size_t size, bool copy_values)
    {
      v.resize(size,copy_values);
    }

    template<typename block_type, int block_size>
    void resize_istl_vector(FieldVector<block_type,block_size>& v, std::size_t size, bool copy_values)
    {
    }

    template<typename DI, typename CI, typename Container>
    typename enable_if<!is_same<typename Container::block_type,double>::value>::type
    allocate_istl_vector(const OrderingBase<DI,CI>& ordering, Container& c)
    {
      for (std::size_t i = 0; i < ordering.dynamic_child_count(); ++i)
        {
          if (ordering.container_blocked())
            {
              resize_istl_vector(c[i],ordering.dynamic_child(i).blockCount(),false);
              allocate_istl_vector(ordering.dynamic_child(i),c[i]);
            }
          else
            allocate_istl_vector(ordering.dynamic_child(i),c);
        }
    }

    template<typename DI, typename CI, typename Container>
    typename enable_if<is_same<typename Container::block_type,double>::value>::type
    allocate_istl_vector(const OrderingBase<DI,CI>& ordering, Container& c)
    {
    }

    template<typename GFS, typename C>
    class ISTLBlockVectorContainer
    {

    public:
      typedef typename C::field_type ElementType;
      typedef ElementType E;
      typedef C ContainerType;
      typedef GFS GridFunctionSpace;
      typedef ContainerType BaseT;
      typedef typename ContainerType::field_type field_type;
      typedef typename ContainerType::iterator iterator;
      typedef typename ContainerType::const_iterator const_iterator;
      typedef typename ContainerType::block_type block_type;
      typedef typename ContainerType::size_type size_type;

      typedef typename GFS::Ordering::Traits::ContainerIndex ContainerIndex;

      template<typename LFSCache>
      struct LocalView
      {

        dune_static_assert((is_same<typename LFSCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFS>::value),
                           "The LocalFunctionSpace passed to LocalView must belong to the underlying GridFunctionSpace.");

        typedef E ElementType;
        typedef typename LFSCache::LocalFunctionSpace LFS;
        typedef typename LFS::Traits::DOFIndex DOFIndex;
        typedef typename LFS::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex ContainerIndex;

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
              accessBaseContainer(local_container)[i] = (*_container)[_lfs_cache->container_index(i)];
            }
        }

        template<typename LC>
        void write(const LC& local_container)
        {
          for (size_type i = 0; i < size(); ++i)
            {
              (*_container)[_lfs_cache->container_index(i)] = accessBaseContainer(local_container)[i];
            }
        }

        template<typename LC>
        void add(const LC& local_container)
        {
          for (size_type i = 0; i < size(); ++i)
            {
              (*_container)[_lfs_cache->container_index(i)] += accessBaseContainer(local_container)[i];
            }
        }

        void commit()
        {
        }

        ElementType& operator[](size_type i)
        {
          return (*_container)[_lfs_cache->container_index(i)];
        }

        const ElementType& operator[](size_type i) const
        {
          return (*_container)[_lfs_cache->container_index(i)];
        }

        ElementType& operator[](const DOFIndex& di)
        {
          return (*_container)[_lfs_cache->container_index(di)];
        }

        const ElementType& operator[](const DOFIndex& di) const
        {
          return (*_container)[_lfs_cache->container_index(di)];
        }

        ElementType& operator[](const ContainerIndex& ci)
        {
          return (*_container)[ci];
        }

        const ElementType& operator[](const ContainerIndex& ci) const
        {
          return (*_container)[ci];
        }


      private:

        ISTLBlockVectorContainer* _container;
        const LFSCache* _lfs_cache;

      };


      template<typename LFSCache>
      struct ConstLocalView
      {

        dune_static_assert((is_same<typename LFSCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFS>::value),
                           "The LocalFunctionSpace passed to LocalView must belong to the underlying GridFunctionSpace.");

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
              accessBaseContainer(local_container)[i] = (*_container)[_lfs_cache->container_index(i)];
            }
        }

        const ElementType& operator[](size_type i) const
        {
          return _container[_lfs_cache->container_index(i)];
        }

        const ElementType& operator[](const DOFIndex& di) const
        {
          return _container[_lfs_cache->container_index(di)];
        }

        const ElementType& operator[](const ContainerIndex& ci) const
        {
          return _container[ci];
        }


      private:

        const ISTLBlockVectorContainer* _container;
        const LFSCache* _lfs_cache;

      };


      ISTLBlockVectorContainer (const GFS& gfs_)
        : container(gfs_.ordering()->blockCount())
      {
        allocate_istl_vector(*gfs_.ordering(),container);
      }

      ISTLBlockVectorContainer (const GFS& gfs_, const E& e)
        : container(gfs_.ordering()->blockCount())
      {
        allocate_istl_vector(*gfs_.ordering(),container);
        container=e;
      }

      size_type N() const
      {
        return container.N();
      }


      ISTLBlockVectorContainer& operator= (const E& e)
      {
        container=e;
        return *this;
      }

      ISTLBlockVectorContainer& operator*= (const E& e)
      {
        container*=e;
        return *this;
      }


      ISTLBlockVectorContainer& operator+= (const E& e)
      {
        container+=e;
        return *this;
      }

      ISTLBlockVectorContainer& operator+= (const ISTLBlockVectorContainer& e)
      {
        container+=e;
        return *this;
      }

      ISTLBlockVectorContainer& operator-= (const ISTLBlockVectorContainer& e)
      {
        container-=e;
        return *this;
      }

      block_type& operator[](std::size_t i)
      {
        return container[i];
      }

      const block_type& operator[](std::size_t i) const
      {
        return container[i];
      }

      E& operator[](const ContainerIndex& ci)
      {
        return access_istl_vector_element(container,ci,ci.size()-1);
      }

      const E& operator[](const ContainerIndex& ci) const
      {
        return access_istl_vector_element(container,ci,ci.size()-1);
      }

      typename Dune::template FieldTraits<E>::real_type two_norm() const
      {
        return container.two_norm();
      }

      typename Dune::template FieldTraits<E>::real_type one_norm() const
      {
        return container.one_norm();
      }

      typename Dune::template FieldTraits<E>::real_type infinity_norm() const
      {
        return container.infinity_norm();
      }

      E operator*(const ISTLBlockVectorContainer& y) const
      {
        return container*y.base();
      }

      ISTLBlockVectorContainer& axpy(const E& a, const ISTLBlockVectorContainer& y)
      {
        container.axpy(a, y.base());
        return *this;
      }

      // for debugging and AMG access
      ContainerType& base ()
      {
        return container;
      }

      const ContainerType& base () const
      {
        return container;
      }

      operator ContainerType&()
      {
        return container;
      }

      operator const ContainerType&() const
      {
        return container;
      }

      iterator begin()
      {
        return container.begin();
      }


      const_iterator begin() const
      {
        return container.begin();
      }

      iterator end()
      {
        return container.end();
      }


      const_iterator end() const
      {
        return container.end();
      }

      size_t flatsize() const
      {
        return container.dim();
      }

      template<typename X>
      void std_copy_to (std::vector<X>& x) const
      {
        // FIXME: do this hierachically
        size_t n = flatsize();
        x.resize(n);
        for (size_t i=0; i<n; i++)
          x[i] = container[i][i];
      }

      template<typename X>
      void std_copy_from (const std::vector<X>& x)
      {
        // FIXME: do this hierachically
        //test if x has the same size as the container
        assert (x.size() == flatsize());
        for (size_t i=0; i<flatsize(); i++)
          container[i][i] = x[i];
      }

    private:
      ContainerType container;
    };


    template<typename GFS, typename E>
    struct ISTLVectorSelectorHelper
    {

      typedef typename TypeTree::AccumulateType<
        GFS,
        extract_istl_vector<E>,
        combine_istl_vector_siblings,
        void,
        combine_istl_vector_parent
        >::type istl_container_type;

      typedef ISTLBlockVectorContainer<GFS,istl_container_type> Type;

    };

    template<std::size_t blockSize, typename GFS, typename E>
    struct BackendVectorSelectorHelper<ISTLFieldVectorBackend<blockSize>, GFS, E>
      : public ISTLVectorSelectorHelper<GFS,E>
    {};

    template<bool blocked, typename GFS, typename E>
    struct BackendVectorSelectorHelper<ISTLVectorBackend<blocked>, GFS, E>
      : public ISTLVectorSelectorHelper<GFS,E>
    {};


  } // namespace PDELab
} // namespace Dune

#endif
