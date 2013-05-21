// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTLVECTORBACKEND_HH
#define DUNE_ISTLVECTORBACKEND_HH

#include <vector>
#include <stack>

#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>

#include <dune/istl/bvector.hh>

#include <dune/pdelab/backend/istl/tags.hh>
#include <dune/pdelab/backend/istl/vectoriterator.hh>

#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>


#include "backendselector.hh"
#include "istlmatrixbackend.hh"


namespace Dune {
  namespace PDELab {

    // Recursive accessors for vector entries using tag dispatch

#ifndef DOXYGEN // All of the following functions are mere implementation details

    namespace istl {

      template<typename CI, typename Block>
      typename Block::field_type&
      access_vector_element(tags::field_vector_1, Block& b, const CI& ci, int i)
      {
        assert(i == -1);
        return b[0];
      }

      template<typename CI, typename Block>
      typename Block::field_type&
      access_vector_element(tags::field_vector_n, Block& b, const CI& ci, int i)
      {
        assert(i == 0);
        return b[ci[0]];
      }

      template<typename CI, typename Block>
      typename Block::field_type&
      access_vector_element(tags::block_vector, Block& b, const CI& ci, int i)
      {
        return access_vector_element(container_tag(b[ci[i]]),b[ci[i]],ci,i-1);
      }


      template<typename CI, typename Block>
      const typename Block::field_type&
      access_vector_element(tags::field_vector_1, const Block& b, const CI& ci, int i)
      {
        assert(i == -1);
        return b[0];
      }

      template<typename CI, typename Block>
      const typename Block::field_type&
      access_vector_element(tags::field_vector_n, const Block& b, const CI& ci, int i)
      {
        assert(i == 0);
        return b[ci[0]];
      }

      template<typename CI, typename Block>
      const typename Block::field_type&
      access_vector_element(tags::block_vector, const Block& b, const CI& ci, int i)
      {
        return access_istl_vector_element(container_tag(b[ci[i]]),b[ci[i]],ci,i-1);
      }


      template<typename Vector>
      void resize_vector(tags::block_vector, Vector& v, std::size_t size, bool copy_values)
      {
        v.resize(size,copy_values);
      }

      template<typename Vector>
      void resize_vector(tags::field_vector, Vector& v, std::size_t size, bool copy_values)
      {
      }

      template<typename DI, typename GDI, typename CI, typename Container>
      void allocate_vector(tags::field_vector, const OrderingBase<DI,GDI,CI>& ordering, Container& c)
      {
      }

      template<typename DI, typename GDI, typename CI, typename Container>
      void allocate_vector(tags::block_vector, const OrderingBase<DI,GDI,CI>& ordering, Container& c)
      {
        for (std::size_t i = 0; i < ordering.childOrderingCount(); ++i)
          {
            if (ordering.containerBlocked())
              {
                resize_vector(container_tag(c[i]),c[i],ordering.childOrdering(i).blockCount(),false);
                allocate_vector(container_tag(c[i]),ordering.childOrdering(i),c[i]);
              }
            else
              allocate_vector(container_tag(c),ordering.childOrdering(i),c);
          }
      }

      template<typename Ordering, typename Container>
      void dispatch_vector_allocation(const Ordering& ordering, Container& c, HierarchicContainerAllocationTag tag)
      {
        allocate_vector(container_tag(c),ordering,c);
      }

      template<typename Ordering, typename Container>
      void dispatch_vector_allocation(const Ordering& ordering, Container& c, FlatContainerAllocationTag tag)
      {
        resize_vector(container_tag(c),c,ordering.blockCount(),false);
      }



    } // namespace istl

#endif // DOXYGEN

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

      ISTLBlockVectorContainer (const GFS& gfs)
        : _gfs(gfs)
        , _container(make_shared<Container>(gfs.ordering().blockCount()))
      {
        istl::dispatch_vector_allocation(gfs.ordering(),*_container,typename GFS::Ordering::ContainerAllocationTag());
      }

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
        (*_container) = r.base();
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


    namespace istl {

      // ********************************************************************************
      // Helper functions for uniform access to ISTL containers
      //
      // The following suite of raw() functions should be used in places where an
      // algorithm might work on either the bare ISTL container or the PDELab
      // wrapper and has to access the bare container.
      // ********************************************************************************

      //! Returns the raw ISTL object associated with v, or v itself it is already an ISTL object.
      template<typename V>
      V& raw(V& v)
      {
        return v;
      }

      //! Returns the raw ISTL object associated with v, or v itself it is already an ISTL object.
      template<typename V>
      const V& raw(const V& v)
      {
        return v;
      }

      //! Returns the raw ISTL type associated with C, or C itself it is already an ISTL type.
      template<typename C>
      struct raw_type
      {
        typedef C type;
      };

#ifndef DOXYGEN

      template<typename GFS, typename C>
      typename ISTLBlockVectorContainer<GFS,C>::Container&
      raw(ISTLBlockVectorContainer<GFS,C>& v)
      {
        return v.base();
      }

      template<typename GFS, typename C>
      const typename ISTLBlockVectorContainer<GFS,C>::Container&
      raw(const ISTLBlockVectorContainer<GFS,C>& v)
      {
        return v.base();
      }

      template<typename GFSU, typename GFSV, typename C>
      typename ISTLMatrixContainer<GFSU,GFSV,C>::Container&
      raw(ISTLMatrixContainer<GFSU,GFSV,C>& m)
      {
        return m.base();
      }

      template<typename GFSU, typename GFSV, typename C>
      const typename ISTLMatrixContainer<GFSU,GFSV,C>::Container&
      raw(const ISTLMatrixContainer<GFSU,GFSV,C>& m)
      {
        return m.base();
      }

      template<typename GFS, typename C>
      struct raw_type<ISTLBlockVectorContainer<GFS,C> >
      {
        typedef C type;
      };

      template<typename GFSU, typename GFSV, typename C>
      struct raw_type<ISTLMatrixContainer<GFSU,GFSV,C> >
      {
        typedef C type;
      };

#endif // DOXYGEN


      // ********************************************************************************
      // TMPs for deducing ISTL block structure from GFS backends
      // ********************************************************************************


#ifndef DOXYGEN // All of the following TMP magic is about irrelevant to the user as can be...

      // tag dispatch switch on GFS tag for per-node functor - general version
      template<typename E,typename Node, typename Tag>
      struct vector_descriptor_helper
      {
        // export backend type, as the actual TMP is in the parent reduction functor
        typedef typename Node::Traits::Backend type;
      };

      // descriptor for backends of leaf spaces collecting various information about
      // possible blocking structures
      template<typename E, typename Backend>
      struct leaf_vector_descriptor
      {

        dune_static_assert(Backend::Traits::block_type != ISTLParameters::dynamic_blocking,
                           "Dynamcially blocked leaf spaces are not supported by this backend.");

        // flag for sibling reduction - always true in the leaf case
        static const bool support_no_blocking = true;

        // flag indicating whether the associated vector type supports cascading
        // the static blocking further up the tree (i.e. create larger static blocks
        // at the parent node level. Due to ISTL limitations, this only works once in
        // the hierarchy, so we only support cascading if we don't already do static
        // blocking at the current level.
        static const bool support_cascaded_blocking =
          Backend::Traits::block_type == ISTLParameters::no_blocking; // FIXME

        // The static block size of the associated vector
        static const std::size_t block_size =
          Backend::Traits::block_type == ISTLParameters::static_blocking
          ? Backend::Traits::block_size
          : 1;

        // The cumulative block size is used by the algorithm to calculate total block
        // size over several children for cascaded blocking. Right now, this will always be set to
        // the block size passed in by the user, but it might also be possible to provide this
        // information in the FiniteElementMap and provide automatic blocking detection.
        static const std::size_t cumulative_block_size = Backend::Traits::block_size;

        // The element type for the vector.
        typedef E element_type;

        // The ISTL vector type associated with the current subtree.
        typedef BlockVector<FieldVector<E,block_size> > vector_type;

      };

      // Tag dispatch for leaf spaces - extract leaf descriptor.
      template<typename E, typename Node>
      struct vector_descriptor_helper<E,Node,LeafGridFunctionSpaceTag>
      {
        typedef leaf_vector_descriptor<E,typename Node::Traits::Backend> type;
      };

      // the actual functor
      template<typename E>
      struct extract_vector_descriptor
      {

        template<typename Node, typename TreePath>
        struct doVisit
        {
          // visit all nodes
          static const bool value = true;
        };

        template<typename Node, typename TreePath>
        struct visit
        {
          // forward to actual implementation via tag dispatch
          typedef typename vector_descriptor_helper<E,Node,typename Node::ImplementationTag>::type type;
        };

      };

      // Descriptor for combining sibling nodes in the tree
      template<typename Sibling, typename Child>
      struct cascading_vector_descriptor
      {

        // We only support cascaded blocking if all children support it
        static const bool support_cascaded_blocking =
          Sibling::support_cascaded_blocking &&
          Child::support_cascaded_blocking;

        // ISTL requires a single, globally constant blocking structure
        // for its containers, so we make sure the siblings don't disagree
        // on it.
        static const bool support_no_blocking =
          (Sibling::support_no_blocking &&
          is_same<
            typename Sibling::vector_type,
            typename Child::vector_type
           >::value);

        // block size
        static const std::size_t block_size =
          support_no_blocking ? Sibling::block_size : 1;

        // The element type for the vector.
        typedef typename Sibling::element_type element_type;

        // Accumulate total block size of all siblings
        static const std::size_t cumulative_block_size =
          Sibling::cumulative_block_size + Child::cumulative_block_size;

        // The ISTL vector type associated with the current subtree.
        typedef BlockVector<FieldVector<element_type,block_size> > vector_type;

      };


      // Switch that turns off standard reduction for the first child of a node.
      // Default case: do the standard reduction.
      template<typename D1, typename D2>
      struct initial_reduction_switch
      {
        typedef cascading_vector_descriptor<D1,D2> type;
      };

      // specialization for first child
      template<typename D2>
      struct initial_reduction_switch<void,D2>
      {
        typedef D2 type;
      };

      // sibling reduction functor
      struct combine_vector_descriptor_siblings
      {

        template<typename D1, typename D2>
        struct reduce
          : public initial_reduction_switch<D1,D2>
        {};

      };

      // Data part of child -> parent reduction descriptor
      template<typename Child, typename Backend>
      struct parent_child_vector_descriptor_data
      {

        // If all our have a common blocking structure, we can just
        // concatenate them without doing any blocking
        static const bool support_no_blocking =
          Child::support_no_blocking;

        // We support cascaded blocking if neither we nor any of our
        // children are blocked yet.
        static const bool support_cascaded_blocking =
          Child::support_cascaded_blocking &&
          Backend::Traits::block_type == ISTLParameters::no_blocking;

        // Throw an assertion if the user requests static blocking at this level,
        // but we cannot support it.
        dune_static_assert((Backend::Traits::block_type != ISTLParameters::static_blocking) ||
                           Child::support_cascaded_blocking,
                           "invalid blocking structure.");

        // If we block statically, we create bigger blocks, otherwise the
        // block size doesn't change.
        static const std::size_t block_size =
          Backend::Traits::block_type == ISTLParameters::static_blocking
          ? Child::cumulative_block_size
          : Child::block_size;

        // Just forward this...
        static const std::size_t cumulative_block_size =
          Child::cumulative_block_size;

        // The element type for the vector.
        typedef typename Child::element_type element_type;

        // The ISTL vector type associated with our subtrees.
        typedef typename Child::vector_type child_vector_type;

      };

      // dispatch switch on blocking type - prototype
      template<typename Data, ISTLParameters::Blocking>
      struct parent_child_vector_descriptor;

      // dispatch switch on blocking type - no blocking case
      template<typename Data>
      struct parent_child_vector_descriptor<
        Data,
        ISTLParameters::no_blocking
        >
        : public Data
      {
        dune_static_assert(Data::support_no_blocking,
                           "Cannot combine incompatible child block structures without static blocking. "
                           "Did you want to apply static blocking at this level?");

        // Just forward the child vector type
        typedef typename Data::child_vector_type vector_type;
      };

      // dispatch switch on blocking type - dynamic blocking case
      template<typename Data>
      struct parent_child_vector_descriptor<
        Data,
        ISTLParameters::dynamic_blocking
        >
        : public Data
      {
        dune_static_assert(Data::support_no_blocking,
                           "Incompatible child block structures detected, cannot perform dynamic blocking. "
                           "Did you want to apply static blocking at this level?");

        // Wrap the child vector type in another BlockVector
        typedef BlockVector<typename Data::child_vector_type> vector_type;
      };

      // dispatch switch on blocking type - static blocking case
      template<typename Data>
      struct parent_child_vector_descriptor<
        Data,
        ISTLParameters::static_blocking
        >
        : public Data
      {
        // build new block vector with large field block size
        typedef BlockVector<
          FieldVector<
            typename Data::element_type,
            Data::block_size
            >
          > vector_type;
      };

      // Child - parent reduction functor
      struct combine_vector_descriptor_parent
      {

        template<typename Child, typename Backend>
        struct reduce
        {

          struct type
            : public parent_child_vector_descriptor<parent_child_vector_descriptor_data<
                                                      Child,
                                                      Backend>,
                                                    Backend::Traits::block_type
                                                    >
          {};
        };

      };

      // policy describing the GFS tree -> ISTL vector reduction
      template<typename E>
      struct vector_creation_policy
        : public TypeTree::TypeAccumulationPolicy<extract_vector_descriptor<E>,
                                                  combine_vector_descriptor_siblings,
                                                  void,
                                                  combine_vector_descriptor_parent,
                                                  TypeTree::bottom_up_reduction>
      {};

#endif // DOXYGEN

    } // namespace istl


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

#endif
