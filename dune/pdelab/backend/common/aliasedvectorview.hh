// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_COMMON_ALIASEDVECTORVIEW_HH
#define DUNE_PDELAB_BACKEND_COMMON_ALIASEDVECTORVIEW_HH

#include <dune/common/typetraits.hh>
#include <dune/common/deprecated.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>

namespace Dune {
  namespace PDELab {


    template<typename V, typename LFSC>
    struct ConstAliasedVectorView
    {

      typedef typename std::remove_const<V>::type Container;
      typedef LFSC LFSCache;

      typedef typename Container::E ElementType;
      typedef typename Container::size_type size_type;
      typedef typename LFSCache::DOFIndex DOFIndex;
      typedef typename LFSCache::ContainerIndex ContainerIndex;

      using value_type = ElementType;


      ConstAliasedVectorView()
        : _container(nullptr)
        , _lfs_cache(nullptr)
        , _data(nullptr)
      {}

      ConstAliasedVectorView(V& container)
        : _container(&container)
        , _lfs_cache(nullptr)
        , _data(nullptr)
      {}

      void attach(V& container)
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
        _data = _container->data(lfs_cache);
      }

      const ElementType* data() const
      {
        return _data;
      }

      void unbind()
      {
        _lfs_cache = nullptr;
        _data = nullptr;
      }

      size_type size() const
      {
        return cache().size();
      }

      const ElementType& operator[](size_type i) const
      {
        return _data[i];
      }

      const ElementType& operator[](const ContainerIndex& ci) const
      {
        return container()[ci];
      }

      template<typename LFS>
      const ElementType& operator()(const LFS& lfs, size_type i) const
      {
        return this->_data[lfs.localIndex(i)];
      }

      const Container& container() const
      {
        return *_container;
      }

      const LFSCache& cache() const
      {
        return *_lfs_cache;
      }

    protected:

      V* _container;
      const LFSCache* _lfs_cache;
      typename std::conditional<
        std::is_const<V>::value,
        const ElementType*,
        ElementType*
        >::type _data;

    };


    template<typename V, typename LFSC>
    struct AliasedVectorView
      : public ConstAliasedVectorView<V,LFSC>
    {

      typedef V Container;
      typedef typename Container::ElementType ElementType;
      typedef typename Container::size_type size_type;

      typedef LFSC LFSCache;
      typedef typename LFSCache::DOFIndex DOFIndex;
      typedef typename LFSCache::ContainerIndex ContainerIndex;

      using value_type = ElementType;
      using weight_type = ElementType;

      using ConstAliasedVectorView<V,LFSC>::cache;
      using ConstAliasedVectorView<V,LFSC>::size;

      // Explicitly pull in operator[] from the base class to work around a problem
      // with clang not finding the const overloads of the operator from the base class.
      using ConstAliasedVectorView<V,LFSC>::operator[];

      // pull in const version of data access
      using ConstAliasedVectorView<V,LFSC>::data;

      AliasedVectorView()
        : weight_(1.0)
      {}

      AliasedVectorView(Container& container)
        : ConstAliasedVectorView<V,LFSC>(container)
        , weight_(1.0)
      {}

      void commit()
      {}

      template<typename LFS>
      void accumulate(const LFS& lfs, size_type n, value_type v)
      {
        this->_data[lfs.localIndex(n)] += v;
      }

      template<typename LFS>
      void rawAccumulate(const LFS& lfs, size_type n, value_type v)
      {
        accumulate(lfs,n,v);
      }

      ElementType& operator[](size_type i)
      {
        return this->_data[i];
      }

      ElementType& operator[](const ContainerIndex& ci)
      {
        return container()[ci];
      }

      ElementType* data()
      {
        return this->_data;
      }

      const ElementType* data() const
      {
        return this->_data;
      }

      Container& container()
      {
        return *(this->_container);
      }

      void setWeight(weight_type weight)
      {
        weight_ = weight;
      }

      weight_type weight()
      {
        return weight_;
      }

    private :
      weight_type weight_;
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_COMMON_ALIASEDVECTORVIEW_HH
