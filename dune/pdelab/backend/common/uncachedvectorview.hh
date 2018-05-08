// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_COMMON_UNCACHEDVECTORVIEW_HH
#define DUNE_PDELAB_BACKEND_COMMON_UNCACHEDVECTORVIEW_HH

#include <dune/common/typetraits.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/pdelab/common/globalVariable.hh>

namespace Dune {
  namespace PDELab {


    template<typename V, typename LFSC>
    struct ConstUncachedVectorView
    {

      typedef typename std::remove_const<V>::type Container;
      typedef LFSC LFSCache;

      typedef typename Container::E ElementType;
      typedef typename Container::size_type size_type;
      typedef typename LFSCache::DOFIndex DOFIndex;
      typedef typename LFSCache::ContainerIndex ContainerIndex;


      ConstUncachedVectorView()
        : _container(nullptr)
        , _lfs_cache(nullptr)
      {}

      ConstUncachedVectorView(V& container)
        : _container(&container)
        , _lfs_cache(nullptr)
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
      }

      void unbind()
      {
      }

      size_type size() const
      {
        return cache().size();
      }

      template<typename LC>
      void read(LC& local_container) const
      {
        if (evilGlobalVariable) {

          auto refEl = Dune::ReferenceElements<double, 2>::general(Dune::GeometryTypes::cube(2));

          auto& coeffs = cache().localFunctionSpace().finiteElement().localCoefficients();

          for (int c = 0; c < refEl.dimension + 1; ++c) {
            for (int s = 0; s < refEl.size(c); ++s) {
              // evaluate consecutive index of subentity
              auto container_index = cache().containerIndex(s, c);
              auto local_index = coeffs.localDOF(Dune::LocalKey(s, c, 0));
              auto stride = coeffs.stride(s, c);
              auto chunk_size = coeffs.chunk_size(s, c);
              for (int i = 0, j = 0; i < coeffs.size_index(s, c); i+=chunk_size) {
                for (; j < chunk_size; ++j) {
                  // store data
                  accessBaseContainer(local_container)[local_index++] = container()[container_index];
                  container_index[0]++;
                }
                local_index += stride - 1;
              }
            }
          }
        } else {
          for (size_type i = 0; i < size(); ++i)
          {
            accessBaseContainer(local_container)[i] = container()[cache().containerIndex(i)];
          }
        }
      }

      template<typename ChildLFS, typename LC>
      void read(const ChildLFS& child_lfs, LC& local_container) const
      {
        for (size_type i = 0; i < child_lfs.size(); ++i)
          {
            const size_type local_index = child_lfs.localIndex(i);
            accessBaseContainer(local_container)[local_index] = container()[cache().containerIndex(local_index)];
          }
      }

      template<typename ChildLFS, typename LC>
      void read_sub_container(const ChildLFS& child_lfs, LC& local_container) const
      {
        for (size_type i = 0; i < child_lfs.size(); ++i)
          {
            const size_type local_index = child_lfs.localIndex(i);
            accessBaseContainer(local_container)[i] = container()[cache().containerIndex(local_index)];
          }
      }


      const ElementType& operator[](size_type i) const
      {
        return container()[cache().containerIndex(i)];
      }


      // disable this function if DOFIndex and ContainerIndex have the same type - required for interoperability
      // with function spaces based on dune-functions bases
      template<typename DI>
      std::enable_if_t<
        (std::is_same<DI,DOFIndex>{} and not std::is_same<DI,ContainerIndex>{}),
        const ElementType&
        >
      operator[](const DI& di) const
      {
        return container()[cache().containerIndex(di)];
      }


      const ElementType& operator[](const ContainerIndex& ci) const
      {
        return container()[ci];
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

    };


    template<typename V, typename LFSC>
    struct UncachedVectorView
      : public ConstUncachedVectorView<V,LFSC>
    {

      typedef V Container;
      typedef typename Container::ElementType ElementType;
      typedef typename Container::size_type size_type;

      typedef LFSC LFSCache;
      typedef typename LFSCache::DOFIndex DOFIndex;
      typedef typename LFSCache::ContainerIndex ContainerIndex;

      using ConstUncachedVectorView<V,LFSC>::cache;
      using ConstUncachedVectorView<V,LFSC>::size;

      // Explicitly pull in operator[] from the base class to work around a problem
      // with clang not finding the const overloads of the operator from the base class.
      using ConstUncachedVectorView<V,LFSC>::operator[];

      UncachedVectorView()
      {}

      UncachedVectorView(Container& container)
        : ConstUncachedVectorView<V,LFSC>(container)
      {}

      template<typename LC>
      void write(const LC& local_container)
      {
        for (size_type i = 0; i < size(); ++i)
          {
            container()[cache().containerIndex(i)] = accessBaseContainer(local_container)[i];
          }
      }

      template<typename LC>
      void add(const LC& local_container)
      {
        if (evilGlobalVariable) {

          auto refEl = Dune::ReferenceElements<double, 2>::general(Dune::GeometryTypes::cube(2));

          auto& coeffs = cache().localFunctionSpace().finiteElement().localCoefficients();

          for (int c = 0; c < refEl.dimension + 1; ++c) {
            for (int s = 0; s < refEl.size(c); ++s) {
              // evaluate consecutive index of subentity
              auto container_index = cache().containerIndex(s, c);
              auto local_index = coeffs.localDOF(Dune::LocalKey(s, c, 0));
              auto stride = coeffs.stride(s, c);
              auto chunk_size = coeffs.chunk_size(s, c);
              for (int i = 0, j = 0; i < coeffs.size_index(s, c); i+=chunk_size) {
                for (; j < chunk_size; ++j) {
                  // store data
                  container()[container_index] += accessBaseContainer(local_container)[local_index++];
                  container_index[0]++;
                }
                local_index += stride - 1;
              }
            }
          }
        } else
          for (size_type i = 0; i < size(); ++i)
          {
            container()[cache().containerIndex(i)] += accessBaseContainer(local_container)[i];
          }
      }



      template<typename ChildLFS, typename LC>
      void write(const ChildLFS& child_lfs, const LC& local_container)
      {
        for (size_type i = 0; i < child_lfs.size(); ++i)
          {
            const size_type local_index = child_lfs.localIndex(i);
            container()[cache().containerIndex(local_index)] = accessBaseContainer(local_container)[local_index];
          }
      }

      template<typename ChildLFS, typename LC>
      void add(const ChildLFS& child_lfs, const LC& local_container)
      {
        for (size_type i = 0; i < child_lfs.size(); ++i)
          {
            const size_type local_index = child_lfs.localIndex(i);
            container()[cache().containerIndex(local_index)] += accessBaseContainer(local_container)[local_index];
          }
      }




      template<typename ChildLFS, typename LC>
      void write_sub_container(const ChildLFS& child_lfs, const LC& local_container)
      {
        for (size_type i = 0; i < child_lfs.size(); ++i)
          {
            const size_type local_index = child_lfs.localIndex(i);
            container()[cache().containerIndex(local_index)] = accessBaseContainer(local_container)[i];
          }
      }

      template<typename ChildLFS, typename LC>
      void add_sub_container(const ChildLFS& child_lfs, const LC& local_container)
      {
        for (size_type i = 0; i < child_lfs.size(); ++i)
          {
            const size_type local_index = child_lfs.localIndex(i);
            container()[cache().containerIndex(local_index)] += accessBaseContainer(local_container)[i];
          }
      }

      void commit()
      {
      }


      ElementType& operator[](size_type i)
      {
        return container()[cache().containerIndex(i)];
      }

      // disable this function if DOFIndex and ContainerIndex have the same type - required for interoperability
      // with function spaces based on dune-functions bases
      template<typename DI>
      std::enable_if_t<
        (std::is_same<DI,DOFIndex>{} and not std::is_same<DI,ContainerIndex>{}),
        ElementType&
        >
      operator[](const DOFIndex& di)
      {
        return container()[cache().containerIndex(di)];
      }


      ElementType& operator[](const ContainerIndex& ci)
      {
        return container()[ci];
      }


      Container& container()
      {
        return *(this->_container);
      }


    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_COMMON_UNCACHEDVECTORVIEW_HH
