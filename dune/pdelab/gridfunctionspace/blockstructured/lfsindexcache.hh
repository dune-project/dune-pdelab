//
// Created by marckoch on 09/05/18.
//

#ifndef DUNE_PDELAB_LFSINDEXCACHE_HH
#define DUNE_PDELAB_LFSINDEXCACHE_HH

#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

namespace Dune{
  namespace Blockstructured{

    template<typename LFS, typename C>
    class LFSIndexCache
        : public Dune::PDELab::LFSIndexCacheBase<LFS,C,typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag,false>
    {
    public:

      typedef typename LFS::Traits::GridFunctionSpace GFS;
      typedef typename GFS::Ordering Ordering;
      typedef typename Ordering::Traits::ContainerIndex ContainerIndex;
      typedef ContainerIndex CI;
      typedef typename Ordering::Traits::DOFIndex DOFIndex;
      typedef DOFIndex DI;
      typedef std::size_t size_type;

      typedef std::vector<CI> CIVector;

      using Base = Dune::PDELab::LFSIndexCacheBase<LFS,C,typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag,false>;


      LFSIndexCache(const LFS& lfs, const C& constraints, bool enable_constraints_caching)
          : Base(lfs, constraints, enable_constraints_caching)
          ,_lfs(lfs)
      {
      }

      void update()
      {
        auto refEl = Dune::ReferenceElements<double,2>::general(Dune::GeometryTypes::cube(2));

        for (int c = 0; c < refEl.dimension + 1; ++c) {
          _container_indices_sc[c].resize(refEl.size(c));
          for (int s = 0; s < refEl.size(c); ++s) {
            // evaluate consecutive index of subentity
            _container_indices_sc[c][s].clear();
            _lfs.gridFunctionSpace().ordering().mapIndex(_lfs._dof_indices_sc[c][s].view(), _container_indices_sc[c][s]);
          }
        }
      }

      const CI& containerIndex(size_type s, size_type c) const
      {
        return _container_indices_sc[c][s];
      }


      const CI& containerIndex(size_type i) const
      {
        DUNE_THROW(Dune::NotImplemented, "");
        return CI();
      }


      const CI& containerIndex(const DI& i) const
      {
        DUNE_THROW(Dune::NotImplemented, "");
        return CI();
      }

    private:

      const LFS& _lfs;
      std::array<CIVector, 3> _container_indices_sc;

    };
  }
}


#endif //DUNE_PDELAB_LFSINDEXCACHE_HH
