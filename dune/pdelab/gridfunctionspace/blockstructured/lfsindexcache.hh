//
// Created by marckoch on 09/05/18.
//

#ifndef DUNE_PDELAB_LFSINDEXCACHE_HH
#define DUNE_PDELAB_LFSINDEXCACHE_HH

#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/blockstructured/indexWrapper.hh>

namespace Dune{
  namespace Blockstructured{

    class BlockstructuredLFSCBase{};


    template<typename LFS, typename C>
    class LFSIndexCache
        : public Dune::PDELab::LFSIndexCacheBase<LFS,C,typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag,false>,
            public BlockstructuredLFSCBase
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

        _container_index_storage_subentity_wise.clear();
        _container_index_storage_subentity_wise.resize((*_lfs._dof_index_storage_subentity_wise_ptr).size());

        for (int leaf = 0; leaf < _container_index_storage_subentity_wise.size(); ++leaf)
          for (int c = 0; c < refEl.dimension + 1; ++c)
            for (int s = 0; s < refEl.size(c); ++s)
              // evaluate consecutive index of subentity
              _lfs.gridFunctionSpace().ordering().mapIndex((*_lfs._dof_index_storage_subentity_wise_ptr)[leaf].indexView(s, c),
                                                           _container_index_storage_subentity_wise[leaf].index(s, c));
      }

      const CI& containerIndex(size_type leaf, size_type s, size_type c) const
      {
        return _container_index_storage_subentity_wise[leaf].index(s, c);
      }


      const CI& containerIndex(size_type i) const
      {
        DUNE_THROW(Dune::NotImplemented, "Use this index cache in a structured way, by iterating over the reference element subentities");
        return CI();
      }


      const CI& containerIndex(const DI& i) const
      {
        DUNE_THROW(Dune::NotImplemented, "Use this index cache in a structured way, by iterating over the reference element subentities");
        return CI();
      }

      size_type numberOfLeafs() const
      {
        return _container_index_storage_subentity_wise.size();
      }

    private:

      const LFS& _lfs;
      std::vector<SubentityWiseIndexWrapper<CI>> _container_index_storage_subentity_wise;

    };
  }
}


#endif //DUNE_PDELAB_LFSINDEXCACHE_HH
