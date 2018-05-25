//
// Created by marckoch on 09/05/18.
//

#ifndef DUNE_PDELAB_LFSINDEXCACHE_HH
#define DUNE_PDELAB_LFSINDEXCACHE_HH

#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/blockstructured/indexWrapper.hh>

namespace Dune{
  namespace Blockstructured{

    template<int d>
    using SubentityWiseLocalIndexContainer = std::array<std::vector<std::vector<std::size_t>>, d + 1>;

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

        auto& subentityWiseDOFs = *_lfs._dof_index_storage_subentity_wise_ptr;

        const std::size_t nLeafs = subentityWiseDOFs.size();

        _container_index_storage_subentity_wise.clear();
        _container_index_storage_subentity_wise.resize(nLeafs);

        _local_index_storage_subentity_wise.clear();
        _local_index_storage_subentity_wise.resize(nLeafs);

        offset.resize(nLeafs);

        TypeTree::forEachLeafNode(_lfs, [this,&refEl,&subentityWiseDOFs] (auto& Node, auto& TreePath){
          const auto leaf = Node.offsetLeafs;

          this->offset[leaf] = Node.offset;

          this->_local_index_storage_subentity_wise[leaf] = &Node.finiteElement().localCoefficients().getLocalIndexContainer();
          for (int c = 0; c < refEl.dimension + 1; ++c)
            for (int s = 0; s < refEl.size(c); ++s)
              // evaluate consecutive index of subentity
              this->_lfs.gridFunctionSpace().ordering().mapIndex(subentityWiseDOFs[leaf].indexView(s, c),
                                                                 this->_container_index_storage_subentity_wise[leaf].index(s, c));
        });
      }

      const CI& containerIndex(size_type leaf, size_type s, size_type c) const
      {
        return _container_index_storage_subentity_wise[leaf].index(s, c);
      }


      const CI& containerIndex(size_type i) const
      {
        DUNE_THROW(Dune::NotImplemented, "Use this index cache by iterating over the reference element subentities");
        return {};
      }


      const CI& containerIndex(const DI& i) const
      {
        DUNE_THROW(Dune::NotImplemented, "Use this index cache by iterating over the reference element subentities");
        return {};
      }

      std::size_t numberOfLeafs() const
      {
        return _container_index_storage_subentity_wise.size();
      }

      std::size_t sizeOfLocalDOFs(size_type leaf, size_type s, size_type c) const
      {
        return (*_local_index_storage_subentity_wise[leaf])[c][s].size();
      }

      std::size_t localIndex(size_type leaf, size_type s, size_type c, size_type i) const
      {
        return (*_local_index_storage_subentity_wise[leaf])[c][s][i] + offset[leaf];
      }

    private:

      const LFS& _lfs;
      std::vector<Dune::Blockstructured::SubentityWiseIndexWrapper<CI>> _container_index_storage_subentity_wise;
      std::vector<std::size_t> offset;
      std::vector<const Dune::Blockstructured::SubentityWiseLocalIndexContainer<GFS::Traits::GridView::dimension>*>
          _local_index_storage_subentity_wise;

    };
  }
}


#endif //DUNE_PDELAB_LFSINDEXCACHE_HH
