//
// Created by marckoch on 09/05/18.
//

#ifndef DUNE_PDELAB_LFSINDEXCACHE_HH
#define DUNE_PDELAB_LFSINDEXCACHE_HH

#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/blockstructured/indexWrapper.hh>
#include <dune/common/power.hh>
#include "inversecoefficients.hh"

namespace Dune{
  namespace Blockstructured{

    struct BlockstructuredQkDescriptor{
      int k;
      int blocks;

      BlockstructuredQkDescriptor() = default;

      template <typename FE>
      BlockstructuredQkDescriptor(const FE&)
          : k(FE::k), blocks(FE::blocks)
      {}

      bool operator<(const BlockstructuredQkDescriptor& rhs) const
      {
        return (k < rhs.k) or ((k == rhs.k) and (blocks < rhs.blocks));
      }

    };


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
          , localDOFsOffset(Dune::TypeTree::TreeInfo<LFS>::leafCount)
          , qkDescriptors(Dune::TypeTree::TreeInfo<LFS>::leafCount)
      {
      }

      void update()
      {
        initializeInverseLocalCoefficientsMap();

        auto refEl = Dune::ReferenceElements<double,2>::general(Dune::GeometryTypes::cube(2));

        auto& subentityWiseDOFs = *_lfs._subentityWiseDOFs_ptr;

        globalContainerIndices.clear();
        globalContainerIndices.resize(numberOfLeafs());

        TypeTree::forEachLeafNode(_lfs, [this,&refEl,&subentityWiseDOFs] (auto& Node, auto& TreePath){
          const auto leaf = Node.offsetLeafs;

          localDOFsOffset[leaf] = Node.offset;

          for (int c = 0; c < refEl.dimension + 1; ++c)
            for (int s = 0; s < refEl.size(c); ++s)
              // evaluate consecutive index of subentity
              _lfs.gridFunctionSpace().ordering().mapIndex(subentityWiseDOFs[leaf].indexView(s, c),
                                                           globalContainerIndices[leaf].index(s, c));
        });
      }

      const CI& containerIndex(size_type leaf, size_type s, size_type c) const
      {
        return globalContainerIndices[leaf].index(s, c);
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
        return Dune::TypeTree::TreeInfo<LFS>::leafCount;
      }

      std::size_t sizeOfLocalDOFs(size_type leaf, size_type s, size_type c) const
      {
        return inverseLocalCoefficientsMap.at(qkDescriptors[leaf]).container[c][s].size();
      }

      std::size_t localIndex(size_type leaf, size_type s, size_type c, size_type i) const
      {
        return inverseLocalCoefficientsMap.at(qkDescriptors[leaf]).container[c][s][i] + localDOFsOffset[leaf];
      }

    private:

      void initializeInverseLocalCoefficientsMap()
      {
        TypeTree::forEachLeafNode(_lfs, [this] (auto& Node, auto& TreePath){
          const auto& fe = Node.finiteElement();
          this->qkDescriptors[Node.offsetLeafs] = BlockstructuredQkDescriptor(fe);

          inverseLocalCoefficientsMap.try_emplace(this->qkDescriptors[Node.offsetLeafs], fe);
        });
      }

      const LFS& _lfs;
      std::vector<Dune::Blockstructured::SubentityWiseIndexWrapper<CI, GFS::Traits::GridView::dimension>> globalContainerIndices;
      std::map<BlockstructuredQkDescriptor, InverseQkLocalCoefficients<GFS::Traits::GridView::dimension>> inverseLocalCoefficientsMap;
      std::vector<std::size_t> localDOFsOffset;
      std::vector<BlockstructuredQkDescriptor> qkDescriptors;
    };
  }
}


#endif //DUNE_PDELAB_LFSINDEXCACHE_HH
