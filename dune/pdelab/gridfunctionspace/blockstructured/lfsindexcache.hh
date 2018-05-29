// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_LFSINDEXCACHE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_LFSINDEXCACHE_HH

#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/blockstructured/indexWrapper.hh>
#include <dune/common/power.hh>
#include <map>
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
          , localCoefficients(Dune::TypeTree::TreeInfo<LFS>::leafCount)
      {
      }

      void update()
      {
        initializeLocalCoefficients();

        auto refEl = Dune::ReferenceElements<double, d>::general(Dune::GeometryTypes::cube(d));

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

      constexpr std::size_t numberOfLeafs() const
      {
        return Dune::TypeTree::TreeInfo<LFS>::leafCount;
      }

      constexpr std::size_t codims() const
      {
        return d + 1;
      }

      std::size_t subentities(std::size_t codim) const
      {
        const auto refEl = Dune::ReferenceElements<double, d>::general(Dune::GeometryTypes::cube(d));
        return refEl.size(codim);
      }

      std::size_t sizeOfLocalDOFs(size_type leaf, size_type s, size_type c) const
      {
        return localCoefficients[leaf]->container[c][s].size();
      }

      std::size_t localIndex(size_type leaf, size_type s, size_type c, size_type i) const
      {
        return localCoefficients[leaf]->container[c][s][i] + localDOFsOffset[leaf];
      }

    private:

      void initializeLocalCoefficients()
      {
        TypeTree::forEachLeafNode(_lfs, [this] (auto& Node, auto& TreePath){
          const auto& fe = Node.finiteElement();
          const auto qkDescriptor = BlockstructuredQkDescriptor(fe);
          inverseLocalCoefficientsMap.try_emplace(qkDescriptor, fe);
        });

        TypeTree::forEachLeafNode(_lfs, [this] (auto& Node, auto& TreePath){
          const auto& fe = Node.finiteElement();
          const auto qkDescriptor = BlockstructuredQkDescriptor(fe);
          localCoefficients[Node.offsetLeafs] = &inverseLocalCoefficientsMap.at(qkDescriptor);
        });
      }

      constexpr static std::size_t d = GFS::Traits::GridView::dimension;

      const LFS& _lfs;
      std::vector<Dune::Blockstructured::SubentityWiseIndexWrapper<CI, d>> globalContainerIndices;
      std::map<BlockstructuredQkDescriptor, InverseQkLocalCoefficients<d>> inverseLocalCoefficientsMap;
      std::vector<std::size_t> localDOFsOffset;
      std::vector<const InverseQkLocalCoefficients<d>*> localCoefficients;
    };
  }
}


#endif //DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_LFSINDEXCACHE_HH
