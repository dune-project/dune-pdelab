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
  namespace PDELab {
    namespace Blockstructured{

      template<typename LFS, typename C>
      class LFSIndexCache
          : public Dune::PDELab::LFSIndexCacheBase<LFS, C, typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag, false> {
      public:

        typedef typename LFS::Traits::GridFunctionSpace GFS;
        typedef typename GFS::Ordering Ordering;
        typedef typename Ordering::Traits::ContainerIndex ContainerIndex;
        typedef ContainerIndex CI;
        typedef typename Ordering::Traits::DOFIndex DOFIndex;
        typedef DOFIndex DI;
        typedef std::size_t size_type;

        typedef std::vector<CI> CIVector;

        using Base = Dune::PDELab::LFSIndexCacheBase<LFS, C, typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag, false>;


        LFSIndexCache(const LFS &lfs, const C &constraints, bool enable_constraints_caching)
            : Base(lfs, constraints, enable_constraints_caching),
              localDOFsOffset(Dune::TypeTree::TreeInfo<LFS>::leafCount),
              localCoefficients(Dune::TypeTree::TreeInfo<LFS>::leafCount) {
        }

        void update() {
          initializeLocalCoefficients();

          const auto& lfs = Base::localFunctionSpace();

          auto &subentityWiseDOFs = *lfs._subentityWiseDOFs_ptr;

          globalContainerIndices.resize(numberOfLeafs());
          for(auto& containerIndices: globalContainerIndices)
            containerIndices.clear();

          TypeTree::forEachLeafNode(lfs, [this, &lfs, &subentityWiseDOFs](auto &Node, auto &TreePath) {
            auto refEl = Dune::ReferenceElements<double, d>::general(Node.finiteElement().type());

            const auto leaf = Node.offsetLeafs;

            localDOFsOffset[leaf] = Node.offset;

            for (int c = 0; c < refEl.dimension + 1; ++c)
              for (int s = 0; s < refEl.size(c); ++s)
                // evaluate consecutive index of subentity
                lfs.gridFunctionSpace().ordering().mapIndex(subentityWiseDOFs[leaf].indexView(s, c),
                                                             globalContainerIndices[leaf].index(s, c));
          });
        }

        const CI &containerIndex(size_type leaf, size_type s, size_type c) const {
          return globalContainerIndices[leaf].index(s, c);
        }

        constexpr std::size_t numberOfLeafs() const {
          return Dune::TypeTree::TreeInfo<LFS>::leafCount;
        }

        constexpr std::size_t codims() const {
          return d + 1;
        }

        std::size_t subentities(std::size_t codim) const {
          // all leafs *should* have the same geometry, therefore we can use the first one
          return localCoefficients[0]->container[codim].size();
        }

        std::size_t sizeOfLocalDOFs(size_type leaf, size_type s, size_type c) const {
          return localCoefficients[leaf]->container[c][s].size();
        }

        std::size_t localIndex(size_type leaf, size_type s, size_type c, size_type i) const {
          return localCoefficients[leaf]->container[c][s][i] + localDOFsOffset[leaf];
        }

      private:

        void initializeLocalCoefficients() {
          TypeTree::forEachLeafNode(Base::localFunctionSpace(), [this](auto &Node, auto &TreePath) {
            const auto &fe = Node.finiteElement();
            inverseLocalCoefficientsMap.try_emplace(fe.size(), fe);
          });

          TypeTree::forEachLeafNode(Base::localFunctionSpace(), [this](auto &Node, auto &TreePath) {
            const auto &fe = Node.finiteElement();
            localCoefficients[Node.offsetLeafs] = &inverseLocalCoefficientsMap.at(fe.size());
          });
        }

        constexpr static std::size_t d = GFS::Traits::GridView::dimension;

        std::vector<SubentityWiseIndexWrapper<CI, d>> globalContainerIndices;
        std::map<std::size_t, InverseQkLocalCoefficients<d>> inverseLocalCoefficientsMap;
        std::vector<std::size_t> localDOFsOffset;
        std::vector<const InverseQkLocalCoefficients<d> *> localCoefficients;
      };
    }
  }
}


#endif //DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_LFSINDEXCACHE_HH
