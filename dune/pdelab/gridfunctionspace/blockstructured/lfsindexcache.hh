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

      template<typename DOFLeafIterator,
          typename ContainerLeafIterator,
          std::size_t tree_depth>
      struct map_dof_indices_to_container_indices
          : public TypeTree::TreeVisitor
              , public TypeTree::DynamicTraversal
      {

        template<typename Ordering, typename TreePath>
        void leaf(const Ordering& ordering, TreePath tp)
        {
          ordering.map_lfs_indices(dof_leaf->cbegin(), dof_leaf->cend(), container_leaf->begin());
          dof_leaf++;
          container_leaf++;
        }

        template<typename Ordering, typename TreePath>
        void post(const Ordering& ordering, TreePath tp)
        {
          if (Ordering::consume_tree_index)
          {
            dof_tail_size--;
          }
          while (dof_stack.top() != dof_leaf) {
            using const_iterator = decltype(dof_stack.top()->cbegin());
            Dune::PDELab::DOFIndexViewIterator<const_iterator> dof_pos(dof_stack.top()->cbegin(), dof_tail_size);
            Dune::PDELab::DOFIndexViewIterator<const_iterator> dof_end(dof_stack.top()->cend(), dof_tail_size);
            ordering.map_lfs_indices(dof_pos, dof_end, container_stack.top()->begin());
            dof_stack.top()++;
            container_stack.top()++;
          }
          dof_stack.pop();
          container_stack.pop();
        }

        template<typename Ordering, typename TreePath>
        void pre(const Ordering& ordering, TreePath tp)
        {
          dof_stack.push(dof_leaf);
          container_stack.push(container_leaf);
          if (Ordering::consume_tree_index)
          {
            dof_tail_size++;
          }
        }

        map_dof_indices_to_container_indices(DOFLeafIterator dof_leaf_begin,
                                             ContainerLeafIterator container_leaf_begin,
                                             std::size_t dof_index_tail_length = 0)
            : dof_leaf(dof_leaf_begin), dof_tail_size(0)
            , container_leaf(container_leaf_begin)
        {}


        DOFLeafIterator dof_leaf;
        std::size_t dof_tail_size;
        ContainerLeafIterator container_leaf;
        std::stack<DOFLeafIterator,ReservedVector<DOFLeafIterator,tree_depth> > dof_stack;
        std::stack<ContainerLeafIterator,ReservedVector<ContainerLeafIterator,tree_depth> > container_stack;

      };

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

          auto &subentityWiseDOFs = *lfs.subentityWiseDOFs_ptr;

          globalContainerIndices.resize(numberOfLeafs());

          TypeTree::forEachLeafNode(lfs, [this, &lfs, &subentityWiseDOFs](auto &Node, auto &TreePath) {
            const auto leaf = Node.offsetLeafs;
            localDOFsOffset[leaf] = Node.offset;

            globalContainerIndices[leaf].clear();
            globalContainerIndices[leaf].setup(Node.gridFunctionSpace().finiteElementMap(), Node.finiteElement().type());
          });

          map_dof_indices_to_container_indices<
              typename LFS::DOFIndexSubentityWiseContainer::const_iterator,
              typename decltype(globalContainerIndices)::iterator,
              TypeTree::TreeInfo<Ordering>::depth
          > index_mapper(subentityWiseDOFs.begin(), globalContainerIndices.begin(), lfs.subSpaceDepth());
          TypeTree::applyToTree(lfs.gridFunctionSpace().ordering(),index_mapper);
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
          return localCoefficients[0]->size(codim);
        }

        std::size_t sizeOfLocalDOFs(size_type leaf, size_type s, size_type c) const {
          return localCoefficients[leaf]->size(s, c);
        }

        std::size_t localIndex(size_type leaf, size_type s, size_type c, size_type i) const {
          return localCoefficients[leaf]->localDOF(Dune::LocalKey(s, c, i)) + localDOFsOffset[leaf];
        }

      private:

        void initializeLocalCoefficients() {
          TypeTree::forEachLeafNode(Base::localFunctionSpace(), [this](auto &Node, auto &TreePath) {
            const auto &fe = Node.finiteElement();
            if(inverseLocalCoefficientsMap.find(fe.size()) == inverseLocalCoefficientsMap.end())
              inverseLocalCoefficientsMap.emplace(fe.size(), InverseQkLocalCoefficients<d>(fe));
          });

          TypeTree::forEachLeafNode(Base::localFunctionSpace(), [this](auto &Node, auto &TreePath) {
            const auto &fe = Node.finiteElement();
            localCoefficients[Node.offsetLeafs] = &inverseLocalCoefficientsMap.at(fe.size());
          });
        }

        constexpr static std::size_t d = GFS::Traits::GridView::dimension;

        std::vector<SubentityIndexWrapper<CI, d>> globalContainerIndices;
        std::map<std::size_t, InverseQkLocalCoefficients<d>> inverseLocalCoefficientsMap;
        std::vector<std::size_t> localDOFsOffset;
        std::vector<const InverseQkLocalCoefficients<d> *> localCoefficients;
      };
    }
  }
}


#endif //DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_LFSINDEXCACHE_HH
