#ifndef DUNE_PDELAB_CONCEPTS_LOCAL_INDEX_SET_HH
#define DUNE_PDELAB_CONCEPTS_LOCAL_INDEX_SET_HH

#include <dune/pdelab/concepts/treenode.hh>
#include <dune/pdelab/concepts/multiindex.hh>
#include <dune/pdelab/concepts/lockable.hh>

#include <dune/pdelab/common/tree_traversal.hh>
#include <dune/pdelab/common/partition/region.hh>

#include <concepts>

namespace Dune::PDELab::inline Experimental::Concept {

  template<class Leaf>
  concept LocalIndexSetLeaf = requires(const Leaf leaf, typename Leaf::size_type dof)
  {
    requires std::integral<typename Leaf::size_type>;
    requires MultiIndex<typename Leaf::MultiIndex>;
    { leaf.size() }           -> std::convertible_to<typename Leaf::size_type>;
    { leaf.path() }           -> std::convertible_to<typename Leaf::Path>;
    { leaf.localIndex(dof) }  -> std::convertible_to<typename Leaf::size_type>;
    { leaf.index(dof) }       -> std::convertible_to<typename Leaf::MultiIndex>;
  };

  namespace Impl {
    template<LocalIndexSetLeaf LeafNode>
    void requireLocalIndexSetLeaf(const LeafNode& leaf);
  }

  template<class LIS>
  concept LocalIndexSet = requires(const LIS lindex_set, typename LIS::size_type dof)
  {
    requires std::integral<typename LIS::size_type>;
    requires TreeNode<typename LIS::Tree>;
    requires MultiIndex<typename LIS::MultiIndex>;
    { lindex_set.size() }           -> std::convertible_to<typename LIS::size_type>;
    { lindex_set.maxSize() }        -> std::convertible_to<typename LIS::size_type>;
    { lindex_set.conforming() }     -> std::convertible_to<bool>;
    { lindex_set.partitionRegion() }-> std::convertible_to<EntitySetPartitioner::Region>;
    { lindex_set.tree() }           -> std::convertible_to<const typename LIS::Tree&>;
    { lindex_set.index(dof)  }      -> std::convertible_to<typename LIS::MultiIndex>;
    { lindex_set.globalBasis() }    -> std::convertible_to<const typename LIS::GlobalBasis&>;
    Dune::PDELab::forEachLeafNode(lindex_set.tree(), [](const auto& leaf){
      Impl::requireLocalIndexSetLeaf(leaf);
    });
    { Lockable<LIS> }       -> std::convertible_to<bool>;
  };

} // namespace Dune::PDELab::inline Experimental::Concept

#endif // DUNE_PDELAB_CONCEPTS_LOCAL_INDEX_SET_HH
