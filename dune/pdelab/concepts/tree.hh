#ifndef DUNE_PDELAB_CONCEPT_TREE_HH
#define DUNE_PDELAB_CONCEPT_TREE_HH

#include <dune/pdelab/concepts/treenode.hh>
#include <dune/pdelab/concepts/multiindex.hh>
#include <dune/pdelab/common/tree_traversal.hh>

#include <concepts>
#include <type_traits>


namespace Dune::PDELab::inline Experimental::Concept {

  template<class T>
  concept Tree = requires(T tree)
  {
    // a tree is a traversable object where each node is a tree node
    Dune::PDELab::forEachNode(tree, []<TreeNode Node, MultiIndex Path>(Node&&, Path){});
  };

} // namespace Dune::PDELab::inline Experimental::Concept

#endif // DUNE_PDELAB_CONCEPT_TREE_HH
