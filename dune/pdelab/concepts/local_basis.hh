#ifndef DUNE_PDELAB_CONCEPTS_LOCAL_BASIS_HH
#define DUNE_PDELAB_CONCEPTS_LOCAL_BASIS_HH

#include <dune/pdelab/concepts/tree.hh>
#include <dune/pdelab/concepts/treenode.hh>
#include <dune/pdelab/concepts/multiindex.hh>
#include <dune/pdelab/concepts/local_index_set.hh>

#include <dune/pdelab/common/tree_traversal.hh>

#include <dune/grid/concepts/entity.hh>

#include <concepts>
#include <ranges>

namespace Dune::PDELab::inline Experimental::Concept {

  template<class Leaf>
  concept LocalBasisLeaf = requires(const Leaf leaf, typename Leaf::size_type dof)
  {
    requires Dune::Concept::Entity<typename Leaf::Element>;
    { leaf.element() }          -> std::convertible_to<const typename Leaf::Element&>;
    { leaf.finiteElement() }    -> std::convertible_to<const typename Leaf::FiniteElement&>;
    // help to be more compatible with legacy PDELab local operators
    requires std::same_as<typename Leaf::Traits::FiniteElementType, typename Leaf::FiniteElement>;
    requires std::same_as<typename Leaf::Traits::FiniteElement, typename Leaf::FiniteElement>;
  };


  namespace Impl {
    void requireLocalBasisLeaf(LocalBasisLeaf auto& leaf);
  }

  template<class LST>
  concept LocalBasisTree = requires(LST lbasis_tree)
  {
    requires Tree<LST>;
    Dune::PDELab::forEachLeafNode(lbasis_tree, [](const auto& leaf){
      Impl::requireLocalBasisLeaf(leaf);
    });
  };

  template<class LS>
  concept LocalBasis = requires(LS lbasis, typename LS::Element entity)
  {
    requires Dune::Concept::Entity<typename LS::Element>;
    requires std::same_as<typename LS::GlobalBasis::EntitySet::template Codim<0>::Entity, typename LS::Element>;
    requires MultiIndex<typename LS::MultiIndex>;
    requires LocalIndexSet<LS>;
    requires LocalBasisTree<typename LS::Tree>;
    { lbasis.bind(entity)   } -> std::convertible_to<LS&>;
    { lbasis.unbind()       } -> std::convertible_to<LS&>;
    bind(entity, lbasis, lbasis);
    unbind(lbasis, lbasis);
    requires requires(const LS clbasis) {
      { clbasis.element() } -> std::convertible_to<const typename LS::Element&>;
      { clbasis.tree() } -> std::convertible_to<const typename LS::Tree&>;
    };
  };

} // namespace Dune::PDELab::inline Experimental::Concept

#endif // DUNE_PDELAB_CONCEPTS_LOCAL_BASIS_HH
