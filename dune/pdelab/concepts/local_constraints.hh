#ifndef DUNE_PDELAB_CONCEPTS_LOCAL_CONSTRAINTS_HH
#define DUNE_PDELAB_CONCEPTS_LOCAL_CONSTRAINTS_HH

#include <dune/pdelab/concepts/tree.hh>
#include <dune/pdelab/concepts/treenode.hh>
#include <dune/pdelab/concepts/multiindex.hh>

#include <dune/pdelab/common/tree_traversal.hh>

#include <dune/grid/concepts/entity.hh>

#include <concepts>
#include <ranges>

namespace Dune::PDELab::inline Experimental::Concept {

  namespace Impl {
    template<class Range>
    concept LocalLinearConstraintsRange = requires {
      requires std::ranges::range<Range>;
      requires MultiIndex<typename std::ranges::range_value_t<Range>::first_type>;
      requires std::is_arithmetic_v<typename std::ranges::range_value_t<Range>::second_type>;
    };


    template<class Leaf>
    concept LocalConstraintsLeaf = requires(Leaf leaf, typename Leaf::size_type dof)
    {
      { leaf.isConstrained(dof) }         -> std::convertible_to<bool>;
      { leaf.translationValue(dof) }; //-> std::convertible_to<typename Leaf::Value>; // TODO change name of "Vaule" ??
      { leaf.linearCoefficients(dof) }      -> LocalLinearConstraintsRange;
    };

    void requireLocalConstraintsLeaf(LocalConstraintsLeaf auto& leaf);
  }

  template<class LCT>
  concept LocalConstraintsTree = requires(LCT lconstraints_tree)
  {
    requires Tree<LCT>;
    Dune::PDELab::forEachLeafNode(lconstraints_tree, [](const auto& leaf){
      Impl::requireLocalConstraintsLeaf(leaf);
    });
  };

  template<class LC>
  concept LocalConstraints = requires(LC lconstraints)
  {
    // requires Dune::Concept::Entity<typename LS::Element>;
    requires LocalConstraintsTree<typename LC::Tree>;
    // lbasis.bind(entity);
    lconstraints.unbind();
    // bind(entity, lbasis, lbasis);
    // unbind(lbasis, lbasis);
    { lconstraints.tree() } -> std::convertible_to<const typename LC::Tree&>;
  };

} // namespace Dune::PDELab::inline Experimental::Concept

#endif // DUNE_PDELAB_CONCEPTS_LOCAL_CONSTRAINTS_HH
