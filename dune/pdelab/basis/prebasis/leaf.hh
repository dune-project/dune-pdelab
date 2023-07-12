#ifndef DUNE_PDELAB_BASIS_PREBASIS_LEAF_HH
#define DUNE_PDELAB_BASIS_PREBASIS_LEAF_HH

#include <dune/pdelab/basis/prebasis/node.hh>
#include <dune/pdelab/basis/constraints/unconstrained.hh>

#include <dune/typetree/leafnode.hh>

#include <memory>
#include <utility>

namespace Dune::PDELab::inline Experimental {

/**
 * @brief Traits for leaf Space
 *
 * @tparam MS   A merging strategy
 * @tparam FEM  A Finite element map
 * @tparam CON  Constraint operator
 */
template<class MS, class FEM, class CO>
struct PreBasisLeafTraits
  : public PreBasisNodeTraits<MS>
{
  //! Finite element map
  using FiniteElementMap = FEM;
  //! Constraints
  using ConstraintsOperator = CO;
};

/**
 * @brief Leaf Space node
 *
 * @tparam MS   A merging strategy
 * @tparam FEM  A Finite element map
 * @tparam CON  Constraint
 */
template<class MergingStrategy, class FiniteElementMap, class ConstraintsOperator = Unconstrained>
class PreBasis
  : public TypeTree::LeafNode
  , public PreBasisNode<PreBasisLeafTraits<MergingStrategy,FiniteElementMap,ConstraintsOperator>>
{
  using BaseNode = PreBasisNode<PreBasisLeafTraits<MergingStrategy,FiniteElementMap,ConstraintsOperator>>;
  using TreeNode = TypeTree::LeafNode;

public:
  //! Node traits
  using typename BaseNode::Traits;

  /**
   * @brief Construct a new Leaf Discrete Function Space object
   * @note Use with value semantics
   *
   * @param merging_strategy  Node merging strategy
   * @param fem               Finite element mapper
   * @param constraints_op       ConstraintsOperator
   */
  PreBasis(const MergingStrategy& merging_strategy,
            std::shared_ptr<const FiniteElementMap> fem,
            const ConstraintsOperator& constraints_op = {})
    : BaseNode{ merging_strategy }
    , _finite_element_map{ std::move(fem) }
    , _constraints_op{ constraints_op }
  {
  }

  //! Copy constructor
  PreBasis(const PreBasis&) = default;

  template<std::same_as<void> = void>
  auto makeOrdering() const {
    return BaseNode::mergingStrategy().makeOrdering(*this);
  }

  template<std::same_as<void> = void>
  auto makeLocalOrdering() const {
    return BaseNode::mergingStrategy().makeLocalOrdering(*this);
  }

  //! Get finite element map
  const FiniteElementMap& finiteElementMap() const
  {
    return *_finite_element_map;
  }

  ConstraintsOperator constraintsOperator() const
  {
    return _constraints_op;
  }


  //! Compare if two leaf nodes are the same
  bool operator==(const PreBasis& other) const
  {
    return (other._finite_element_map == _finite_element_map);
  }

  //! Compare if two leaf nodes are different
  bool operator!=(const PreBasis& other) const
  {
    return not(*this == other);
  }

private:
  std::shared_ptr<const FiniteElementMap> _finite_element_map;
  ConstraintsOperator _constraints_op;
};

/**
 * @brief Makes a leaf discrete function space node
 *
 * @tparam MergingStrategy            Merging strategy
 * @tparam FiniteElementMapStorage    Shared pointer to a finite element map
 * @tparam ConstraintsOperator                ConstraintsOperator
 * @param merging_strategy            Merging strategy
 * @param fem                         Shared pointer to a finite element map
 * @param con                         ConstraintsOperator
 * @return auto                       The new discrete function space node
 */
template<class MergingStrategy,
         class FiniteElementMapStorage,
         class ConstraintsOperator>
auto
makePreBasis(const MergingStrategy& merging_strategy,
              FiniteElementMapStorage fem,
              const ConstraintsOperator& constraints_op)
{
  using FEM = std::remove_const_t<typename FiniteElementMapStorage::element_type>;
  using LeafDFS = PreBasis<MergingStrategy, FEM, ConstraintsOperator>;
  return LeafDFS{ merging_strategy, std::move(fem), constraints_op };
}

template<class MergingStrategy,
         class FiniteElementMapStorage>
auto
makePreBasis(const MergingStrategy& merging_strategy,
              FiniteElementMapStorage fem)
{
  return makePreBasis(merging_strategy, std::move(fem), Unconstrained{});
}

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_BASIS_PREBASIS_LEAF_HH
