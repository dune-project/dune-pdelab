#ifndef DUNE_PDELAB_BASIS_PREBASIS_NODE_HH
#define DUNE_PDELAB_BASIS_PREBASIS_NODE_HH

#include <string>
#include <utility>

namespace Dune::PDELab::inline Experimental {

/**
 * @brief Traits for Space nodes
 *
 * @tparam MS  a merging strategy
 */
template<class MS>
struct PreBasisNodeTraits
{
  using MergingStrategy = MS;
};

/**
 * @brief Base node for Space
 *
 * @tparam NodeTraits Traits for this node
 */
template<class NodeTraits>
class PreBasisNode
{
public:
  using Traits = NodeTraits;

  /**
   * @brief Construct a new Discrete Function Space Node object
   *
   * @param merging_strategy  rvalue or lvalue merging strategy
   */
  template<class MergingStrategy>
  PreBasisNode(MergingStrategy&& merging_strategy)
    : _merging_strategy{ std::forward<MergingStrategy>(merging_strategy) }
  {
  }

  //! Get name of the discrete function space node
  const std::string& name() const { return _name; }

  //! Set a name for this node
  void name(const std::string& name) { _name = name; }

  //! Get merging strategy
  typename Traits::MergingStrategy& mergingStrategy()
  {
    return _merging_strategy;
  }

  //! Get merging strategy
  const typename Traits::MergingStrategy& mergingStrategy() const
  {
    return _merging_strategy;
  }

private:
  std::string _name;
  typename Traits::MergingStrategy _merging_strategy;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_BASIS_PREBASIS_NODE_HH
