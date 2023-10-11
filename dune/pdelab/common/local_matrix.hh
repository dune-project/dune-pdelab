#ifndef DUNE_PDELAB_COMMON_LOCAL_MATRIX_HH
#define DUNE_PDELAB_COMMON_LOCAL_MATRIX_HH

#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/concepts/local_index_set.hh>
#include <dune/pdelab/concepts/container.hh>

#include <dune/pdelab/common/local_container_entry.hh>
#include <dune/pdelab/common/container_entry.hh>

#include <dune/typetree/treecontainer.hh>

#include <functional>
#include <mutex>
#include <utility>
#include <vector>

namespace Dune::PDELab::inline Experimental {


// TODO move to concepts folder
template<class T>
concept Coefficient = requires(T v)
{
  { v += v } -> std::convertible_to<T>;
  requires std::constructible_from<T,double>;
  requires std::copyable<T>;
};

template<Concept::Basis _TestBasis, Concept::Basis _TrialBasis, class Container>
class LocalMatrixView
{
public:
  using TestBasis = _TestBasis;
  using TrialBasis = _TrialBasis;
  using Weight = int;

  LocalMatrixView(const auto&, const auto&, Container& container)
    : _container{ std::ref(container) }
  {
  }

  void accumulate(const Concept::TreeNode auto& test_node, std::size_t test_dof,
                  const Concept::TreeNode auto& trial_node, std::size_t trial_dof,
                  auto val)
  {
    localContainerEntry(_container.get(), test_node, test_dof, trial_node, trial_dof) += val;
  }

  void fetch_add(Concept::LocalBasis auto& ltest, Concept::LocalBasis auto& ltrial) {}

  void clear(const Concept::LocalBasis auto& ltest, const Concept::LocalBasis auto& ltrial) {}

  [[nodiscard]] static constexpr Weight weight() { return 1; }

private:

  std::reference_wrapper<Container> _container;
};


template<Concept::Basis _TestBasis, Concept::Basis _TrialBasis, class Container>
class LocalMatrixBuffer
{
  static auto makeContainer(const Concept::LocalBasis auto& ltest, const Concept::LocalBasis auto& ltrial) {
    return TypeTree::makeTreeContainer(ltest.tree(), [&](const auto& test_leaf){
      return TypeTree::makeTreeContainer(ltrial.tree(), [&](const auto& trial_leaf){
        using value_type = std::remove_cvref_t<decltype(localContainerEntry(std::declval<const Container&>(), test_leaf, 0, trial_leaf, 0))>;
        static_assert(Coefficient<value_type>, "Container does not provide types that fullfil the coefficient requirement");
        return std::vector<value_type>{};
      });
    });
  }

  using LocalTreeContainer = decltype(makeContainer(std::declval<_TestBasis>().localView(), std::declval<_TrialBasis>().localView()));

public:

  using TestBasis = _TestBasis;
  using TrialBasis = _TrialBasis;
  using Weight = int;

  LocalMatrixBuffer(const TestBasis& test_basis, const TrialBasis& trial_basis, Container& container)
    : _container{ std::ref(container) }
    , _ltest_constraints{ test_basis.localConstraints() }
    , _ltrial_constraints{ trial_basis.localConstraints() }
  {
    auto it = test_basis.entitySet().template begin<0>();
    if (test_basis.entitySet().size(0) == 0)
      DUNE_THROW(Dune::InvalidStateException, "This function does not work with empty entity sets");

    auto ltest = test_basis.localView();
    auto ltrial = trial_basis.localView();
    bind(*it, ltest, ltrial);
    _lcontainer = makeContainer(ltest, ltrial);
  }

  LocalMatrixBuffer(const Concept::LocalBasis auto& ltest, const Concept::LocalBasis auto& ltrial, Container& container)
    : _lcontainer{ makeContainer(ltest, ltrial) }
    , _container{ std::ref(container) }
    , _ltest_constraints{ ltest.globalBasis().localConstraints() }
    , _ltrial_constraints{ ltrial.globalBasis().localConstraints() }
  {}

  void accumulate(const Concept::TreeNode auto& test_node, std::size_t test_dof,
                  const Concept::TreeNode auto& trial_node, std::size_t trial_dof,
                  auto val)
  {
    _lcontainer[test_node.path()][trial_node.path()][trial_node.size()*test_dof + trial_dof] += val;
  }

  template<Concept::LocalBasis LocalBasisTest, Concept::LocalBasis LocalBasisTrial>
  void fetch_add(LocalBasisTest& ltest, const LocalBasisTrial& ltrial) {
    _ltest_constraints.bind(ltest.element());
    _ltrial_constraints.bind(ltrial.element());

    {
      [[maybe_unused]] auto scope_guard = [&](){
        if constexpr (Concept::Lockable<LocalBasisTest>) {
          if (ltest.partitionRegion() == EntitySetPartitioner::shared_region)
            return std::unique_lock{ltest};
          else
            return std::unique_lock{ltest, std::defer_lock};
        } else {
          return nullptr;
        }
      }();

      forEachLeafNode(ltest.tree(), [&](const auto& ltest_node, auto test_path) {
        const auto& ltest_constrain = containerEntry(_ltest_constraints.tree(), test_path);
        for (std::size_t test_dof = 0; test_dof != ltest_node.size(); ++test_dof) {
          if (ltest_constrain.isConstrained(test_dof)) {
            if (ltest_constrain.linearCoefficients(test_dof).empty()) {
              // Dirichlet constraint case:
              // Technically, we need to set 1. on the trial DOF corresponding to the test DOF.
              // However, note that trees may be different (e.g. operator splitting), how to constrain the right space?
              // TODO the following is only valid if the trial tree is the same as the test tree
              // the problem is that we don't know how to identify the "diagonal" entry!
              const auto& ltrial_node = containerEntry(ltrial.tree(), test_path);
              const auto& ltrial_constrain = containerEntry(_ltrial_constraints.tree(), test_path);
              assert(ltrial_node.size() == ltest_node.size());
              assert(ltrial_constrain.isConstrained(test_dof));
              assert(ltrial_constrain.linearCoefficients(test_dof).empty());
              localContainerEntry(_container.get(), ltest_node, test_dof, ltrial_node, test_dof) = 1.;
            } else {
              DUNE_THROW(NotImplemented, "Hanging nodes not implemented yet");
            }
          } else {
            forEachLeafNode(ltrial.tree(), [&](const auto& ltrial_node, auto trial_path) {
              for (std::size_t trial_dof = 0; trial_dof != ltrial_node.size(); ++trial_dof) {
                const auto& val = _lcontainer[ltest_node.path()][ltrial_node.path()][ltrial_node.size()*test_dof + trial_dof];
                using value_type = std::remove_cvref_t<decltype(val)>;
                if (val != value_type{0.})
                  localContainerEntry(_container.get(), ltest_node, test_dof, ltrial_node, trial_dof) += val;
              }
            });
          }
        }
      });
    }

    _ltest_constraints.unbind();
    _ltrial_constraints.unbind();
  }

  void clear(const Concept::LocalBasis auto& ltest, const Concept::LocalBasis auto& ltrial) {
    forEachLeafNode(ltest.tree(), [&](const auto& ltest_node, auto) {
      forEachLeafNode(ltrial.tree(), [&](const auto& ltrial_node, auto) {
        auto& container = _lcontainer[ltest_node.path()][ltrial_node.path()];
        using value_type = std::remove_cvref_t<decltype(container[0])>;
        container.assign(ltest_node.size()*ltrial_node.size(), value_type{0.});
      });
    });
  }

  [[nodiscard]] static constexpr Weight weight() { return 1; }

private:

  LocalTreeContainer _lcontainer;
  std::reference_wrapper<Container> _container;
  typename TestBasis::LocalConstraints _ltest_constraints;
  typename TrialBasis::LocalConstraints _ltrial_constraints;
};


template<Concept::LocalMutableMatrix LocalMatrix, class WeightType>
class WeightedLocalMatrixView {
public:
  using TestBasis = typename LocalMatrix::TestBasis;
  using TrialBasis = typename LocalMatrix::TrialBasis;
  using Weight = decltype(typename LocalMatrix::Weight{} * WeightType{});

  WeightedLocalMatrixView(LocalMatrix& lmatrix, WeightType weight)
    : _weight{weight}
    , _lmatrix{lmatrix}
  {}

  void accumulate(const Concept::LeafTreeNode auto& ltest_node, auto test_dof, const Concept::LeafTreeNode auto& ltrial_node, auto trial_dof, auto val) {
    _lmatrix.accumulate(ltest_node, test_dof, ltrial_node, trial_dof, weight() * val);
  }

  [[nodiscard]] const auto& operator()(const Concept::LeafTreeNode auto& ltest_node, auto test_dof, const Concept::LeafTreeNode auto& ltrial_node, auto trial_dof) const {
    return _lmatrix(ltest_node, test_dof, ltrial_node, trial_dof);
  }

  [[nodiscard]] auto& operator()(const Concept::LeafTreeNode auto& ltest_node, auto test_dof, const Concept::LeafTreeNode auto& ltrial_node, auto trial_dof) {
    return _lmatrix(ltest_node, test_dof, ltrial_node, trial_dof);
  }

  [[nodiscard]] Weight weight() const {
    return _lmatrix.weight() * _weight;
  }

private:
  WeightType _weight;
  LocalMatrix& _lmatrix;
};


// applies matrix-vector product on the fly
// TODO: make this a proper local matrix interface
template<class TestContaienr, class TrialContainer>
struct LocalMatrixApplyView {

  using TestBasis = typename TestContaienr::Basis;
  using TrialBasis = typename TrialContainer::Basis;

  void accumulate(const auto& ltest, auto test_dof, const auto& ltrial, auto trial_dof, auto value) {
    ltest_container.accumulate(ltest, test_dof, ltria_container(ltrial, trial_dof) * value);
  }

  // [[nodiscard]] static constexpr Weight weight() { return ltest_container.weight(); ??? }

  TestContaienr& ltest_container;
  const TrialContainer& ltria_container;
};



} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_LOCAL_MATRIX_HH
