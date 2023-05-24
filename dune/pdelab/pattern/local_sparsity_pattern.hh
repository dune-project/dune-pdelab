#ifndef DUNE_PDELAB_PATTERN_LOCAL_SPARSITY_PATTERN_HH
#define DUNE_PDELAB_PATTERN_LOCAL_SPARSITY_PATTERN_HH

#include <dune/pdelab/concepts/basis.hh>

#include <dune/pdelab/common/tree_traversal.hh>

#include <dune/typetree/treecontainer.hh>

namespace Dune::PDELab::inline Experimental {

// add links between two local basis degrees of freedom
template<class Pattern>
class LocalSparsityPattern {

  using TestBasis = typename Pattern::RowSizeProvider;
  using TrialBasis = typename Pattern::ColSizeProvider;
  static_assert(Concept::Basis<TestBasis>);
  static_assert(Concept::Basis<TrialBasis>);

  using LocalTestBasis = typename TestBasis::LocalView;
  using LocalTrialBasis = typename TrialBasis::LocalView;

  using TestSizeType = typename LocalTestBasis::size_type;
  using TrialSizeType = typename LocalTrialBasis::size_type;


  static auto makeLinkContainer(const LocalTestBasis& ltest, const LocalTrialBasis& ltrial) {
    return TypeTree::makeTreeContainer(ltest.tree(), [&](const auto&){
      return TypeTree::makeTreeContainer(ltrial.tree(), [&](const auto&){
        return std::vector<std::tuple<TestSizeType,TrialSizeType>>{};
      });
    });
  }

  using LinkContainer = decltype(makeLinkContainer(std::declval<LocalTestBasis>(), std::declval<LocalTrialBasis>()));

  static auto makeLeafBasisPtrs(const auto& lbasis) {
    return TypeTree::makeTreeContainer(lbasis.tree(), [&](const auto& leaf){
      return &leaf;
    });
  }


  using LocalTestPtrs = decltype(makeLeafBasisPtrs(std::declval<LocalTestBasis>()));
  using LocalTrialPtrs = decltype(makeLeafBasisPtrs(std::declval<LocalTrialBasis>()));


public:
  LocalSparsityPattern(std::reference_wrapper<Pattern> pattern, const LocalTestBasis& ltest, const LocalTrialBasis& ltrial)
    : _links{makeLinkContainer(ltest, ltrial)}
    , _ltest_ptr{makeLeafBasisPtrs(ltest)}
    , _ltrial_ptr{makeLeafBasisPtrs(ltrial)}
    , _pattern{ pattern }
  {
    PDELab::forEachLeafNode(ltest.tree(), [&](auto& leaf_ltest, auto test_path){
      PDELab::forEachLeafNode(ltrial.tree(), [&](auto& leaf_ltrial, auto trial_path){
        _ltest_ptr[test_path] = nullptr;
        _ltrial_ptr[trial_path] = nullptr;
      });
    });
  }

  void addLink(
      const Concept::LeafTreeNode auto& leaf_ltest, TestSizeType test_dof,
      const Concept::LeafTreeNode auto& leaf_ltrial, TrialSizeType trial_dof)
  {
#ifndef NDEBUG
    assert(test_dof < leaf_ltest.size());
    assert(trial_dof < leaf_ltrial.size());
    auto old_ltest_ptr = std::exchange(_ltest_ptr[leaf_ltest.path()], &leaf_ltest);
    auto old_ltrial_ptr = std::exchange(_ltrial_ptr[leaf_ltrial.path()], &leaf_ltrial);
    assert(old_ltest_ptr == nullptr or old_ltest_ptr == &leaf_ltest);
    assert(old_ltrial_ptr == nullptr or old_ltrial_ptr == &leaf_ltrial);
#endif
    _links[leaf_ltest.path()][leaf_ltrial.path()].emplace_back(test_dof, trial_dof);
  }

  // transform local links to global links
  void commit(LocalTestBasis& ltest, LocalTrialBasis& ltrial) {

    PDELab::forEachLeafNode(ltest.tree(), [&](auto& leaf_ltest, auto test_path){
      PDELab::forEachLeafNode(ltrial.tree(), [&](auto& leaf_ltrial, auto trial_path){

        assert((_ltest_ptr[test_path] == &leaf_ltest or _ltest_ptr[test_path] == nullptr) &&
          "The test basis used to add a pattern entry is incorrect");
        assert((_ltrial_ptr[trial_path] == &leaf_ltrial or _ltrial_ptr[trial_path] == nullptr )&&
          "The trial basis used to add a pattern entry is incorrect");

        for (auto [test_dof, trial_dof] : _links[test_path][trial_path] ) {

          // TODO: fix pattern constraints when basis constraints get implemented
          bool test_constrained = false /*not empty(leaf_ltest.constraints(test_dof))*/;
          bool trial_constrained = false /*not empty(leaf_ltrial.constraints(trial_dof))*/;

          if (not (trial_constrained or test_constrained)) {
            // when no constraints are applied, we honor what local operator asked for
            assert(test_dof < leaf_ltest.size());
            assert(trial_dof < leaf_ltrial.size());
            _pattern.get().addLink(leaf_ltest.index(test_dof), leaf_ltrial.index(trial_dof));
          } else {
            DUNE_THROW(NotImplemented, "");
          //   // otherwise, we make a tensor product of constraints

          //   // in the following, you may understand the unconstrained and
          //   // dirichlet cases to be the same as weighted constrains to self
          //   // dof (only one constraint) with a unit weight
          //   // (not important for pattern). Then, the links will form a tensor
          //   // product of weights over the pattern.
          //   _test_clinks.clear();
          //   _trial_clinks.clear();

          //   // first fill rows
          //   if (not empty(leaf_ltest.constraints(test_dof)))
          //     for (auto [test_cdof, w] : leaf_ltest.constraints(test_dof))
          //       _test_clinks.push_back(test_cdof);
          //   else // for both dirichlet and unconstrained dofs
          //     _test_clinks.push_back(leaf_ltest.index(test_dof));

          //   // then fill columns
          //   if (not empty(leaf_ltrial.constraints(trial_dof)))
          //     for (auto [trial_cdof, w] : leaf_ltrial.constraints(trial_dof))
          //       _trial_clinks.push_back(trial_cdof);
          //   else // for both dirichlet and unconstrained dofs
          //     _trial_clinks.push_back(leaf_ltrial.index(trial_dof));

          //   // finally, add links on the tensor product of rows and cols
          //   for (auto test_clink : _test_clinks)
          //     for (auto trial_clink : _trial_clinks)
          //       _pattern.get().addLink(test_clink, trial_clink);
          }
        }
        _links[test_path][trial_path].clear();
#ifndef NDEBUG
        _ltest_ptr[test_path] = nullptr;
        _ltrial_ptr[trial_path] = nullptr;
#endif
      });
    });

  }

private:
  LinkContainer _links; // regular links from local pattern
  LocalTestPtrs _ltest_ptr;
  LocalTrialPtrs _ltrial_ptr;
  std::vector<typename LocalTestBasis::MultiIndex> _test_clinks; // row constrain links
  std::vector<typename LocalTrialBasis::MultiIndex> _trial_clinks; // col constrain links
  std::reference_wrapper<Pattern> _pattern;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_PATTERN_LOCAL_SPARSITY_PATTERN_HH
