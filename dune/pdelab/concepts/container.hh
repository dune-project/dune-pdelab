#ifndef DUNE_PDELAB_CONCEPT_CONTAINER_HH
#define DUNE_PDELAB_CONCEPT_CONTAINER_HH

#include <dune/pdelab/concepts/basis.hh>

#include <dune/pdelab/common/tree_traversal.hh>
#include <dune/pdelab/common/container_entry.hh>
#include <dune/pdelab/common/local_container_entry.hh>

#include <concepts>
#include <type_traits>


namespace Dune::PDELab::inline Experimental::Concept {

  namespace Impl {

    template<class C, class Leaf>
    concept ContainerLeaf = requires(C container, const Leaf leaf, const typename Leaf::size_type dof)
    {
      { Dune::PDELab::localContainerEntry(container, leaf, dof) } -> std::same_as<decltype(Dune::PDELab::containerEntry(container, leaf.index(dof)))>;
    };


    template<class C, class Leaf>
    requires ContainerLeaf<C, Leaf>
    void requireContainerLeaf(const Leaf& leaf);

    template<class LC, class Leaf>
    concept LocalConstContainerLeaf = requires(const LC lcontainer, const Leaf leaf, const typename Leaf::size_type dof)
    {
      lcontainer(leaf, dof); // access value
    };

    template<class LC, class Leaf>
    requires LocalConstContainerLeaf<LC, Leaf>
    void requireLocalConstContainerLeaf(const Leaf& leaf);

    template<class LC, class Leaf>
    concept LocalMutableContainerLeaf = requires(LC lcontainer, const Leaf leaf, const typename Leaf::size_type dof)
    {
      { lcontainer.weight() } -> std::convertible_to<typename LC::Weight>;  // access weight
      lcontainer.accumulate(leaf, dof, lcontainer(leaf, dof));              // weighted accumulation
      lcontainer(leaf, dof) += lcontainer.weight() * lcontainer(leaf, dof); // raw weighted accumulation
    };

    template<class LC, class Leaf>
    requires LocalMutableContainerLeaf<LC, Leaf>
    void requireLocalMutableContainerLeaf(const Leaf& leaf);



    template<class LM, class LeafTest, class LeafTrial>
    concept LocalConstMatrixLeaf = requires(const LM lmatrix, const LeafTest ltest, const typename LeafTest::size_type test_dof, const LeafTrial ltrial, const typename LeafTrial::size_type trial_dof)
    {
      lmatrix(ltest, test_dof, ltrial, trial_dof); // access value
    };

    template<class LM, class LeafTest, class LeafTrial>
    requires LocalConstMatrixLeaf<LM, LeafTest, LeafTrial>
    void requireLocalConstMatrixLeaf(const LeafTest& ltest, const LeafTrial& ltrial);

    template<class LM, class LeafTest, class LeafTrial>
    concept LocalMutableMatrixLeaf = requires(const LM lmatrix, const LeafTest ltest, const typename LeafTest::size_type test_dof, const LeafTrial ltrial, const typename LeafTrial::size_type trial_dof)
    {
      { lmatrix.weight() } -> std::convertible_to<typename LM::Weight>;  // access weight
      lmatrix.accumulate(ltest, test_dof, ltrial, trial_dof, lmatrix(ltest, test_dof, ltrial, trial_dof));              // weighted accumulation
      lmatrix(ltest, test_dof, ltrial, trial_dof) += lmatrix.weight() * lmatrix(ltest, test_dof, ltrial, trial_dof); // raw weighted accumulation
    };

    template<class LM, class LeafTest, class LeafTrial>
    requires LocalMutableMatrixLeaf<LM, LeafTest, LeafTrial>
    void requireLocalMutableMatrixLeaf(const LeafTest& ltest, const LeafTrial& ltrial);

  } // namespace Impl

  template<class C, class S>
  concept Container = requires(C container, S space)
  {
    requires Basis<S>;
    Dune::PDELab::forEachLeafNode(space.localView().tree(), [](const auto& leaf){
      Impl::requireContainerLeaf<C>(leaf);
    });
  };

  template<class LC>
  concept LocalConstContainer = requires(LC lcontainer, typename LC::Basis space)
  {
    requires Basis<typename LC::Basis>;
    Dune::PDELab::forEachLeafNode(space.localView().tree(), [](const auto& leaf){
      Impl::requireLocalConstContainerLeaf<LC>(leaf);
    });
  };

  template<class LC>
  concept LocalMutableContainer = requires(LC lcontainer, typename LC::Basis space)
  {
    requires LocalConstContainer<LC>;
    Dune::PDELab::forEachLeafNode(space.localView().tree(), [](const auto& leaf){
      Impl::requireLocalMutableContainerLeaf<LC>(leaf);
    });
  };


  template<class LM>
  concept LocalConstMatrix = requires(typename LM::TestBasis test, typename LM::TrialBasis trial)
  {
    requires Basis<typename LM::TestBasis>;
    requires Basis<typename LM::TrialBasis>;
    Dune::PDELab::forEachLeafNode(test.localView().tree(), [](const auto& ltest){
      // Dune::PDELab::forEachLeafNode(trial.localView().tree(), [](const auto& ltrial){
      //   Impl::requireLocalConstMatrixLeaf<LM>(ltest, ltrial);
      // });
    });
  };

  template<class LM>
  concept LocalMutableMatrix = requires(typename LM::TestBasis test, typename LM::TrialBasis trial)
  {
    requires LocalConstMatrix<LM>;
    Dune::PDELab::forEachLeafNode(test.localView().tree(), [](const auto& ltest){
      // Dune::PDELab::forEachLeafNode(trial.localView().tree(), [](const auto& ltrial){
      //   Impl::requireLocalMutableMatrixLeaf<LM>(ltest, ltrial);
      // });
    });
  };

} // namespace Dune::PDELab::inline Experimental::Concept

#endif // DUNE_PDELAB_CONCEPT_CONTAINER_HH
