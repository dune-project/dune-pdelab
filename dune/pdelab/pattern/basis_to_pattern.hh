#ifndef DUNE_PDELAB_PATTERN_BASIS_TO_PATTERN_HH
#define DUNE_PDELAB_PATTERN_BASIS_TO_PATTERN_HH

#include <dune/pdelab/concepts/treenode.hh>
#include <dune/pdelab/concepts/basis.hh>

#include <dune/pdelab/operator/local_assembly/interface.hh>

#include <dune/pdelab/pattern/local_sparsity_pattern.hh>

#include <dune/typetree/treecontainer.hh>

#include <dune/grid/common/mcmgmapper.hh>

namespace Dune::PDELab::inline Experimental {

template<class LocalOperator, class Pattern>
void basisToPattern(LocalOperator lop, Pattern& pattern) {

  using TestBasis = typename Pattern::RowSizeProvider;
  using TrialBasis = typename Pattern::ColSizeProvider;
  static_assert(Concept::Basis<TestBasis>);
  static_assert(Concept::Basis<TrialBasis>);

  using LocalTestBasis = typename TestBasis::LocalView;
  using LocalTrialBasis = typename TrialBasis::LocalView;

  using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<typename TrialBasis::EntitySet>;

  pattern.clear();

  TestBasis test{ pattern.rowSizeProvider() };
  TrialBasis trial{ pattern.colSizeProvider() };
  Mapper  _mapper{ trial.entitySet(), Dune::mcmgElementLayout() };

  LocalTrialBasis ltrial_in = trial.localView();
  LocalTrialBasis ltrial_out = trial.localView();
  LocalTestBasis ltest_in = test.localView();
  LocalTestBasis ltest_out = test.localView();

  const auto& es = trial.entitySet();

  auto it = es.template begin<0>();
  if (it == es.template end<0>())
    return;

  bind(*it, ltrial_in, ltest_in);
  bind(*it, ltrial_out, ltest_out);

  LocalSparsityPattern<Pattern> lpattern_ii{std::ref(pattern), ltest_in, ltrial_in};
  LocalSparsityPattern<Pattern> lpattern_io{std::ref(pattern), ltest_in, ltrial_out};
  LocalSparsityPattern<Pattern> lpattern_oi{std::ref(pattern), ltest_out, ltrial_in};
  LocalSparsityPattern<Pattern> lpattern_oo{std::ref(pattern), ltest_out, ltrial_out};

  unbind(ltrial_in, ltest_in);
  unbind(ltrial_out, ltest_out);

  for (const auto& entity : elements(trial.entitySet())) {
    if (LocalAssembly::skipEntity(lop, entity))
      continue;

    bind(entity, ltest_in, ltrial_in);

    if (LocalAssembly::doVolume(lop)) {
      LocalAssembly::patternVolume(lop, ltrial_in, ltest_in, lpattern_ii);
      lpattern_ii.commit(ltest_in, ltrial_in);
    }


    if (LocalAssembly::doSkeleton(lop) | LocalAssembly::doBoundary(lop)) {

      auto id_in = _mapper.index(entity);

      for (const auto& is : intersections(es, entity)) {

        if (LocalAssembly::skipIntersection(lop, is))
          continue;

        if (is.neighbor()) { // interior and periodic cases
          if (not LocalAssembly::doSkeleton(lop))
            continue;

          const auto& entity_out = is.outside();

        // visit face only the first time we see it
          auto id_out = _mapper.index(entity_out);
          // we have to be sure that outside entity was not skipped,
          // otherwise, neither side will be visited
          if (id_in < id_out and not LocalAssembly::skipEntity(lop, entity_out))
            continue;

          bind(is.outside(), ltest_out, ltrial_out);

          LocalAssembly::patternSkeleton(lop,
                          is,
                          ltrial_in,
                          ltest_in,
                          ltrial_out,
                          ltest_out,
                          lpattern_ii,
                          lpattern_io,
                          lpattern_oi,
                          lpattern_oo);

          lpattern_ii.commit(ltest_in, ltrial_in);
          lpattern_io.commit(ltest_in, ltrial_out);
          lpattern_oi.commit(ltest_out, ltrial_in);
          lpattern_oo.commit(ltest_out, ltrial_out);

          unbind(ltest_out, ltrial_out);
        } else if (is.boundary()) {
          LocalAssembly::patternBoundary(lop,
                          is,
                          ltrial_in,
                          ltest_in,
                          lpattern_ii);
          lpattern_ii.commit(ltest_in, ltrial_in);
        }

      }
    }

    unbind(ltest_in, ltrial_in);
  }

}

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_PATTERN_BASIS_TO_PATTERN_HH
