#ifndef DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_ASSEMBLER_FORWARD_HH
#define DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_ASSEMBLER_FORWARD_HH

#include <dune/pdelab/operator/forward/instationary/coefficients.hh>
#include <dune/pdelab/operator/forward/instationary/traits.hh>

#include <dune/pdelab/operator/local_assembly/interface.hh>

#include <dune/pdelab/common/local_container.hh>
#include <dune/pdelab/common/tree_traversal.hh>
// #include <dune/pdelab/common/entityset.hh>

#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/concepts/container.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/common/float_cmp.hh>

#include <vector>
#include <ranges>

namespace Dune::PDELab::inline Experimental {

/**
 * @brief Residual of a MassStiffness operator using PDELab local operators
 *
 * @tparam Coefficients
 * @tparam Residual
 * @tparam TrialBasis
 * @tparam TestBasis
 * @tparam MassLocalOperator
 * @tparam StiffnessLocalOperator
 * @tparam TimePoint
 * @tparam Duration
 * @tparam dt_position
 */
template<class          Coefficients,
         class          Residual,
         Concept::Basis TrialBasis,
         Concept::Basis TestBasis,
         class          MassLocalOperator,
         class          StiffnessLocalOperator,
         class TimePoint = double,
         class Duration  = double,
         DurationPosition dt_position = DurationPosition::StiffnessNumerator>
class InstationaryForwardAssembler : public Operator<Coefficients, Residual>
{
  using Base = Operator<Coefficients, Residual>;

  using StageCoefficients = typename std::ranges::range_value_t<Coefficients>;
  using StageResidual     = typename std::ranges::range_value_t<Residual>;

  static_assert(Concept::Container<StageCoefficients,TrialBasis>);
  static_assert(Concept::Container<StageResidual,    TestBasis>);

  static_assert(std::is_copy_constructible_v<MassLocalOperator>);
  static_assert(std::is_copy_constructible_v<StiffnessLocalOperator>);

  using LocalTestBasis  = typename TestBasis::LocalView;
  using LocalTrialBasis = typename TrialBasis::LocalView;

  using MassFactor      = typename InstationaryTraits<dt_position>::MassFactor<Duration>;
  using StiffnessFactor = typename InstationaryTraits<dt_position>::StiffnessFactor<Duration>;

  using LocalResidual           = LocalContainerBuffer<TestBasis, StageResidual>;
  using LocalCoefficients       = LocalContainerBuffer<TrialBasis, const StageCoefficients>;
  using LocalMassResidual       = WeightedLocalContainerView<LocalResidual, MassFactor>;
  using LocalStiffnessResidual  = WeightedLocalContainerView<LocalResidual, StiffnessFactor>;

  using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<typename TrialBasis::EntitySet>;

public:
  InstationaryForwardAssembler(const TrialBasis& trial,
                        const TestBasis& test,
                        const MassLocalOperator& mass_lop,
                        const StiffnessLocalOperator& stiff_lop)
    : _trial{ trial }
    , _test{ test }
    , _mass_lop{ mass_lop }
    , _stiff_lop{ stiff_lop }
    , _mapper{ _trial.entitySet(), Dune::mcmgElementLayout() }
  {
  }

  ErrorCondition apply(const Coefficients& coefficients, Residual& residuals) override
  {
    TRACE_EVENT("dune", "MassStiffness::Residual");
    bool static_dispatch_done = false;
    // unroll the loop for small sizes (important for small local operators)
    using namespace Dune::Indices;
    constexpr auto res_unfold = _3;
    constexpr auto coeff_unfold = _3;
    Dune::Hybrid::forEach(Dune::range(_1, res_unfold), [&](auto stages) {
      if (residuals.size() == stages) {

        Dune::Hybrid::forEach(Dune::range(_1, coeff_unfold), [&](auto steps) {
          if (coefficients.size() == steps) {
            applyImpl(coefficients, steps, residuals, stages);
            static_dispatch_done = true;
          }
        });

        if (not static_dispatch_done)
          applyImpl(coefficients, coefficients.size(), residuals, stages);

        static_dispatch_done = true;
      }
    });

    if (not static_dispatch_done)
      applyImpl(coefficients, coefficients.size(), residuals, residuals.size());
    return {};
  }

private:

  void applyImpl(const Coefficients& coefficients, auto steps, Residual& residuals, auto stages) const
  {
    LocalTrialBasis ltrial_in  = _trial.localView();
    LocalTrialBasis ltrial_out = _trial.localView();

    LocalTestBasis ltest_in  = _test.localView();
    LocalTestBasis ltest_out = _test.localView();

    assert(coefficients.size() == steps);
    assert(residuals.size() == stages);

    const InstationaryCoefficients& icoeff = getInstationaryCoefficients();
    const TimePoint& time_point = getTimePoint();
    const Duration& duration = getDuration();

    std::vector<std::vector<bool>> do_mass(stages);
    std::vector<std::vector<bool>> do_stiff(stages);
    for (std::size_t stage = 0; stage != stages; ++stage) {
      do_mass[stage].resize(steps);
      do_stiff[stage].resize(steps);
      for (std::size_t step = 0; step != steps; ++step) {
        do_mass[stage][step] = icoeff.doMass(stage, step);
        do_stiff[stage][step] = icoeff.doStiffness(stage, step);
      }
    }

    if (_trial.entitySet().size(0) == 0)
      return;
    const auto& it = _trial.entitySet().template begin<0>();

    bind(*it, ltest_in, ltrial_in);

    std::vector<LocalCoefficients> lcoeff_in, lcoeff_out;
    for (std::size_t step = 0; step != steps; ++step) {
      lcoeff_in.emplace_back(LocalCoefficients{ ltrial_in, std::cref(coefficients[step]) });
      lcoeff_out.emplace_back(LocalCoefficients{ ltrial_in, std::cref(coefficients[step]) });
    }

    std::vector<LocalResidual> lres_in, lres_out;
    for (std::size_t stage = 0; stage != stages; ++stage) {
      lres_in.emplace_back(LocalResidual{ ltest_in, std::ref(residuals[stage]) });
      lres_out.emplace_back(LocalResidual{ ltest_in, std::ref(residuals[stage]) });
    }

    unbind(ltest_in, ltrial_in);

    const auto& es = _trial.entitySet();

    auto sub_time_step = [&](std::size_t step) {
      return time_point + duration * icoeff.timeWeight(step);
    };

    auto mass_weight = [&](std::size_t stage, std::size_t step) {
      return icoeff.massWeight(stage, step) * InstationaryTraits<dt_position>::massFactor(duration);
    };

    auto stiff_weight = [&](std::size_t stage, std::size_t step) {
      return icoeff.stiffnessWeight(stage, step) * InstationaryTraits<dt_position>::stiffnessFactor(duration);
    };

    // PDELab::forEachElement(es, [=, this, mlop = _mass_lop, slop = _stiff_lop, &es, &_tableau](const auto& entity) mutable {
    auto mlop = _mass_lop;
    auto slop = _stiff_lop;
    for (auto&& entity : elements(es)) {

      if (LocalAssembly::skipEntity(mlop, entity)) {
        if (not LocalAssembly::skipEntity(slop, entity))
          DUNE_THROW(InvalidStateException, "skip methods should yiled the same result");
        return;
      }

      bind(entity, ltest_in, ltrial_in);

      if (LocalAssembly::doVolume(slop) | LocalAssembly::doVolume(mlop)) {
        for (std::size_t stage = 0; stage != stages; ++stage) {
          lres_in[stage].clear(ltest_in);

          for (std::size_t step = 0; step != steps; ++step) {
            TimePoint tp = sub_time_step(step);

            if (do_stiff[stage][step] | do_mass[stage][step])
              lcoeff_in[step].load(ltrial_in);

            if (do_mass[stage][step]) {
              LocalMassResidual lmass_res_in_rhs{lres_in[stage], mass_weight(stage, step)};
              LocalAssembly::volume(mlop, tp, ltrial_in, lcoeff_in[step], ltest_in, lmass_res_in_rhs);
            }

            if (do_stiff[stage][step]) {
              LocalStiffnessResidual lstiff_res_in_rhs{lres_in[stage], stiff_weight(stage, step)};
              LocalAssembly::volume(slop, tp, ltrial_in, lcoeff_in[step], ltest_in, lstiff_res_in_rhs);
            }
          }
        }
      }

      if (  LocalAssembly::doSkeleton(slop) | LocalAssembly::doBoundary(slop)
          | LocalAssembly::doSkeleton(mlop) | LocalAssembly::doBoundary(mlop)) {

        auto id_in = _mapper.index(entity);

        for (const auto& is : intersections(es, entity)) {
          if (LocalAssembly::skipIntersection(mlop, is)) {
            if (not LocalAssembly::skipIntersection(slop, is))
              DUNE_THROW(InvalidStateException, "skip methods should yiled the same result");
            continue;
          }
          if (is.neighbor()) { // interior and periodic cases

            if (LocalAssembly::doSkeleton(slop) or LocalAssembly::doSkeleton(mlop))
              continue;

            const auto& entity_out = is.outside();
            auto id_out = _mapper.index(entity_out);

            // visit face only the first time we see it
            bool skip = id_in < id_out;

            // we have to be sure that outside entity was not skipped,
            // otherwise, neither side will be visited
            if (skip and not LocalAssembly::skipEntity(mlop, entity_out))
              continue;

            bind(entity_out, ltest_out, ltrial_out);

            for (std::size_t stage = 0; stage != stages; ++stage) {
              lres_out[stage].clear(ltest_out);

              for (std::size_t step = 0; step != steps; ++step) {

                TimePoint tp = sub_time_step(step);
                if (do_mass[stage][step] | do_stiff[stage][step])
                  lcoeff_out[step].load(ltrial_out);

                if (do_mass[stage][step]) {
                  LocalMassResidual lmass_res_in_rhs{lres_in[stage],   mass_weight(stage, step)};
                  LocalMassResidual lmass_res_out_rhs{lres_out[stage], mass_weight(stage, step)};
                  LocalAssembly::skeleton(mlop, is, tp, ltrial_in, lcoeff_in[step], ltest_in, ltrial_out, lcoeff_out[step], ltest_out, lmass_res_in_rhs, lmass_res_out_rhs);
                }

                if (do_stiff[stage][step]) {
                  LocalStiffnessResidual lstiff_res_in_rhs{lres_in[stage],   stiff_weight(stage, step)};
                  LocalStiffnessResidual lstiff_res_out_rhs{lres_out[stage], stiff_weight(stage, step)};
                  LocalAssembly::skeleton(slop, is, tp, ltrial_in, lcoeff_in[step], ltest_in, ltrial_out, lcoeff_out[step], ltest_out, lstiff_res_in_rhs, lstiff_res_out_rhs);
                }
              }

              lres_out[stage].fetch_add(ltest_out, std::true_type());
            }
            unbind(ltest_out, ltrial_out);
          } else if (is.boundary() & (LocalAssembly::doBoundary(slop) or LocalAssembly::doBoundary(mlop))) {

            for (std::size_t stage = 0; stage != stages; ++stage) {
              for (std::size_t step = 0; step != steps; ++step) {
                TimePoint tp = sub_time_step(step);
                if (do_mass[stage][step]) {
                  LocalMassResidual lmass_res_in_rhs{lres_in[stage], mass_weight(stage, step)};
                  LocalAssembly::boundary(mlop, is, tp, ltrial_in, lcoeff_in[step], ltest_in, lmass_res_in_rhs);
                }

                if (do_stiff[stage][step]) {
                  LocalStiffnessResidual lstiff_res_in_rhs{lres_in[stage], stiff_weight(stage, step)};
                  LocalAssembly::boundary(slop, is, tp, ltrial_in, lcoeff_in[step], ltest_in, lstiff_res_in_rhs);
                }
              }
            }
          }
        }
      }

      for (std::size_t stage = 0; stage != stages; ++stage)
        lres_in[stage].fetch_add(ltest_in, std::true_type());

      unbind(ltest_in, ltrial_in);
    }
    // });
  }

  const InstationaryCoefficients& getInstationaryCoefficients() const {
    return this->template get<InstationaryCoefficients>("instationary_coefficients");
  }

  const TimePoint& getTimePoint() const {
    return this->template get<TimePoint>("time_point");
  }

  const Duration& getDuration() const {
    return this->template get<Duration>("duration");
  }

protected:
  TrialBasis _trial;
  TestBasis _test;
  MassLocalOperator _mass_lop;
  StiffnessLocalOperator _stiff_lop;

private:
  Mapper _mapper;
};


template<class          Coefficients,
         class          Residual,
         Concept::Basis TrialBasis,
         Concept::Basis TestBasis,
         class          MassLocalOperator,
         class          StiffnessLocalOperator,
         class TimePoint = double,
         class Duration  = double,
         DurationPosition dt_position = DurationPosition::StiffnessNumerator>
auto makeInstationaryForwardAssembler(
                        const TrialBasis& trial,
                        const TestBasis& test,
                        const MassLocalOperator& mass_lop,
                        const StiffnessLocalOperator& stiff_lop)
{
  using Type = InstationaryForwardAssembler<Coefficients,Residual, TrialBasis,TestBasis, MassLocalOperator, StiffnessLocalOperator, TimePoint, Duration, dt_position>;
  return std::make_unique<Type>(trial, test, mass_lop, stiff_lop);
}

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_OPERATOR_FORWARD_INSTATIONARY_ASSEMBLER_FORWARD_HH
