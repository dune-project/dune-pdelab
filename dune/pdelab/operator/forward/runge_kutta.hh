#ifndef DUNE_PDELAB_OPERATOR_FORWARD_RUNGE_KUTTA_HH
#define DUNE_PDELAB_OPERATOR_FORWARD_RUNGE_KUTTA_HH

#include <dune/pdelab/operator/operator.hh>

#include <dune/pdelab/operator/forward/instationary/coefficients.hh>
#include <dune/pdelab/common/trace.hh>

#include <memory>
#include <vector>

namespace Dune::PDELab::inline Experimental {

/**
 * @brief OneStep operator @f$f_{\delta_t}@f$ stepping a Differentiable Instationary operator
 *
 * @tparam DomainStages
 * @tparam RangeStages
 * @tparam TimePoint
 * @tparam Duration
 */
template<class DomainStages,
         class RangeStages,
         class TimePoint = double,
         class Duration  = double>
requires (std::ranges::range<DomainStages>  && std::copyable<DomainStages> && std::constructible_from<DomainStages, std::size_t> &&
          std::ranges::range<RangeStages>   && std::copyable<RangeStages>  && std::constructible_from<RangeStages,  std::size_t>)
class RungeKutta : public OneStep<std::ranges::range_value_t<DomainStages>> {
  using Inverse = Dune::PDELab::Inverse<DomainStages, RangeStages>;
  using Forward = Dune::PDELab::Operator<DomainStages, RangeStages>;
  using Domain = std::ranges::range_value_t<DomainStages>;
  using RangeStage = std::ranges::range_value_t<RangeStages>;
public:

  ErrorCondition apply(const Domain& start_stage, Domain& final_stage) override
  {
    Domain start_stage_copy = start_stage;
    return this->apply(start_stage_copy, final_stage);
  }

  ErrorCondition apply(Domain& start_stage, Domain& final_stage) override
  {
    auto& forward = getForward();
    auto& tableau = getInstationaryCoefficients();
    if (forward["instationary_coefficients"].has_value())
      DUNE_THROW(RangeError, "No other parameters[\"inverse.forward.instationary_coefficients\"] must be set in the forward model before applying the runge-kutta stepping");

    ErrorCondition ec;
    if (tableau.type() == InstationaryCoefficients::Type::FullyImplicit)
      ec = applyFullyImplicit(start_stage, final_stage);
    else
      ec = applySemiImplicit(start_stage, final_stage);

    if (not ec)
      forward.template get<TimePoint>("time_point") += forward.template get<Duration>("duration");
    forward["instationary_coefficients"] = nullptr;
    return ec;
  }

private:
  ErrorCondition applyFullyImplicit(Domain& start_stage, Domain& final_stage) {
    TRACE_EVENT("dune", "RungeKutta");
    ErrorCondition error_condition;
    auto& tableau = getInstationaryCoefficients();
    auto& inverse = getInverse();
    auto& forward = getForward();

    assert(tableau.extent(0) > 0);

    // set up entire tableau to be solved at once
    forward["instationary_coefficients"] = tableau;

    DomainStages stage_coeffs(tableau.extent(1));
    stage_coeffs[0] = std::move(start_stage);
    try {
      for (std::size_t i = 1; i != stage_coeffs.size(); ++i)
        stage_coeffs[i] = stage_coeffs[0];

      RangeStages step_residuals(tableau.extent(0));
      for (auto& step_residual : step_residuals)
        step_residual = getInitialResidual();

      error_condition = inverse.apply(step_residuals, stage_coeffs);

      if (not error_condition)
        final_stage = stage_coeffs[tableau.extent(1)-1];
    } catch(...) {
      start_stage = std::move(stage_coeffs[0]);
      throw;
    }

    start_stage = std::move(stage_coeffs[0]);
    return error_condition;
  }

  ErrorCondition applySemiImplicit(Domain& start_stage, Domain& final_stage) {
    TRACE_EVENT("dune", "RungeKutta");
    auto& tableau = getInstationaryCoefficients();
    auto& forward = getForward();
    auto& inverse = getInverse();

    assert(tableau.type() != InstationaryCoefficients::Type::FullyImplicit);
    assert(tableau.extent(0) > 0);

    // set up coefficients for this stage with a copy of the intial coefficients
    DomainStages stage_coeff(1);
    stage_coeff[0] = start_stage;

    DomainStages pre_stage_coeff(1);
    pre_stage_coeff[0] = std::move(start_stage);

    ErrorCondition error_condition;

    try {
      RangeStages pre_stage_residual(1);

      // each stage solves one row of the tableau
      for (std::size_t stage = 0; stage != tableau.extent(0); ++stage) {
        TRACE_EVENT("dune", "RungeKutta::Stage");

        assert(pre_stage_coeff.size() == stage+1);
        assert(pre_stage_residual.size() == 1);

        // set up explicit part of tableau row
        forward["instationary_coefficients"] = tableau.slice(stage, 1, 0, stage+1);

        // reset residual the initial value (operators are additive)
        pre_stage_residual[0] = getInitialResidual();
        // compute residual for pre-stages. Since previous stages are already
        // known (i.e. semi-implicit), we don't need to re-compute their residual
        // on every iteration of the inverse operator.
        error_condition = forward.apply(pre_stage_coeff, pre_stage_residual);
        if (error_condition) break;

        // set up coefficients for this stage with a copy of the last used coefficients
        // (only one entry because we are solving the diagonal part of the tableau row)
        stage_coeff[0] = (stage==0) ? final_stage : pre_stage_coeff[stage];

        // set up implicit part (i.e. diagonal) of tableau row
        forward["instationary_coefficients"] = tableau.slice(stage, 1, stage+1, 1);
        // compute operator inverse
        error_condition = inverse.apply(pre_stage_residual, stage_coeff);
        if (error_condition) break;

        // move solved coefficients to the pre_stage coeffs
        if (stage+1 != tableau.extent(0)) {
          DomainStages tmp(stage+2);
          for (std::size_t i = 0; i != stage+1; ++i)
            tmp[i] = std::move(pre_stage_coeff[i]);
          pre_stage_coeff = std::move(tmp);
          pre_stage_coeff[stage+1] = std::move(stage_coeff[0]);
        }
      }

      // set coefficients to the return objects
      start_stage = std::move(pre_stage_coeff[0]);
      final_stage = std::move(stage_coeff[0]);

    } catch (...) {
      start_stage = std::move(pre_stage_coeff[0]);
      throw;
    }

    return error_condition;
  }

  Inverse& getInverse() {
    return this->template get<Inverse>("inverse");
  }

  Forward& getForward() {
    return getInverse().template get<Forward>("forward");
  }

  const InstationaryCoefficients& getInstationaryCoefficients() const {
    return this->template get<InstationaryCoefficients>("instationary_coefficients");
  }

  const RangeStage& getInitialResidual() const {
    return this->template get<const RangeStage>("initial_residual");
  }
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_OPERATOR_FORWARD_RUNGE_KUTTA_HH
