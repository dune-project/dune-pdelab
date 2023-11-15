#ifndef DUNE_PDELAB_OPERATOR_INVERSE_NEWTON_HH
#define DUNE_PDELAB_OPERATOR_INVERSE_NEWTON_HH

#include <dune/pdelab/operator/operator.hh>
#include <dune/pdelab/operator/line_search.hh>

#include <dune/pdelab/common/convergence/condition.hh>
#include <dune/pdelab/common/container_entry.hh>
#include <dune/pdelab/common/execution.hh>

#include <memory>
#include <vector>

namespace Dune::PDELab::inline Experimental {

/**
 * @brief Approximation of the Inverse @f$f^{-1}@f$ of a Differentiable operator @f$f@f$ using Newton's method
 *
 * @tparam Domain
 * @tparam Range
 * @tparam RangeField
 */
template<class Domain, class Range, class RangeField = double>
class NewtonOperator : public Inverse<Domain, Range>
{

  using Norm = std::function<RangeField(const Range&)>;
  using Forward = Dune::PDELab::Operator<Domain, Range>;
  using Solver = Dune::PDELab::Inverse<Domain, Range>;
  using ConvergenceCondition = Dune::PDELab::Convergence::Condition<RangeField>;
  using LineSearch = Dune::PDELab::LineSearch<Domain,Range,RangeField>;

public:

  NewtonOperator() {
    PropertyTree& properties = *this;
    properties["verbosity"] = 1;
    properties["convergence_condition"].documentation =
      "ConvergenceCondition object that controls the termination condition of the newton iteration. "
      "If not present, the parameters 'absolute_tolerance' and 'relative_tolerance' will be added.";
    properties["convergence_condition"].setter = [](Property& ppt) {
      if (not ppt.has_value()) return;
      auto& ppt_tree = ppt.as_tree();
      if (not ppt_tree.hasKey("relative_tolerance"))
        ppt_tree["relative_tolerance"] = Convergence::makeRelativeTolerance();
      if (not ppt_tree.hasKey("absolute_tolerance"))
        ppt_tree["absolute_tolerance"] = Convergence::makeAbsoluteTolerance<RangeField>();
    };
    properties["convergence_condition"] = std::shared_ptr<ConvergenceCondition>(std::make_shared<Convergence::DefaultCondition<RangeField>>());

    properties["line_search"].documentation =
      "LineSearch object that controls the line search stage at each newton iteration";
    properties["line_search"] = std::shared_ptr<LineSearch>(std::make_shared<LineSearchNone<Domain,Range,RangeField>>());

    properties["linearization_threshold"].documentation = "Threshold to linearize the derivative ..";
    properties["linearization_threshold"].setter = [](const Property& value){
      if (Dune::FloatCmp::lt(unwrap_property_ref<const double>(value), 0.))
        DUNE_THROW(RangeError, "Linearization threshold has to be bigger or equal than zero");
    };
    properties["linearization_threshold"] = 0.0;

    properties["dx_inverse"].documentation
      = "An operator able solve (approximately) the inverse of the derivative of x. "
      "During the course of newton iterations, dx itself will be set in dx_inverse[\"forward\"]";
    properties["dx_inverse"].setter = [=](Property& dx_inv_pptr){
      if (dx_inv_pptr.has_value())
        [[maybe_unused]] Solver& dx_inv = unwrap_property_ref<Solver>(dx_inv_pptr);
    };

    if (not properties.hasKey("dx_inverse_fixed_tolerance"))
      properties["dx_inverse_fixed_tolerance"] = false;
    properties["dx_inverse_fixed_tolerance"].documentation
      = "Whether dx_inverse relative reduction is fixed between newton iterations. "
        "The newton method may have a better guess for a relative reduction on the for "
        "the dx_inverse correction at each iteration. If this is set to true, the newton "
        "application will modify the `dx_inverse[\"convergence_condition.relative_tolerance\"]` "
        "property before each iteration of the newton method";

    if (not properties.hasKey("dx_inverse_min_relative_tolerance")) {
      properties["dx_inverse_min_relative_tolerance"].setter = [](const Property& min_rl_rd) {
        double value = unwrap_property_ref<const double>(min_rl_rd);
        if (Dune::FloatCmp::lt(value, 0.))
          DUNE_THROW(RangeError, "Derivative inverse minimum relative reduction has to be bigger or equal than zero");
        if (Dune::FloatCmp::gt(value, 1.))
          DUNE_THROW(RangeError, "Derivative inverse minimum relative reduction has to be less or equal than one");
      };
      properties["dx_inverse_min_relative_tolerance"] = 0.1;
    }

    properties["dx_inverse_min_relative_tolerance"].documentation
      = "Minimum relative reduction of the dx_inverse. That is, if "
        "`dx_inverse_fixed_tolerance == true`, then property "
        "`dx_inverse[\"convergence_condition.relative_tolerance\"]` "
        "will always be at least the value set in the parameter";
  }

  NewtonOperator(const NewtonOperator&) = default;

  //! Solve the nonlinear problem using x as initial guess and for storing the result
  virtual ErrorCondition apply(const Range& init_residual, Domain& x) override
  {
    [[maybe_unused]] auto newton_timestamp = perfetto::TrackEvent::GetTraceTimeNs();
    TRACE_EVENT("dune", "Newton", newton_timestamp);
    assert(_norm_op);

    auto& conv_cond = getConvergenceCondition();
    auto& dx_inverse = getDerivativeInverse();
    auto& line_search = getLineSearch();
    auto& forward = getForward();
    const int verbosity = this->template get<int>("verbosity");

    // Calculate initial defect
    Range residual = init_residual;
    forward.apply(std::as_const(x), residual).or_throw();

    RangeField defect = _norm_op(residual);
    std::vector<RangeField> defects{defect};
    TRACE_COUNTER("dune", "Newton::Defect", defect);
    if (verbosity >= 2)
      std::cout << "  Initial defect: "
                << std::setw(12) << std::setprecision(4) << std::scientific
                << defect << std::endl;

    Domain correction = x;
    forEachContainerEntry(PDELab::Execution::par_unseq, correction, []<class T>(T& v){v = T{0};});
    Domain zero = correction;


    ErrorCondition error_condition;
    Convergence::Reason convergence_reason;
    while (Convergence::Reason::Iterating == (convergence_reason = conv_cond.evaluate(defects))){
      [[maybe_unused]] auto it_timestamp = perfetto::TrackEvent::GetTraceTimeNs();
      TRACE_EVENT("dune", "Newton::Iteration", it_timestamp);

      prepareStep(defects, x);
      TRACE_COUNTER("dune",
                    "Newton::InverseTargetRelativeTolerance",
                    it_timestamp,
                    dx_inverse.template get<double>(
                      "convergence_condition.relative_tolerance"));

      if (defects.size() > 1) correction = zero;

      error_condition = dx_inverse.apply(residual, correction);
      if (error_condition) break;

      defect = defects.back();
      residual = init_residual;
      error_condition = line_search.apply(forward, _norm_op, x, correction, residual, defect);
      if (error_condition) break;

      TRACE_COUNTER("dune", "Newton::Defect", it_timestamp, defect);
      if (verbosity == 2)
        std::cout << "  Newton iteration "
                  << std::setw(2)
                  << defects.size()
                  << ".  New defect: "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << defect
                  << ".  Reduction (this): "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << defects.back()/defect
                  << ".  Reduction (total): "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << defects.front()/defect
                  << std::endl;
      defects.push_back(defect);
    }
    dx_inverse["forward"] = nullptr;
    [[maybe_unused]] auto iterations = defects.size() - 1;

    if (not error_condition)
      error_condition = make_error_condition(convergence_reason);
    [[maybe_unused]] auto reduction = defects.back()/defects.front();
    TRACE_COUNTER("dune", "Newton::Iterations",        newton_timestamp, iterations);
    TRACE_COUNTER("dune", "Newton::Reduction",         newton_timestamp, reduction);
    TRACE_COUNTER("dune", "Newton::Converged",         newton_timestamp, bool(not error_condition));
    TRACE_COUNTER("dune", "Newton::ConvergenceRate",   newton_timestamp, std::pow(reduction, 1.0/iterations));
    TRACE_COUNTER("dune", "Newton::ConvergenceReason", newton_timestamp, static_cast<int>(convergence_reason));

    return error_condition;
  }

  void setNorm(const std::function<RangeField(const Range&)>& norm_op) {
    _norm_op = norm_op;
  }

private:

  void prepareStep(std::span<const RangeField> defects, const Domain& x)
  {
    const PropertyTree& properties = *this;
    auto rel_tol                = properties.get<double     >("convergence_condition.relative_tolerance");
    auto dx_inverse_fixed_tol   = properties.get<bool       >("dx_inverse_fixed_tolerance");
    auto min_dx_inverse_rel_tol = properties.get<double     >("dx_inverse_min_relative_tolerance");
    auto lin_threshold          = properties.get<double     >("linearization_threshold");

    // absolute tolerance is optional
    std::optional<RangeField> abs_tol;
    const auto& abs_tol_ppt = properties.get("convergence_condition.absolute_tolerance");
    if (abs_tol_ppt.has_value()) abs_tol = unwrap_property_ref<const RangeField>(abs_tol_ppt);

    assert(not empty(defects));
    auto first_defect = defects.front();
    auto current_defect = defects.back();
    auto previous_defect = *std::next(rbegin(defects), size(defects) > 1);
    auto defect_rate = current_defect/previous_defect;
    std::size_t iteration = defects.size() - 1;

    // if defect is above the linearization threshold, linearize differential
    // operator and reset linear operator for the dx_inverse
    bool linearize = (iteration == 0) or (defect_rate > lin_threshold);
    if (linearize)
      _derivative = getForward().derivative(
        x, _derivative ? std::move(_derivative) : nullptr);
    auto& dx_inverse = getDerivativeInverse();
    dx_inverse["forward"] = std::shared_ptr<Forward>(_derivative);

    using std::min, std::max, std::clamp;

    // case where we are allowed to set a relative tolerance
    if (not dx_inverse_fixed_tol) {

      double it_dx_inverse_rel_tol = min_dx_inverse_rel_tol;
      double max_dx_inverse_rel_tol = std::numeric_limits<double>::infinity();

      // in case of the first iteration we stick to min_dx_inverse_rel_tol to force maximum reduction
      if (iteration !=0) {
        // to achieve second order convergence of newton we need a linear
        // reduction of at least defect_rate^2.
        it_dx_inverse_rel_tol = defect_rate*defect_rate;

        // Determine maximum defect where Newton is converged
        RangeField stop_defect = max(first_defect * rel_tol, abs_tol.value_or(0.));

        // since we know the target stop defect, we can set an upper bound to the linear reduction
        // this will set a more relaxed reduction requirement on the last newton iteration
        max_dx_inverse_rel_tol = stop_defect/(10*current_defect);
      }
      it_dx_inverse_rel_tol = clamp<double>(it_dx_inverse_rel_tol, min_dx_inverse_rel_tol, max_dx_inverse_rel_tol);

      dx_inverse["convergence_condition.relative_tolerance"] = it_dx_inverse_rel_tol;
    }
  }


  Forward& getForward() {
    return this->template get<Forward>("forward");
  }

  Solver& getDerivativeInverse() {
    return this->template get<Solver>("dx_inverse");
  }

  LineSearch& getLineSearch() {
    return this->template get<LineSearch>("line_search");
  }

  ConvergenceCondition& getConvergenceCondition() {
    return this->template get<ConvergenceCondition>("convergence_condition");
  }

private:
  std::function<RangeField(const Range&)> _norm_op;
  std::shared_ptr<Operator<Domain, Range>> _derivative;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_OPERATOR_INVERSE_NEWTON_HH
