#ifndef DUNE_PDELAB_SOLVER_NEWTON_HH
#define DUNE_PDELAB_SOLVER_NEWTON_HH

#include <chrono>
#include <functional>

#include <dune/logging/logger.hh>
#include <dune/pdelab/backend/common/interface.hh>
#include <dune/pdelab/solver/residual.hh>
#include <dune/pdelab/solver/utility.hh>

namespace Dune::PDELab {


  constexpr inline auto defaultNewtonTermination = [](auto& newton, const ParameterTree& params)
  {
    using Real = typename std::decay_t<decltype(newton)>::Real;

    std::size_t max_iterations = params.get<std::size_t>("iterations",40);
    bool force_iteration = params.get<bool>("force_iteration",false);
    Real reduction = params.get<Real>("reduction",1e-8);
    Real absolute_limit = params.get<Real>("absolute_limit",1e-12);

    newton.logger().info(2,
      "defaultTermination config: max_iterations={} force_iteration={} reduction={} abs_limit={}"_fmt,
      max_iterations,force_iteration,reduction,absolute_limit);

    return [=,&newton](auto& solution)
    {
      auto iterations = newton.iteration();

      using std::max;
      auto defect = newton.defect();
      auto initial_defect = newton.initialDefect();
      Real stop_defect = max(absolute_limit,initial_defect * reduction);

      bool converged = newton.defect() < stop_defect;

      if (force_iteration and iterations == 0)
        converged = false;

      if (not converged and iterations > max_iterations)
        DUNE_THROW(Exception,"not converged");

      return std::make_tuple(converged,stop_defect);
    };
  };


  constexpr inline auto defaultNewtonLineSearch = [](auto& newton, const ParameterTree& params)
  {
    using Real = typename std::decay_t<decltype(newton)>::Real;

    auto strategy = LineSearchStrategy::hackbuschReusken;
    if (params.hasKey("strategy"))
      strategy = lineSearchStrategyFromString(params["strategy"]);

    std::size_t max_iterations = params.get<std::size_t>("iterations",10);
    Real damping_factor = params.get<Real>("damping_factor",0.5);
    Real initial_factor   = params.get<Real>("initial_factor",1.0);

    newton.logger().info(2,"defaultLinesearch config: strategy={} max_iterations={} damping_factor={} initial_factor={}"_fmt,
                         strategy,max_iterations,damping_factor,initial_factor);

    return [=,&newton](auto& solution, auto& correction)
    {
      auto log = newton.logger();
      log.indent(4);

      if (strategy == LineSearchStrategy::none)
      {
        log.detail("line search disabled"_fmt);
        solution -= correction;
        newton.updateDefect(solution);
        return true;
      }

      Real lambda = initial_factor;
      Real best_lambda = 0.0;
      Real best_defect = newton.defect();

      // create a temporary copy of the solution vector
      std::decay_t<decltype(solution)> original(solution);

      std::size_t i = 0;

      bool converged = false;

      for (; i < max_iterations ; ++i)
      {

        log.debug("trying line search damping factor: {:#16.4}"_fmt,lambda);

        solution.axpy(-lambda,correction);
        Real defect = newton.updateDefect(solution);

        if (defect <= (1.0 - lambda/4) * newton.previousDefect())
        {
          log.debug("line search converged"_fmt);
          converged = true;
          break;
        }

        if (defect < best_defect)
        {
          best_defect = defect;
          best_lambda = lambda;
        }

        lambda *= damping_factor;
        solution = original;
      }

      if (not converged)
      {
        log.debug("max line search iterations ({}) exceeded"_fmt,max_iterations);
        switch (strategy)
        {
        case LineSearchStrategy::hackbuschReusken:
          solution = original;
          newton.updateDefect(solution);
          break;

        case LineSearchStrategy::hackbuschReuskenAcceptBest:
          if (best_lambda == 0.0)
          {
            solution = original;
            newton.updateDefect(solution);
          }
          if (lambda != best_lambda)
          {
            solution = original;
            solution.axpy(-best_lambda,correction);
            newton.updateDefect(solution);
          }
          log.debug("best damping factor found: {:#16.4}"_fmt,best_lambda);
          converged = true;

        default:
          // should never get here!
          std::abort();
        }
      }

      if (converged)
        log.detail("line search damping factor: {:#16.4}"_fmt,lambda);

      return converged;

    };
  };


  template<typename Domain_, typename Range_>
  class NewtonPDESolver
    : public PDESolver<Domain_, Range_>
  {

  public:

    using Domain = Domain_;
    using Range  = Range_;
    using Real   = typename PDESolver<Domain_,Range_>::Real;

    using LinearSolver      = Backend::LinearSolver<Domain,Range>;
    using ResidualEvaluator = Dune::PDELab::ResidualEvaluator<Domain,Range>;
    using Termination       = std::function<std::tuple<bool,Real>(Domain&)>;
    using LineSearch        = std::function<bool(Domain&,Domain&)>;

    void solve(Domain& solution) override
    {
      _iteration = 0;

      using Clock = std::chrono::steady_clock;
      using Duration = Clock::duration;

      auto assembler_time = Duration::zero();
      auto linear_solver_time = Duration::zero();
      auto line_search_time   = Duration::zero();

      auto to_seconds = [](Duration duration)
      {
        return std::chrono::duration<double>(duration).count();
      };

      auto start_solve = Clock::now();

      updateDefect(solution);
      _initial_defect = _prev_defect = _defect;
      _log.notice("Initial defect: {:12.4e}"_fmt,_initial_defect);

      for (auto [converged,stop_defect] = _terminate(solution) ;
           not converged ;
           std::tie(converged,stop_defect) = _terminate(solution)
        )
      {
        using std::abs;
        using std::max;
        using std::min;
        using std::pow;

        _log.info("Newton iteration {:3}"_fmt,_iteration);
        _log.info("{:-<30}"_fmt,"");

        auto start = Clock::now();
        _solver->setLinearizationPoint(solution, abs(_defect / _prev_defect) > _reassemble_threshold);
        auto end = Clock::now();
        assembler_time += time;

        _log.info(2,"assembly time: {:30.4e}"_fmt,to_seconds(end-start));

        Real linear_reduction = _min_linear_reduction;
        if (not _fixed_linear_reduction)
        {

          /*
           * To achieve second order convergence of newton
           * we need a linear reduction of at least
           * current_defect^2/prev_defect^2.
           * For the last newton step a linear reduction of
           * 1/10*end_defect/current_defect
           * is sufficient for convergence.
           */

          if (stop_defect / (10 * _defect) > _defect * _defect / _prev_defect * _prev_defect)
            linear_reduction = stop_defect / (10 * _defect);
          else
            linear_reduction = min(_min_linear_reduction,_defect * _defect / (_prev_defect * _prev_defect));

        }

        _prev_defect = _defect;

        _correction = 0.0;

        _log.info(2,"Solving linear system..."_fmt);
        start = Clock::now();
        _solver->solve(_correction,_residual,linear_reduction);
        end = Clock::now();
        linear_solver_time += end - start;

        _log.detail(2,"linear solver iterations: {:16}"_fmt,_iteration);
        _log.detail(2,"linear defect reduction:  {:16.4e}"_fmt,_reduction);

        _log.info(2,"Performing line search..."_fmt);
        start = Clock::now();
        bool linesearch_converged = _line_search(solution,_correction);
        end = Clock::now();
        line_search_time += end - start;

        if (not linesearch_converged)
          _log.warning("Line search did not converge"_fmt);

        _reduction = _defect / _initial_defect;
        ++_iteration;
        _convergence_rate = pow(_reduction, 1.0/_iteration);

      }


      (*_eval_residual)(solution,_residual);
      _solver->solve(_correction,_residual);
      solution -= _correction;
    }

    Logging::Logger logger() const
    {
      return _log;
    }

    Real defect() const
    {
      return _defect;
    }

    Real initialDefect() const
    {
      return _initial_defect;
    }

    Real previousDefect() const
    {
      return _prev_defect;
    }

    Real updateDefect(const Domain& solution)
    {
      _residual = 0.0;
      (*_eval_residual)(solution,_residual);
      return _defect = _solver->norm(_residual);
    }

    std::size_t iteration() const
    {
      return _iteration;
    }

    int startStep(Real time, Real dt) override
    {
      _solver->startStep(time,dt);
      return _eval_residual->startStep(time,dt);
    }

    bool acceptStage(int stage, const Domain& solution) override
    {
      bool changed = _solver->acceptStage(stage,solution);
      changed |= _eval_residual->acceptStage(stage,solution);
      return changed;
    }

    NewtonPDESolver(
      std::shared_ptr<LinearSolver> solver,
      std::shared_ptr<ResidualEvaluator> eval_residual,
      const ParameterTree& params
      )
      : _log(Dune::PDELab::Logging::componentLogger(params,"newton"))
      , _solver(std::move(solver))
      , _eval_residual(std::move(eval_residual))
      , _correction(_eval_residual->makeDomain())
      , _residual(_eval_residual->makeRange())
      , _reassemble_threshold(params.get<Real>("reassemble_threshold",0.0))
      , _min_linear_reduction(params.get<Real>("min_linear_reduction",1e-3))
      , _fixed_linear_reduction(params.get("fixed_linear_reduction",false))
      , _log(Logging::componentLogger(params,"newton"))
    {
      _log.info("Newton config: reassemble_threshold={} min_linear_reduction={} fixed_linear_reduction={}"_fmt,
                _reassemble_threshold,_min_linear_reduction,_fixed_linear_reduction);
      // we have to initialize these two in the body, otherwise we get undefined behavior because we
      // are passing in a reference to ourselves
      _terminate = defaultNewtonTermination(*this,params.sub("terminate"));
      _line_search = defaultNewtonLineSearch(*this,params.sub("line_search"));
    }

  private:

    Logger _log;

    std::shared_ptr<LinearSolver> _solver;
    std::shared_ptr<ResidualEvaluator> _eval_residual;
    Termination _terminate;
    LineSearch _line_search;
    Domain _correction;
    Range _residual;

    Real _defect = 0.0;
    Real _initial_defect = 0.0;
    Real _prev_defect = 0.0;
    std::size_t _iteration = 0;
    Real _convergence_rate = 0.0;
    Real _reduction = 0.0;

    Real _reassemble_threshold;
    bool _fixed_linear_reduction;
    Real _min_linear_reduction;

    Logging::Logger _log;

  };

} // end namespace Dune::PDELab

#endif // DUNE_PDELAB_SOLVER_NEWTON_HH
