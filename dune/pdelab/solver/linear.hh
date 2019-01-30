#ifndef DUNE_PDELAB_SOLVER_LINEAR_HH
#define DUNE_PDELAB_SOLVER_LINEAR_HH

#include <chrono>

#include <dune/logging/logger.hh>
#include <dune/pdelab/backend/common/interface.hh>

namespace Dune::PDELab {

  template<typename Domain_, typename Range_>
  struct PDESolver
    : public Backend::TimestepAwareComponent<Domain_>
  {

    using Domain = Domain_;
    using Range  = Range_;

    virtual void solve(Domain& domain) = 0;

  };

  template<typename Domain_, typename Range_>
  struct ResidualEvaluator
    : public Backend::TimestepAwareComponent<Domain_>
  {

    using Domain = Domain_;
    using Range  = Range_;

    virtual void operator()(const Domain& domain, Range& range) = 0;

    virtual Range makeRange() = 0;
    virtual Domain makeDomain() = 0;

  };

  template<typename GO>
  class GridOperatorBasedResidual
    : public ResidualEvaluator<typename GO::Traits::Domain,typename GO::Traits::Range>
  {

    using Base = ResidualEvaluator<typename GO::Traits::Domain,typename GO::Traits::Range>;

  public:

    using GridOperator = GO;
    using Domain = typename GridOperator::Traits::Domain;
    using Range  = typename GridOperator::Traits::Range;
    using Real   = typename Base::Real;

    void operator()(const Domain& domain, Range& range) override
    {
      _grid_operator->residual(domain,range);
    }

    int startStep(Real time, Real dt) override
    {
      return _grid_operator->startStep(time,dt);
    }

    bool acceptStage(int stage, const Domain& solution) override
    {
      return _grid_operator->acceptStage(stage,solution);
    }

    Range makeRange() override
    {
      return Range(_grid_operator->testSpace());
    }

    Domain makeDomain() override
    {
      return Domain(_grid_operator->trialSpace());
    }

    GridOperatorBasedResidual(std::shared_ptr<GridOperator> grid_operator)
      : _grid_operator(std::move(grid_operator))
    {}

  private:

    std::shared_ptr<GridOperator> _grid_operator;

  };


  template<typename Domain_, typename Range_>
  class LinearPDESolver
    : public PDESolver<Domain_, Range_>
  {

  public:

    using Domain = Domain_;
    using Range  = Range_;

    using LinearSolver      = Backend::LinearSolver<Domain,Range>;
    using ResidualEvaluator = Dune::PDELab::ResidualEvaluator<Domain,Range>;

    using Real = typename LinearSolver::Real;

    void solve(Domain& solution) override
    {
      (*_eval_residual)(solution,_residual);
      _solver->solve(_correction,_residual);
      solution -= _correction;
    }

    LinearPDESolver(
      std::shared_ptr<LinearSolver> solver,
      std::shared_ptr<ResidualEvaluator> eval_residual,
      const ParameterTree& params
      )
      : _solver(std::move(solver))
      , _eval_residual(std::move(eval_residual))
      , _correction(_eval_residual->makeDomain())
      , _residual(_eval_residual->makeRange())
    {}

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

  private:

    std::shared_ptr<LinearSolver> _solver;
    std::shared_ptr<ResidualEvaluator> _eval_residual;
    Domain _correction;
    Range _residual;

  };

  auto defaultNewtonTermination = [](auto& newton, const ParameterTree& params)
  {
    using Real = typename std::decay_t<decltype(newton)>::Real;

    std::size_t max_iterations = params.get<std::size_t>("iterations",40);
    bool force_iteration = params.get<bool>("force_iteration",false);
    Real reduction = params.get<Real>("reduction",1e-8);
    Real absolute_limit = params.get<Real>("absolute_limit",1e-12);

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

  enum class LineSearchStrategy
  {
    none,
    hackbuschReusken,
    hackbuschReuskenAcceptBest
  };

  LineSearchStrategy lineSearchStrategyFromString(const std::string& name)
  {
    if (name == "none")
      return LineSearchStrategy::none;
    if (name == "hackbusch_reusken")
      return LineSearchStrategy::hackbuschReusken;
    if (name == "hackbusch_reusken_accept_best")
      return LineSearchStrategy::hackbuschReuskenAcceptBest;
    DUNE_THROW(Exception,"Unkown line search strategy: " << name);
  }

  auto defaultNewtonLineSearch = [](auto& newton, const ParameterTree& params)
  {
    using Real = typename std::decay_t<decltype(newton)>::Real;

    auto strategy = LineSearchStrategy::hackbuschReusken;
    if (params.hasKey("strategy"))
      strategy = lineSearchStrategyFromString(params["strategy"]);

    std::size_t max_iterations = params.get<std::size_t>("iterations",10);
    Real damping_factor = params.get<Real>("damping_factor",0.5);
    Real initial_factor   = params.get<Real>("initial_factor",1.0);

    return [=,&newton](auto& solution, auto& correction)
    {
      auto log = newton.logger();
      log.indent(4);

      if (strategy == LineSearchStrategy::none)
      {
        log.info("line search disabled"_fmt);
        solution -= correction;
        newton.updateDefect(solution);
        return;
      }

      Real lambda = initial_factor;
      Real best_lambda = 0.0;
      Real best_defect = newton.defect();

      // create a temporary copy of the solution vector
      std::decay_t<decltype(solution)> original(solution);

      std::size_t i = 0;

      for (;;)
      {

        log.info("trying line search damping factor: {:#16.4}"_fmt,lambda);

        solution.axpy(-lambda,correction);
        Real defect = newton.updateDefect(solution);

        if (defect <= (1.0 - lambda/4) * newton.previousDefect())
        {
        log.info("line search converged"_fmt);
          break;
        }

        if (defect < best_defect)
        {
          best_defect = defect;
          best_lambda = lambda;
        }

        if (++i >= max_iterations)
        {
          log.info("max line search iterations ({}) exceeded"_fmt,max_iterations);
          switch (strategy)
          {
          case LineSearchStrategy::hackbuschReusken:
            solution = original;
            newton.updateDefect(solution);
            DUNE_THROW(Exception,"line search failed");

          case LineSearchStrategy::hackbuschReuskenAcceptBest:
            if (best_lambda == 0.0)
            {
              solution = original;
              newton.updateDefect(solution);
              DUNE_THROW(Exception,"line search failed");
            }
            if (lambda != best_lambda)
            {
              solution = original;
              solution.axpy(-best_lambda,correction);
              newton.updateDefect(solution);
            }
            log.info("best damping factor found: {:#16.4}"_fmt,best_lambda);
            return;

          default:
            // should never get here!
            std::abort();
          }
        }

        lambda *= damping_factor;
        solution = original;
      }

      log.info("line search damping factor: {:#16.4}"_fmt,lambda);

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
    using LineSearch        = std::function<void(Domain&,Domain&)>;

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
        assembler_time += end - start;

        _log.info(2,"assembly time: {:30.4e}"_fmt,to_seconds(end-start));

        Real linear_reduction = _min_linear_reduction;
        if (not _fixed_linear_reduction)
        {

          /*
            To achieve second order convergence of newton
            we need a linear reduction of at least
            current_defect^2/prev_defect^2.
            For the last newton step a linear reduction of
            1/10*end_defect/current_defect
            is sufficient for convergence.
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
        _line_search(solution,_correction);
        end = Clock::now();
        line_search_time += end - start;

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
      : _solver(std::move(solver))
      , _eval_residual(std::move(eval_residual))
      , _terminate(defaultNewtonTermination(*this,params.sub("terminate")))
      , _line_search(defaultNewtonLineSearch(*this,params.sub("line_search")))
      , _correction(_eval_residual->makeDomain())
      , _residual(_eval_residual->makeRange())
      , _reassemble_threshold(params.get<Real>("reassemble_threshold",0.0))
      , _min_linear_reduction(params.get<Real>("min_linear_reduction",1e-3))
      , _fixed_linear_reduction(params.get("fixed_linear_reduction",false))
      , _log(Logging::logger(params))
    {}

  private:

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

#endif // DUNE_PDELAB_SOLVER_LINEAR_HH
