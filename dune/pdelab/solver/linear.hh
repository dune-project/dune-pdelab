#ifndef DUNE_PDELAB_SOLVER_LINEAR_HH
#define DUNE_PDELAB_SOLVER_LINEAR_HH

#include <chrono>
#include <functional>

#include <dune/logging/logger.hh>
#include <dune/pdelab/backend/common/interface.hh>
#include <dune/pdelab/solver/residual.hh>

namespace Dune::PDELab {

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

} // end namespace Dune::PDELab

#endif // DUNE_PDELAB_SOLVER_LINEAR_HH
