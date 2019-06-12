#ifndef DUNE_PDELAB_SOLVER_RESIDUAL_HH
#define DUNE_PDELAB_SOLVER_RESIDUAL_HH

#include <memory>

#include <dune/pdelab/backend/common/interface.hh>

namespace Dune::PDELab {

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

} // end namespace Dune::PDELab

#endif // DUNE_PDELAB_SOLVER_RESIDUAL_HH
