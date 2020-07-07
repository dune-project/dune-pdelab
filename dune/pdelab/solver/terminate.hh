// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_TERMINATE_HH
#define DUNE_PDELAB_SOLVER_TERMINATE_HH

namespace Dune::PDELab
{

  class TerminateError : public Exception {};


  class TerminateInterface
  {
  public:
    //! Every abstract base class should have a virtual destructor
    virtual ~TerminateInterface () {}

    virtual bool terminate() = 0;

    virtual void setParameters(const ParameterTree&) = 0;
  };


  template <typename Solver>
  class DefaultTerminate : public TerminateInterface
  {
  public:
    using Real = typename Solver::Real;

    DefaultTerminate(Solver& solver) : _solver(solver) {}

    virtual bool terminate() override
    {
      if (_force_iteration && _solver.result().iterations == 0)
        return false;
      auto converged = _solver.result().defect < _solver.getAbsoluteLimit() || _solver.result().defect < _solver.result().first_defect * _solver.getReduction();
      if (_solver.result().iterations >= _maxIterations && not _solver.result().converged)
        DUNE_THROW(TerminateError,
                   "Terminate::terminate(): Maximum iteration count reached");
      return converged;
    }

    virtual void setParameters(const ParameterTree& parameterTree) override
    {
      _maxIterations = parameterTree.get<unsigned int>("MaxIterations", _maxIterations);
      _force_iteration = parameterTree.get<bool>("ForceIteration", _force_iteration);

      std::cout << _maxIterations << std::endl;
    }

  private:
    Solver& _solver;
    unsigned int _maxIterations = 40;
    bool _force_iteration = false;
  };
}

#endif
