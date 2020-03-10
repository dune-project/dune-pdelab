// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_NEWTONTERMINATE_HH
#define DUNE_PDELAB_SOLVER_NEWTONTERMINATE_HH

namespace Dune::PDELab
{
  class TerminateInterface
  {
  public:
    //! Every abstract base class should have a virtual destructor
    virtual ~TerminateInterface () {}

    virtual bool terminate() = 0;

    virtual void setParameters(const ParameterTree&) = 0;
  };


  template <typename Newton>
  class DefaultTerminate : public TerminateInterface
  {
  public:
    using Real = typename Newton::Real;

    DefaultTerminate(Newton& newton) : _newton(newton) {}

    virtual bool terminate() override
    {
      if (_force_iteration && _newton.result().iterations == 0)
        return false;
      auto converged = _newton.result().defect < _newton.getAbsoluteLimit() || _newton.result().defect < _newton.result().first_defect * _newton.getReduction();
      if (_newton.result().iterations >= _maxIterations && not _newton.result().converged)
        DUNE_THROW(NewtonNotConverged,
                   "NewtonTerminate::terminate(): Maximum iteration count reached");
      return converged;
    }

    virtual void setParameters(const ParameterTree& parameterTree) override
    {
      _maxIterations = parameterTree.get<unsigned int>("max_iterations", _maxIterations);
      _force_iteration = parameterTree.get<bool>("force_iteration", _force_iteration);
    }

  private:
    Newton& _newton;
    unsigned int _maxIterations = 40;
    bool _force_iteration = false;
  };
}

#endif
