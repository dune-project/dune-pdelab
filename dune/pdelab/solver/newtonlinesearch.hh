// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_HH
#define DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_HH

#include <dune/pdelab/solver/newtonerrors.hh>

namespace Dune::PDELab
{
  template <typename Domain>
  class LineSearchInterface
  {
  public:
    //! Every abstract base class should have a virtual destructor
    virtual ~LineSearchInterface () {}

    virtual void lineSearch(Domain&, const Domain&) = 0;

    virtual void setParameters(const ParameterTree&) = 0;
  };


  template <typename Newton>
  class DefaultLineSearch : public LineSearchInterface<typename Newton::Domain>
  {
  public:
    using Domain = typename Newton::Domain;
    using Real = typename Newton::Real;

    DefaultLineSearch(Newton& newton) : _newton(newton) {}

    virtual void lineSearch(Domain& solution, const Domain& correction) override
    {
      if ((_lineSearchStrategy == noLineSearch) || (_newton.result().defect < _newton.getAbsoluteLimit())){
        solution.axpy(-1.0, correction);
        _newton.updateDefect(solution);
        return;
      }

      auto verbosity = _newton.getVerbosityLevel();

      if (verbosity >= 4)
        std::cout << "      Performing line search..." << std::endl;
      Real lambda = 1.0;
      Real bestLambda = 0.0;
      Real bestDefect = _newton.result().defect;
      Real previousDefect = _newton.result().defect;
      bool converged = false;

      if (not _previousSolution)
        _previousSolution = std::make_shared<Domain>(solution);
      else
        *_previousSolution = solution;

      for (unsigned int iteration = 0; iteration < _lineSearchMaxIterations; ++iteration){
        if (verbosity >= 4)
          std::cout << "          trying line search damping factor:   "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << lambda
                    << std::endl;

        solution.axpy(-lambda, correction);
        _newton.updateDefect(solution);
        if (verbosity >= 4){
          if (not std::isfinite(_newton.result().defect))
            std::cout << "          NaNs detected" << std::endl;
        }       // ignore NaNs and try again with lower lambda

        if (_newton.result().defect <= (1.0 - lambda/4) * previousDefect){
          if (verbosity >= 4)
            std::cout << "          line search converged" << std::endl;
          converged = true;
          break;
        }

        if (_newton.result().defect < bestDefect){
          bestDefect = _newton.result().defect;
          bestLambda = lambda;
        }

        lambda *= _lineSearchDampingFactor;
        solution = *_previousSolution;
      }

      if (not converged){
        if (verbosity >= 4)
          std::cout << "          max line search iterations exceeded" << std::endl;

        switch (_lineSearchStrategy){
        case hackbuschReusken:
          solution = *_previousSolution;
          _newton.updateDefect(solution);
          DUNE_THROW(NewtonLineSearchError,
                     "NewtonLineSearch::line_search(): line search failed, "
                     "max iteration count reached, "
                     "defect did not improve enough");
        case hackbuschReuskenAcceptBest:
          if (bestLambda == 0.0){
            solution = *_previousSolution;
            _newton.updateDefect(solution);
            DUNE_THROW(NewtonLineSearchError,
                       "NewtonLineSearch::line_search(): line search failed, "
                       "max iteration count reached, "
                       "defect did not improve in any of the iterations");
          }
          if (bestLambda != lambda){
            solution = *_previousSolution;
            solution.axpy(-bestLambda, correction);
            _newton.updateDefect(solution);
            converged = true;
          }
        }
      }

      if (converged)
        if (verbosity >= 4)
          std::cout << "          line search damping factor:   "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << lambda << std::endl;
    }

    virtual void setParameters(const ParameterTree& parameterTree)
    {
      if (parameterTree.hasKey("line_search_strategy")){
        auto strategy = parameterTree.get<std::string>("line_search_strategy");
        _lineSearchStrategy = lineSearchStrategyFromName(strategy);
      }
      _lineSearchMaxIterations = parameterTree.get<unsigned int>("line_search_max_iterations",
                                                                 _lineSearchMaxIterations);
      _lineSearchDampingFactor = parameterTree.get<Real>("line_search_damping_factor",
                                                         _lineSearchDampingFactor);
    }

  private:
    enum LineSearchStrategy
    {
      noLineSearch,
      hackbuschReusken,
      hackbuschReuskenAcceptBest
    };

    LineSearchStrategy lineSearchStrategyFromName (const std::string & s) {
      if (s == "noLineSearch")
        return noLineSearch;
      if (s == "hackbuschReusken")
        return hackbuschReusken;
      if (s == "hackbuschReuskenAcceptBest")
        return hackbuschReuskenAcceptBest;
      DUNE_THROW(Exception, "unknown line search strategy" << s);
    }

    Newton& _newton;
    std::shared_ptr<Domain> _previousSolution;

    // Line search parameters
    LineSearchStrategy _lineSearchStrategy = LineSearchStrategy::hackbuschReusken;
    unsigned int _lineSearchMaxIterations = 10;
    Real _lineSearchDampingFactor = 0.5;
  };
} // namespace Dune::PDELab

#endif
