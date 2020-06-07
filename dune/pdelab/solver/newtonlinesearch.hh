// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_HH
#define DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_HH

#include <dune/pdelab/solver/newtonerrors.hh>

namespace Dune::PDELab
{
  //! Abstract base class describing the line search interface
  template <typename Domain>
  class LineSearchInterface
  {
  public:
    //! Every abstract base class should have a virtual destructor
    virtual ~LineSearchInterface () {}

    //! Do line search
    virtual void lineSearch(Domain&, const Domain&) = 0;

    //! Set parameters
    virtual void setParameters(const ParameterTree&) = 0;
  };


  //! Class for simply updating the solution without line search
  template <typename Newton>
  class LineSearchNone : public LineSearchInterface<typename Newton::Domain>
  {
  public:
    using Domain = typename Newton::Domain;
    using Real = typename Newton::Real;

    LineSearchNone(Newton& newton) : _newton(newton) {}

    //! Do line search (in this case just update the solution)
    virtual void lineSearch(Domain& solution, const Domain& correction) override
    {
      solution.axpy(-1.0, correction);
      _newton.updateDefect(solution);
    }

    virtual void setParameters(const ParameterTree&) override {}

  private:
    Newton& _newton;
  };


  /** \brief Hackbusch-Reusken line search
   *
   * If the parameter line_search_accept_best is set through the setParameters
   * method this line search will simply return the best result even if it did
   * not converge.
   */
  template <typename Newton>
  class LineSearchHackbuschReusken : public LineSearchInterface<typename Newton::Domain>
  {
  public:
    using Domain = typename Newton::Domain;
    using Real = typename Newton::Real;

    LineSearchHackbuschReusken(Newton& newton) : _newton(newton) {}

    //! Do line search
    virtual void lineSearch(Domain& solution, const Domain& correction) override
    {
      if ((_newton.result().defect < _newton.getAbsoluteLimit())){
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

        if (not _acceptBest){
          solution = *_previousSolution;
          _newton.updateDefect(solution);
          DUNE_THROW(NewtonLineSearchError,
                     "NewtonLineSearch::line_search(): line search failed, "
                     "max iteration count reached, "
                     "defect did not improve enough");
        }
        else{
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


    /* \brief Set parameters
     *
     * Possible parameters are:
     *
     * - line_search_max_iterations: Maximum number of line search iterations.
     *
     * - line_search_damping_factor: Multiplier to line search parameter after each iteration.
     *
     * - line_search_accept_best: Accept the best line search parameter if
     *   there was any improvement, even if the convergence criterion was not
     *   reached.
     */
    virtual void setParameters(const ParameterTree& parameterTree) override
    {
      _lineSearchMaxIterations = parameterTree.get<unsigned int>("line_search_max_iterations",
                                                                 _lineSearchMaxIterations);
      _lineSearchDampingFactor = parameterTree.get<Real>("line_search_damping_factor",
                                                         _lineSearchDampingFactor);
      _acceptBest = parameterTree.get<bool>("line_search_accept_best",
                                            _acceptBest);
    }

  private:
    Newton& _newton;
    std::shared_ptr<Domain> _previousSolution;

    // Line search parameters
    unsigned int _lineSearchMaxIterations = 10;
    Real _lineSearchDampingFactor = 0.5;
    bool _acceptBest = false;
  };

  //! Flags for different line search strategies
  enum class LineSearchStrategy
  {
    noLineSearch,
    hackbuschReusken
  };

  // we put this into an emty namespace, so that we don't violate the one-definition-rule
  namespace {
    /** \brief Get a LineSearchStrategy from a string identifier
     *
     * \param name Identifier used to pick LineSearchStrategy
     *
     * Possible values for name: "noLineSearch", "hackbuschReusken"
     */
    inline
    LineSearchStrategy lineSearchStrategyFromString (const std::string& name)
    {
      if (name == "noLineSearch")
        return LineSearchStrategy::noLineSearch;
      if (name == "hackbuschReusken")
        return LineSearchStrategy::hackbuschReusken;
      DUNE_THROW(Exception,"Unkown line search strategy: " << name);
    }
  }


  /** \brief Get a pointer to a line search
   *
   * \tparam Newton A Newton solver
   *
   * \param newton Newton solver object

   * \param name Identifier to choose line search. Possible values:
   * - "noLineSearch": Return pointer to LineSearchNone
   * - "hackbuschReusken": Return pointer to LineSearchHackbuschReusken
   */
  template <typename Newton>
  std::shared_ptr<LineSearchInterface<typename Newton::Domain>>
  getLineSearch(Newton& newton, const std::string& name)
  {
    auto strategy = lineSearchStrategyFromString(name);
    if (strategy == LineSearchStrategy::noLineSearch){
      auto lineSearch = std::make_shared<LineSearchNone<Newton>> (newton);
      return lineSearch;
    }
    if (strategy == LineSearchStrategy::hackbuschReusken){
      auto lineSearch = std::make_shared<LineSearchHackbuschReusken<Newton>> (newton);
      return lineSearch;
    }
    DUNE_THROW(Exception,"Unkown line search strategy");
  }

} // namespace Dune::PDELab

#endif
