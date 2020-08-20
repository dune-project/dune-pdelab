// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_LINESEARCH_HH
#define DUNE_PDELAB_SOLVER_LINESEARCH_HH


namespace Dune::PDELab
{

  class LineSearchError : public Exception {};

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
  template <typename Solver>
  class LineSearchNone : public LineSearchInterface<typename Solver::Domain>
  {
  public:
    using Domain = typename Solver::Domain;
    using Real = typename Solver::Real;

    LineSearchNone(Solver& solver) : _solver(solver) {}

    //! Do line search (in this case just update the solution)
    virtual void lineSearch(Domain& solution, const Domain& correction) override
    {
      solution.axpy(-1.0, correction);
      _solver.updateDefect(solution);
    }

    virtual void setParameters(const ParameterTree&) override {}

  private:
    Solver& _solver;
  };


  /** \brief Hackbusch-Reusken line search
   *
   * If the parameter line_search_accept_best is set through the setParameters
   * method this line search will simply return the best result even if it did
   * not converge.
   */
  template <typename Solver>
  class LineSearchHackbuschReusken : public LineSearchInterface<typename Solver::Domain>
  {
  public:
    using Domain = typename Solver::Domain;
    using Real = typename Solver::Real;

    LineSearchHackbuschReusken(Solver& solver, bool forceAcceptBest = false) :
      _solver(solver), _forceAcceptBest(forceAcceptBest) {}

    //! Do line search
    virtual void lineSearch(Domain& solution, const Domain& correction) override
    {
      if ((_solver.result().defect < _solver.getAbsoluteLimit())){
        solution.axpy(-1.0, correction);
        _solver.updateDefect(solution);
        return;
      }

      auto verbosity = _solver.getVerbosityLevel();

      if (verbosity >= 4)
        std::cout << "      Performing line search..." << std::endl;
      Real lambda = 1.0;
      Real bestLambda = 0.0;
      Real bestDefect = _solver.result().defect;
      Real previousDefect = _solver.result().defect;
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
        _solver.updateDefect(solution);
        if (verbosity >= 4){
          if (not std::isfinite(_solver.result().defect))
            std::cout << "          NaNs detected" << std::endl;
        }       // ignore NaNs and try again with lower lambda

        if (_solver.result().defect <= (1.0 - lambda/4) * previousDefect){
          if (verbosity >= 4)
            std::cout << "          line search converged" << std::endl;
          converged = true;
          break;
        }

        if (_solver.result().defect < bestDefect){
          bestDefect = _solver.result().defect;
          bestLambda = lambda;
        }

        lambda *= _lineSearchDampingFactor;
        solution = *_previousSolution;
      }

      if (not converged){
        if (verbosity >= 4)
          std::cout << "          max line search iterations exceeded" << std::endl;

        if (not (_acceptBest or _forceAcceptBest)){
          solution = *_previousSolution;
          _solver.updateDefect(solution);
          DUNE_THROW( LineSearchError,
                     "LineSearch::lineSearch(): line search failed, "
                     "max iteration count reached, "
                     "defect did not improve enough");
        }
        else{
          if (bestLambda == 0.0){
            solution = *_previousSolution;
            _solver.updateDefect(solution);
            DUNE_THROW(LineSearchError,
                       "LineSearch::lineSearch(): line search failed, "
                       "max iteration count reached, "
                       "defect did not improve in any of the iterations");
          }
          if (bestLambda != lambda){
            solution = *_previousSolution;
            solution.axpy(-bestLambda, correction);
            _solver.updateDefect(solution);
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
     * - MaxIterations: Maximum number of line search iterations.
     *
     * - DampingFactor: Multiplier to line search parameter after each iteration.
     *
     * - AcceptBest: Accept the best line search parameter if
     *   there was any improvement, even if the convergence criterion was not
     *   reached.
     */
    virtual void setParameters(const ParameterTree& parameterTree) override
    {
      _lineSearchMaxIterations = parameterTree.get<unsigned int>("MaxIterations",
                                                                 _lineSearchMaxIterations);
      _lineSearchDampingFactor = parameterTree.get<Real>("DampingFactor",
                                                         _lineSearchDampingFactor);
      _acceptBest = parameterTree.get<bool>("AcceptBest", _acceptBest);
    }

  private:
    Solver& _solver;
    std::shared_ptr<Domain> _previousSolution;

    // Line search parameters
    unsigned int _lineSearchMaxIterations = 10;
    Real _lineSearchDampingFactor = 0.5;
    bool _acceptBest = false;
    bool _forceAcceptBest;
  };

  //! Flags for different line search strategies
  enum class LineSearchStrategy
  {
    noLineSearch,
    hackbuschReusken,
    hackbuschReuskenAcceptBest
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
      if (name == "hackbuschReuskenAcceptBest")
        return LineSearchStrategy::hackbuschReuskenAcceptBest;
      DUNE_THROW(Exception,"Unkown line search strategy: " << name);
    }
  }


  /** \brief fectory function to create an instace of a line-search
   *
   * \tparam Solver A solver
   *
   * \param solver Solver object

   * \param name Identifier to choose line search. Possible values:
   * - "noLineSearch": Return pointer to LineSearchNone
   * - "hackbuschReusken": Return pointer to LineSearchHackbuschReusken
   */
  template <typename Solver>
  std::shared_ptr<LineSearchInterface<typename Solver::Domain>>
  createLineSearch(Solver& solver, LineSearchStrategy strategy)
  {
    if (strategy == LineSearchStrategy::noLineSearch){
      auto lineSearch = std::make_shared<LineSearchNone<Solver>> (solver);
      return lineSearch;
    }
    if (strategy == LineSearchStrategy::hackbuschReusken){
      auto lineSearch = std::make_shared<LineSearchHackbuschReusken<Solver>> (solver);
      return lineSearch;
    }
    if (strategy == LineSearchStrategy::hackbuschReuskenAcceptBest){
      auto lineSearch = std::make_shared<LineSearchHackbuschReusken<Solver>> (solver, true);
      std::cout << "Warning: linesearch hackbuschReuskenAcceptBest is deprecated and will be removed after PDELab 2.7.\n"
                << "         Please use 'hackbuschReusken' and add the parameter 'LineSearchAcceptBest : true'";
      return lineSearch;
    }
    DUNE_THROW(Exception,"Unkown line search strategy");
  }

} // namespace Dune::PDELab

#endif
