// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_NEWTON_HH
#define DUNE_PDELAB_SOLVER_NEWTON_HH

#include <dune/common/exceptions.hh>
#include <dune/common/ios_state.hh>

#include <dune/pdelab/solver/utility.hh>

namespace Dune::PDELab
{
  namespace Impl
  {
    // Some SFinae magic to execute setReuse(bool) on a backend
    template<typename T1, typename = void>
    struct HasSetReuse
      : std::false_type
    {};

    template<typename T>
    struct HasSetReuse<T, decltype(std::declval<T>().setReuse(true), void())>
      : std::true_type
    {};

    template<typename T>
    inline void setLinearSystemReuse(T& solver_backend, bool reuse, std::true_type)
    {
      if (!solver_backend.getReuse() && reuse)
        dwarn << "WARNING: Newton needed to override your choice to reuse the linear system in order to work!" << std::endl;
    solver_backend.setReuse(reuse);
    }

    template<typename T>
    inline void setLinearSystemReuse(T&, bool, std::false_type)
    {}

    template<typename T>
    inline void setLinearSystemReuse(T& solver_backend, bool reuse)
    {
      setLinearSystemReuse(solver_backend, reuse, HasSetReuse<T>());
    }
  }

  // Exception classes used in NewtonSolver
  class NewtonError : public Exception {};
  class NewtonDefectError : public NewtonError {};
  class NewtonLinearSolverError : public NewtonError {};
  class NewtonLineSearchError : public NewtonError {};
  class NewtonNotConverged : public NewtonError {};

  template <typename GridOperator_, typename LinearSolver_>
  class Newton
  {
  public:
    //! Type of the grid operator
    using GridOperator = GridOperator_;

    //! Type of the linear solver
    using LinearSolver = LinearSolver_;

    //! Type of the domain (solution)
    using Domain = typename GridOperator::Traits::Domain;

    //! Type of the range (residual)
    using Range = typename GridOperator::Traits::Range;

    //! Type of the Jacobian matrix
    using Jacobian = typename GridOperator::Traits::Jacobian;

    //! Number type
    using Real = typename Dune::FieldTraits<typename Domain::ElementType>::real_type;

    //! Type of results
    using Result = PDESolverResult<Real>;

    //! Return results
    const Result& result() const
    {
      if (not _resultValid)
        DUNE_THROW(NewtonError, "NewtonSolver::result() called before NewtonSolver::solve()");
      return _result;
    }

    virtual bool terminate()
    {
      if (_force_iteration && _result.iterations == 0)
        return false;
      _result.converged = _result.defect < _absoluteLimit || _result.defect < _result.first_defect * _reduction;
      if (_result.iterations >= _maxIterations && not _result.converged)
        DUNE_THROW(NewtonNotConverged,
                   "NewtonTerminate::terminate(): Maximum iteration count reached");
      return _result.converged;
    }

    enum LineSearchStrategy
    {
      noLineSearch,
      hackbuschReusken,
      hackbuschReuskenAcceptBest
    };

    virtual void lineSearch(Domain& solution)
    {
      if ((_lineSearchStrategy == noLineSearch) || (_result.defect < _absoluteLimit)){
        solution.axpy(-1.0, _correction);
        updateDefect(solution);
        return;
      }

      if (_verbosity >= 4)
        std::cout << "      Performing line search..." << std::endl;
      Real lambda = 1.0;
      Real best_lambda = 0.0;
      Real best_defect = _result.defect;

      if (not _previousSolution)
        _previousSolution = std::make_shared<Domain>(solution);
      else
        *_previousSolution = solution;

      unsigned int i = 0;
      while (1){
        if (_verbosity >= 4)
          std::cout << "          trying line search damping factor:   "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << lambda
                    << std::endl;

            solution.axpy(-lambda, _correction);
            try {
              updateDefect(solution);
            }
            catch (NewtonDefectError&){
              if (_verbosity >= 4)
                std::cout << "          NaNs detected" << std::endl;
            }       // ignore NaNs and try again with lower lambda

            if (_result.defect <= (1.0 - lambda/4) * _previousDefect){
              if (_verbosity >= 4)
                std::cout << "          line search converged" << std::endl;
              break;
            }

            if (_result.defect < best_defect){
              best_defect = _result.defect;
              best_lambda = lambda;
            }

            if (++i >= _lineSearchMaxIterations){
              if (_verbosity >= 4)
                std::cout << "          max line search iterations exceeded" << std::endl;
              switch (_lineSearchStrategy){
              case hackbuschReusken:
                solution = *_previousSolution;
                updateDefect(solution);
                DUNE_THROW(NewtonLineSearchError,
                           "NewtonLineSearch::line_search(): line search failed, "
                           "max iteration count reached, "
                           "defect did not improve enough");
              case hackbuschReuskenAcceptBest:
                if (best_lambda == 0.0){
                  solution = *_previousSolution;
                  updateDefect(solution);
                  DUNE_THROW(NewtonLineSearchError,
                             "NewtonLineSearch::line_search(): line search failed, "
                             "max iteration count reached, "
                             "defect did not improve in any of the iterations");
                }
                if (best_lambda != lambda){
                  solution = *_previousSolution;
                  solution.axpy(-best_lambda, _correction);
                  updateDefect(solution);
                }
                break;
              case noLineSearch:
                break;
              }
                break;
            }
            lambda *= _lineSearchDampingFactor;
            solution = *_previousSolution;
      }
      if (_verbosity >= 4)
        std::cout << "          line search damping factor:   "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << lambda << std::endl;
    }

    virtual void prepareStep(const Domain& solution)
    {
      _reassembled = false;
      if (_result.defect/_previousDefect > _reassembleThreshold){
        if (_verbosity>=3)
              std::cout << "      Reassembling matrix..." << std::endl;
        *_jacobian = Real(0.0);
        _gridOperator.jacobian(solution, *_jacobian);
        _reassembled = true;
      }

      _linearReduction = _minLinearReduction;
      if (not _fixedLinearReduction){
        // Determine maximum defect, where Newton is converged.
        using std::min;
        using std::max;
        auto stop_defect = max(_result.first_defect*_reduction, _absoluteLimit);

        // To achieve second order convergence of newton we need a linear
        // reduction of at least current_defect^2/prev_defect^2.  For the
        // last newton step a linear reduction of
        // 1/10*end_defect/current_defect is sufficient for convergence.
        if (stop_defect/(10*_result.defect) > _result.defect*_result.defect/(_previousDefect*_previousDefect))
          _linearReduction = stop_defect/(10*_result.defect);
        else
          _linearReduction = min(_minLinearReduction, _result.defect*_result.defect/(_previousDefect*_previousDefect));
        }

        if (_verbosity >= 3)
          std::cout << "      requested linear reduction:       "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << _linearReduction << std::endl;
    }

    virtual void linearSolve()
    {
      if (_verbosity >= 4)
        std::cout << "      Solving linear system..." << std::endl;

      // If the jacobian was not reassembled we might save some work in the solver backend
      Impl::setLinearSystemReuse(_linearSolver, not _reassembled);

      // Solve the linear system
      _correction = 0.0;
      _linearSolver.apply(*_jacobian, _correction, _residual, _linearReduction);

      if (not _linearSolver.result().converged)
        DUNE_THROW(NewtonLinearSolverError,
                   "NewtonSolver::linearSolve(): Linear solver did not converge "
                   "in " << _linearSolver.result().iterations << " iterations");
      if (_verbosity >= 4)
        std::cout << "          linear solver iterations:     "
                  << std::setw(12) << _linearSolver.result().iterations << std::endl
                  << "          linear defect reduction:      "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << _linearSolver.result().reduction << std::endl;
    }

    //! Provide an inital guess and solve inplace
    virtual void apply(Domain& solution)
    {
      // Reset solver statistics
      _result.clear();
      _resultValid = true;

      // Store old ios flags (will be reset when this goes out of scope)
      ios_base_all_saver restorer(std::cout);

      // Prepare time measuring
      using Clock = std::chrono::steady_clock;
      using Duration = Clock::duration;
      auto assembler_time = Duration::zero();
      auto linear_solver_time = Duration::zero();
      auto line_search_time   = Duration::zero();
      auto to_seconds =
        [](Duration duration){
          return std::chrono::duration<double>(duration).count();
        };
      auto start_solve = Clock::now();

      //=========================
      // Calculate initial defect
      //=========================
      updateDefect(solution);
      _result.first_defect = _result.defect;
      _previousDefect = _result.defect;

      if (_verbosity >= 2)
        std::cout << "  Initial defect: "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << _result.defect << std::endl;

      //==========================
      // Calculate Jacobian matrix
      //==========================
      if (not _jacobian)
        _jacobian = std::make_shared<Jacobian>(_gridOperator);

      //=========================
      // Nonlinear iteration loop
      //=========================
      while (not terminate()){
        if(_verbosity >= 3)
          std::cout << "  Newton iteration " << _result.iterations
                    << " --------------------------------" << std::endl;

        //=============
        // Prepare step
        //=============
        auto start = Clock::now();
        prepareStep(solution);
        auto end = Clock::now();
        assembler_time += end -start;

        // Store defect
        _previousDefect = _result.defect;

        //====================
        // Solve linear system
        //====================
        start = Clock::now();
        linearSolve();
        end = Clock::now();
        linear_solver_time += end -start;

        //===================================
        // Do line search and update solution
        //===================================
        start = Clock::now();
        lineSearch(solution);
        end = Clock::now();
        line_search_time += end -start;

        //========================================
        // Store statistics and create some output
        //========================================
        _result.reduction = _result.defect/_result.first_defect;
        _result.iterations++;
        _result.conv_rate = std::pow(_result.reduction, 1.0/_result.iterations);

        if (_verbosity >= 3)
          std::cout << "      linear solver time:               "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << to_seconds(linear_solver_time) << std::endl
                    << "      defect reduction (this iteration):"
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << _result.defect/_previousDefect << std::endl
                    << "      defect reduction (total):         "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << _result.reduction << std::endl
                    << "      new defect:                       "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << _result.defect << std::endl;
        if (_verbosity == 2)
          std::cout << "  Newton iteration "
                    << std::setw(2)
                    << _result.iterations
                    << ".  New defect: "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << _result.defect
                    << ".  Reduction (this): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << _result.defect/_previousDefect
                    << ".  Reduction (total): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << _result.reduction
                    << std::endl;
      } // end while loop of nonlinear iterations

      auto end_solve = Clock::now();
      auto solve_time = end_solve - start_solve;
      _result.elapsed = to_seconds(solve_time);

      if (_verbosity == 1)
        std::cout << "  Newton converged after "
                  << std::setw(2)
                  << _result.iterations
                  << " iterations.  Reduction: "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << _result.reduction
                  << "   (" << std::setprecision(4) << _result.elapsed << "s)"
                  << std::endl;

      if (not _keepMatrix)
        _jacobian.reset();
    }

    //! Update _residual and defect in _result
    virtual void updateDefect(const Domain& solution)
    {
      _residual = 0.0;
      _gridOperator.residual(solution, _residual);
      _result.defect =  _linearSolver.norm(_residual);
      if(not std::isfinite(_result.defect))
        DUNE_THROW(NewtonDefectError,
                   "NewtonSolver::defect(): Non-linear defect is NaN or Inf");
    }

    //! Set how much output you get
    void setVerbosity(unsigned int verbosity)
    {
      if (_gridOperator.trialGridFunctionSpace().gridView().comm().rank()>0)
        _verbosity = 0;
      else
        _verbosity = verbosity;
    }

    //! Set reduction Newton needs to achieve
    void setReduction(Real reduction)
    {
      _reduction = reduction;
    }

    //! Set absolute convergence limit
    void setAbsoluteLimit(Real absoluteLimit)
    {
      _absoluteLimit = absoluteLimit;
    }

    //! Set whether the jacobian matrix should be kept across calls to apply().
    void setKeepMatrix(bool b)
    {
      _keepMatrix = b;
    }

    //! Return whether the jacobian matrix is kept across calls to apply().
    bool keepMatrix() const
    {
      return _keepMatrix;
    }

    //! Discard the stored Jacobian matrix.
    void discardMatrix()
    {
      if(_jacobian)
        _jacobian.reset();
    }

    //! Set maximum number of nonlinear Newton iterations
    void setMaxIterations(unsigned int maxit)
    {
      _maxIterations = maxit;
    }

    //! Set to true if Newton should be forced to do at least one iteration
    void setForceIteration(bool forceIteration)
    {
      _force_iteration = forceIteration;
    }

    //! Set line search strategy
    void setLineSearchStrategy(LineSearchStrategy strategy)
    {
      _lineSearchStrategy = strategy;
    }

    //! Set line search strategy
    void setLineSearchStrategy(std::string strategy)
    {
      _lineSearchStrategy = lineSearchStrategyFromName(strategy);
    }

    //! Set maximum amount of line search iterations
    void setLineSearchMaxIterations(unsigned int maxit)
    {
      _lineSearchMaxIterations = maxit;
    }

    /**\brief Set damping factor in line search
     *
     * This will be used as multiplier if the line search did not yet converge.
     */
    void setLineSearchDampingFactor(Real dampingFactor)
    {
      _lineSearchDampingFactor = dampingFactor;
    }

    /**\brief set the minimal reduction in the linear solver

       \note with min_linear_reduction > 0, the linear reduction will be
       determined as mininum of the min_linear_reduction and the
       linear_reduction needed to achieve second order
       Newton convergence. */
    void setMinLinearReduction(Real minLinearReduction)
    {
      _minLinearReduction = minLinearReduction;
    }

    /**\brief set a fixed reduction in the linear solver (overwrites setMinLinearReduction)

       \note with fixed_linear_reduction > 0, the linear reduction
       rate will always be fixed to min_linear_reduction. */
    void setFixedLinearReduction(bool fixedLinearReduction)
    {
      _fixedLinearReduction = fixedLinearReduction;
    }

    /**\brief set a threshold, when the linear operator is reassembled

       We allow to keep the linear operator over several newton
       iterations. If the reduction in the newton drops below a
       given threshold the linear operator is reassembled to ensure
       convergence.
    */
    void setReassembleThreshold(Real reassembleThreshold)
    {
      _reassembleThreshold = reassembleThreshold;
    }

    /** \brief Interpret a parameter tree as a set of options for the newton solver

        example configuration:

        \code
        [NewtonParameters]
        reassemble_threshold = 0.1
        line_search_max_iterations = 10
        max_iterations = 7
        absolute_limit = 1e-6
        reduction = 1e-4
        min_linear_reduction = 1e-3
        line_search_damping_factor  = 0.9
        \endcode

        and invocation in the code:
        \code
        newton.setParameters(param.sub("NewtonParameters"));
        \endcode
    */
    void setParameters(const ParameterTree& parameterTree){
        if (parameterTree.hasKey("verbosity"))
          setVerbosity(parameterTree.get<unsigned int>("verbosity"));
        if (parameterTree.hasKey("reduction"))
          setReduction(parameterTree.get<Real>("reduction"));
        if (parameterTree.hasKey("absolute_limit"))
          setAbsoluteLimit(parameterTree.get<Real>("absolute_limit"));
        if (parameterTree.hasKey("keeep_matrix"))
          setKeepMatrix(parameterTree.get<bool>("keep_matrix"));

        if (parameterTree.hasKey("max_iterations"))
          setMaxIterations(parameterTree.get<unsigned int>("max_iterations"));
        if (parameterTree.hasKey("force_iteration"))
          setForceIteration(parameterTree.get<bool>("force_iteration"));

        if (parameterTree.hasKey("line_search_strategy"))
          setLineSearchStrategy(parameterTree.get<std::string>("line_search_strategy"));
        if (parameterTree.hasKey("line_search_max_iterations"))
          setLineSearchMaxIterations(parameterTree.get<unsigned int>("line_search_max_iterations"));
        if (parameterTree.hasKey("line_search_damping_factor"))
          setLineSearchDampingFactor(parameterTree.get<Real>("line_search_damping_factor"));

        if (parameterTree.hasKey("min_linear_reduction"))
          setMinLinearReduction(parameterTree.get<Real>("min_linear_reduction"));
        if (parameterTree.hasKey("fixed_linear_reduction"))
          setFixedLinearReduction(parameterTree.get<bool>("fixed_linear_reduction"));
        if (parameterTree.hasKey("reassemble_threshold"))
          setReassembleThreshold(parameterTree.get<Real>("reassemble_threshold"));
    }

    Newton(
      const GridOperator& gridOperator,
      LinearSolver& linearSolver)
      : _gridOperator(gridOperator)
      , _linearSolver(linearSolver)
      , _residual(gridOperator.testGridFunctionSpace())
      , _correction(gridOperator.trialGridFunctionSpace())
    {}

    Newton(
      const GridOperator& gridOperator,
      LinearSolver& linearSolver,
      const ParameterTree& parameterTree)
      : _gridOperator(gridOperator)
      , _linearSolver(linearSolver)
      , _residual(gridOperator.testGridFunctionSpace())
      , _correction(gridOperator.trialGridFunctionSpace())

    {
      setParameters(parameterTree);
    }

  private:
    LineSearchStrategy lineSearchStrategyFromName (const std::string & s) {
      if (s == "noLineSearch")
        return noLineSearch;
      if (s == "hackbuschReusken")
        return hackbuschReusken;
      if (s == "hackbuschReuskenAcceptBest")
        return hackbuschReuskenAcceptBest;
      DUNE_THROW(Exception, "unknown line search strategy" << s);
    }

    const GridOperator& _gridOperator;
    LinearSolver& _linearSolver;

    // Vectors and Jacobi matrix we set up only once
    Range _residual;
    Domain _correction;
    shared_ptr<Jacobian> _jacobian;
    shared_ptr<Domain> _previousSolution;

    // Class for storing results
    Result _result;
    bool _resultValid = false; // result class only valid after calling apply
    Real _previousDefect = 0.0;

    // Remember if jacobian was reassembled in prepareStep
    bool _reassembled = true; // will be set in prepare step
    Real _linearReduction = 0.0; // will be set in prepare step

    // User parameters
    unsigned int _verbosity = 0;
    Real _reduction = 1e-8;
    Real _absoluteLimit = 1e-12;
    bool _keepMatrix = true;

    // User parameters for terminate()
    unsigned int _maxIterations = 40;
    bool _force_iteration = false;

    // User parameters for lineSearch()
    LineSearchStrategy _lineSearchStrategy = LineSearchStrategy::hackbuschReusken;
    int _lineSearchMaxIterations = 10;
    Real _lineSearchDampingFactor = 0.5;

    // User parameters for prepareStep()
    Real _minLinearReduction = 1e-3;
    bool _fixedLinearReduction = false;
    Real _reassembleThreshold = 0.0;
  };

} // namespace Dune::PDELab

#endif
