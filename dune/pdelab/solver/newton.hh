// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_NEWTON_HH
#define DUNE_PDELAB_SOLVER_NEWTON_HH

#include <dune/common/exceptions.hh>
#include <dune/common/ios_state.hh>

#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/solver/newtonerrors.hh>
#include <dune/pdelab/solver/linesearch.hh>
#include <dune/pdelab/solver/terminate.hh>
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


  /** \brief Newton solver for solving non-linear problems
   *
   * - The line search and the termination criterion can be changed at runtime
   *   by the setTerminate() and the setLineSearch() methods.
   *
   * - If Newton is created using the default parameters it is an inexact
   *   Newton since the default reduction for the linear systems is quite
   *   low. You can change this through setMinLinearReduction()
   *
   * \tparam GridOperator_ Grid operator for evaluation of resdidual and Jacobian
   * \tparam LinearSolver_ Solver backend for solving linear system of equations
   */
  template <typename GridOperator_, typename LinearSolver_>
  class NewtonMethod
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
        DUNE_THROW(NewtonError, "NewtonMethod::result() called before NewtonMethod::solve()");
      return _result;
    }

    virtual void prepareStep(Domain& solution)
    {
      _reassembled = false;
      if (_result.defect/_previousDefect > _reassembleThreshold){
        if (_hangingNodeModifications){
          auto dirichletValues = solution;
          // Set all non dirichlet values to zero
          Dune::PDELab::set_shifted_dofs(_gridOperator.localAssembler().trialConstraints(), 0.0, dirichletValues);
          // Set all constrained DOFs to zero in solution
          Dune::PDELab::set_constrained_dofs(_gridOperator.localAssembler().trialConstraints(), 0.0, solution);
          // Copy correct Dirichlet values back into solution vector
          Dune::PDELab::copy_constrained_dofs(_gridOperator.localAssembler().trialConstraints(), dirichletValues, solution);
          // Interpolate periodic constraints / hanging nodes
          _gridOperator.localAssembler().backtransform(solution);
        }
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
                   "NewtonMethod::linearSolve(): Linear solver did not converge "
                   "in " << _linearSolver.result().iterations << " iterations");
      if (_verbosity >= 4)
        std::cout << "          linear solver iterations:     "
                  << std::setw(12) << _linearSolver.result().iterations << std::endl
                  << "          linear defect reduction:      "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << _linearSolver.result().reduction << std::endl;
    }

    //! Solve the nonlinear problem using solution as initial guess and for storing the result
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
      while (not _terminate->terminate()){
        if(_verbosity >= 3)
          std::cout << "  Newton iteration " << _result.iterations
                    << " --------------------------------" << std::endl;

        //=============
        // Prepare step
        //=============
        auto start = Clock::now();
        try{
          prepareStep(solution);
        }
        catch (...)
        {
          // Keep track of statistics when the method fails. We record
          // independently the time spent in non-converging attempts.
          // Check OneStepMethod to see how these data are propagated.
          auto end = Clock::now();
          assembler_time += end-start;
          _result.assembler_time = to_seconds(assembler_time);
          throw;
        }
        auto end = Clock::now();
        assembler_time += end -start;
        _result.assembler_time = to_seconds(assembler_time);

        // Store defect
        _previousDefect = _result.defect;

        //====================
        // Solve linear system
        //====================
        start = Clock::now();
        try{
          linearSolve();
        }
        catch (...)
        {
          // Separately catch statistics for linear solver failures.
          end = Clock::now();
          linear_solver_time += end-start;
          _result.linear_solver_time = to_seconds(linear_solver_time);
          _result.linear_solver_iterations = _linearSolver.result().iterations;
          throw;
        }
        end = Clock::now();
        linear_solver_time += end -start;
        _result.linear_solver_time = to_seconds(linear_solver_time);
        _result.linear_solver_iterations = _linearSolver.result().iterations;

        //===================================
        // Do line search and update solution
        //===================================
        start = Clock::now();
        _lineSearch->lineSearch(solution, _correction);
        // lineSearch(solution);
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
    virtual void updateDefect(Domain& solution)
    {
      if (_hangingNodeModifications){
        auto dirichletValues = solution;
        // Set all non dirichlet values to zero
        Dune::PDELab::set_shifted_dofs(_gridOperator.localAssembler().trialConstraints(), 0.0, dirichletValues);
        // Set all constrained DOFs to zero in solution
        Dune::PDELab::set_constrained_dofs(_gridOperator.localAssembler().trialConstraints(), 0.0, solution);
        // Copy correct Dirichlet values back into solution vector
        Dune::PDELab::copy_constrained_dofs(_gridOperator.localAssembler().trialConstraints(), dirichletValues, solution);
        // Interpolate periodic constraints / hanging nodes
        _gridOperator.localAssembler().backtransform(solution);
      }

      _residual = 0.0;
      _gridOperator.residual(solution, _residual);

      // Use the maximum norm as a stopping criterion. This helps loosen the tolerance
      // when solving for stationary solutions of nonlinear time-dependent problems.
      // The default is to use the Euclidean norm (in the else-block) as before
      if (_useMaxNorm){
        auto rankMax = Backend::native(_residual).infinity_norm();
        _result.defect = _gridOperator.testGridFunctionSpace().gridView().comm().max(rankMax);
      }
      else
        _result.defect =  _linearSolver.norm(_residual);
    }

    //! Set how much output you get
    void setVerbosityLevel(unsigned int verbosity)
    {
      if (_gridOperator.trialGridFunctionSpace().gridView().comm().rank()>0)
        _verbosity = 0;
      else
        _verbosity = verbosity;
    }

    //! Get verbosity level
    unsigned int getVerbosityLevel() const
    {
      return _verbosity;
    }

    //! Set reduction Newton needs to achieve
    void setReduction(Real reduction)
    {
      _reduction = reduction;
    }

    //! Get reduction
    Real getReduction() const
    {
      return _reduction;
    }

    //! Set absolute convergence limit
    void setAbsoluteLimit(Real absoluteLimit)
    {
      _absoluteLimit = absoluteLimit;
    }

    Real getAbsoluteLimit() const
    {
      return _absoluteLimit;
    }

    //! Set whether the jacobian matrix should be kept across calls to apply().
    void setKeepMatrix(bool b)
    {
      _keepMatrix = b;
    }

    //! Set whether to use the maximum norm for stopping criteria.
    void setUseMaxNorm(bool b)
    {
      _useMaxNorm = b;
    }

    //! Does the problem have hanging nodes
    void setHangingNodeModifications(bool b)
    {
      _hangingNodeModifications = b;
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

    /**\brief Set the minimal reduction in the linear solver
     *
     * \note with minLinearReduction > 0, the linear reduction will be
     * determined as mininum of the minLinearReduction and the linear reduction
     * needed to achieve second order Newton convergence. (As long as you are
     * not using a fixed linear reduction)
     */
    void setMinLinearReduction(Real minLinearReduction)
    {
      _minLinearReduction = minLinearReduction;
    }

    /** \brief Set wether to use a fixed reduction in the linear solver
     *
     * \note If fixedLinearReduction is true, the linear reduction rate will
     *  always be fixed to minLinearReduction.
     */
    void setFixedLinearReduction(bool fixedLinearReduction)
    {
      _fixedLinearReduction = fixedLinearReduction;
    }

    /** \brief Set a threshold, when the linear operator is reassembled
     *
     * We allow to keep the linear operator over several newton iterations. If
     * the reduction in the newton drops below a given threshold the linear
     * operator is reassembled to ensure convergence.
     */
    void setReassembleThreshold(Real reassembleThreshold)
    {
      _reassembleThreshold = reassembleThreshold;
    }

    /** \brief Interpret a parameter tree as a set of options for the newton solver
     *
     *  Possible parameters:
     *
     *  example configuration:
     *
     *  \code
     *  [newton_parameters]
     *  ReassembleThreshold = 0.1
     *  AbsoluteLimit = 1e-6
     *  Reduction = 1e-4
     *  MinLinearReduction = 1e-3
     *  MaxIterations = 15
     *  LineSearchDampingFactor = 0.7
     *  \endcode
     *
     *  and invocation in the code:
     *  \code
     *  newton.setParameters(param.sub("NewtonParameters"));
     *  \endcode
     *
     *  This can also be used to set single parameters like this
     *
     *  \code
     *  Dune::ParameterTree ptree;
     *  ptree["verbosity"] = "4";
     *  newton.setParameters(ptree);
     *  \endcode
     */
    void setParameters(const ParameterTree& parameterTree){
      _verbosity = parameterTree.get("VerbosityLevel", _verbosity);
      _reduction = parameterTree.get("Reduction", _reduction);
      _absoluteLimit = parameterTree.get("AbsoluteLimit", _absoluteLimit);
      _keepMatrix = parameterTree.get("KeepMatrix", _keepMatrix);
      _useMaxNorm = parameterTree.get("UseMaxNorm", _useMaxNorm);
      _hangingNodeModifications = parameterTree.get("HangingNodeModifications", _hangingNodeModifications);
      _minLinearReduction = parameterTree.get("MinLinearReduction", _minLinearReduction);
      _fixedLinearReduction = parameterTree.get("FixedLinearReduction", _fixedLinearReduction);
      _reassembleThreshold = parameterTree.get("ReassembleThreshold", _reassembleThreshold);

      // first create the linesearch, depending on the parameter
      std::string lineSearchStrategy = parameterTree.get("LineSearchStrategy","hackbuschReusken");
      auto strategy = lineSearchStrategyFromString(lineSearchStrategy);
      _lineSearch = createLineSearch(*this, strategy);

      // now set parameters
      if (parameterTree.hasSub("Terminate")){
        _terminate->setParameters(parameterTree.sub("Terminate"));
      }
      else{
        ParameterTree terminateTree;
        terminateTree["MaxIterations"] = std::to_string(parameterTree.get("MaxIterations", 40));
        terminateTree["ForceIteration"] = std::to_string(parameterTree.get("ForceIteration", false));
        _terminate->setParameters(terminateTree);
      }
      if (parameterTree.hasSub("LineSearch")){
        _lineSearch->setParameters(parameterTree.sub("LineSearch"));
      }
      else{
        ParameterTree lineSearchTree;
        lineSearchTree["MaxIterations"] = std::to_string(parameterTree.get("LineSearchMaxIterations", 10));
        lineSearchTree["DampingFactor"] = std::to_string(parameterTree.get("LineSearchDampingFactor", 0.5));
        lineSearchTree["AcceptBest"] = std::to_string(parameterTree.get("LineSearchAcceptBest", false));
        _lineSearch->setParameters(lineSearchTree);
      }
    }

    //! Set the termination criterion
    void setTerminate(std::shared_ptr<TerminateInterface> terminate)
    {
      _terminate = terminate;
    }

    /**\brief Set the line search
     *
     * See getLineSearch() for already implemented line searches
     */
    void setLineSearch(std::shared_ptr<LineSearchInterface<Domain>> lineSearch)
    {
      _lineSearch = lineSearch;
    }

    //! Construct Newton using default parameters with default parameters
    /**
       in p
     */
    NewtonMethod(
      const GridOperator& gridOperator,
      LinearSolver& linearSolver)
      : _gridOperator(gridOperator)
      , _linearSolver(linearSolver)
      , _residual(gridOperator.testGridFunctionSpace())
      , _correction(gridOperator.trialGridFunctionSpace())
    {
      _terminate = std::make_shared<DefaultTerminate<NewtonMethod>> (*this);
      _lineSearch = createLineSearch(*this, LineSearchStrategy::hackbuschReusken);
    }

    //! Construct Newton passing a parameter tree
    NewtonMethod(
      const GridOperator& gridOperator,
      LinearSolver& linearSolver,
      const ParameterTree& parameterTree)
      : _gridOperator(gridOperator)
      , _linearSolver(linearSolver)
      , _residual(gridOperator.testGridFunctionSpace())
      , _correction(gridOperator.trialGridFunctionSpace())

    {
      _terminate = std::make_shared<DefaultTerminate<NewtonMethod>> (*this);
      setParameters(parameterTree);
    }

  private:
    const GridOperator& _gridOperator;
    LinearSolver& _linearSolver;

    // Vectors and Jacobi matrix we set up only once
    Range _residual;
    Domain _correction;
    std::shared_ptr<Jacobian> _jacobian;
    std::shared_ptr<Domain> _previousSolution;

    std::shared_ptr<TerminateInterface> _terminate;
    std::shared_ptr<LineSearchInterface<Domain>> _lineSearch;

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
    bool _useMaxNorm = false;

    // Special treatment if we have hanging nodes
    bool _hangingNodeModifications = false;

    // User parameters for prepareStep()
    Real _minLinearReduction = 1e-3;
    bool _fixedLinearReduction = false;
    Real _reassembleThreshold = 0.0;
  };

} // namespace Dune::PDELab

#endif
