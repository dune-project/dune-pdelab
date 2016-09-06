#ifndef DUNE_PDELAB_STATIONARY_LINEARPROBLEM_HH
#define DUNE_PDELAB_STATIONARY_LINEARPROBLEM_HH

#include <iostream>

#include <dune/common/timer.hh>
#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/backend/solver.hh>

namespace Dune {
  namespace PDELab {

    //===============================================================
    // A class for solving linear stationary problems.
    // It assembles the matrix, computes the right hand side and
    // solves the problem.
    // This is only a first vanilla implementation which has to be improved.
    //===============================================================

    // Status information of linear problem solver
    template<class RFType>
    struct StationaryLinearProblemSolverResult : LinearSolverResult<RFType>
    {
      RFType first_defect;       // the first defect
      RFType defect;             // the final defect
      double assembler_time;     // Cumulative time for matrix assembly
      double linear_solver_time; // Cumulative time for linear solver
      int linear_solver_iterations; // Total number of linear iterations

      StationaryLinearProblemSolverResult()
        : first_defect(0.0)
        , defect(0.0)
        , assembler_time(0.0)
        , linear_solver_time(0.0)
        , linear_solver_iterations(0)
      {}

    };

    template<typename GO, typename LS, typename V>
    class StationaryLinearProblemSolver
    {
      typedef typename Dune::template FieldTraits<typename V::ElementType >::real_type Real;
      typedef typename GO::Traits::Jacobian M;
      typedef typename GO::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
      using W = Dune::PDELab::Backend::Vector<TrialGridFunctionSpace,typename V::ElementType>;
      typedef GO GridOperator;

    public:
      typedef StationaryLinearProblemSolverResult<double> Result;

      StationaryLinearProblemSolver(const GO& go, LS& ls, V& x, Real reduction, Real min_defect = 1e-99, int verbose=1)
        : _go(go)
        , _ls(ls)
        , _x(&x)
        , _reduction(reduction)
        , _min_defect(min_defect)
        , _hanging_node_modifications(false)
        , _keep_matrix(true)
        , _verbose(verbose)
      {}

      StationaryLinearProblemSolver (const GO& go, LS& ls, Real reduction, Real min_defect = 1e-99, int verbose=1)
        : _go(go)
        , _ls(ls)
        , _x()
        , _reduction(reduction)
        , _min_defect(min_defect)
        , _hanging_node_modifications(false)
        , _keep_matrix(true)
        , _verbose(verbose)
      {}

      //! Construct a StationaryLinearProblemSolver for the given objects and read parameters from a ParameterTree.
      /**
       * This constructor reads the parameter controlling its operation from a passed-in ParameterTree
       * instead of requiring the user to specify all of them as individual constructor parameters.
       * Currently the following parameters are read:
       *
       * Name                       | Default Value | Explanation
       * -------------------------- | ------------- | -----------
       * reduction                  |               | Required relative defect reduction
       * min_defect                 | 1e-99         | minimum absolute defect at which to stop
       * hanging_node_modifications | false         | perform required transformations for hanging nodes
       * keep_matrix                | true          | keep matrix between calls to apply() (but reassemble values every time)
       * verbosity                  | 1             | control amount of debug output
       *
       * Apart from reduction, all parameters have a default value and are optional.
       * The actual reduction for a call to apply() is calculated as r = max(reduction,min_defect/start_defect),
       * where start defect is the norm of the residual of x.
       */
      StationaryLinearProblemSolver(const GO& go, LS& ls, V& x, const ParameterTree& params)
        : _go(go)
        , _ls(ls)
        , _x(&x)
        , _reduction(params.get<Real>("reduction"))
        , _min_defect(params.get<Real>("min_defect",1e-99))
        , _hanging_node_modifications(params.get<bool>("hanging_node_modifications",false))
        , _keep_matrix(params.get<bool>("keep_matrix",true))
        , _verbose(params.get<int>("verbosity",1))
      {}

      //! Construct a StationaryLinearProblemSolver for the given objects and read parameters from a ParameterTree.
      /**
       * This constructor reads the parameter controlling its operation from a passed-in ParameterTree
       * instead of requiring the user to specify all of them as individual constructor parameters.
       * Currently the following parameters are read:
       *
       * Name                       | Default Value | Explanation
       * -------------------------- | ------------- | -----------
       * reduction                  |               | Required relative defect reduction
       * min_defect                 | 1e-99         | minimum absolute defect at which to stop
       * hanging_node_modifications | false         | perform required transformations for hanging nodes
       * keep_matrix                | true          | keep matrix between calls to apply() (but reassemble values every time)
       * verbosity                  | 1             | control amount of debug output
       *
       * Apart from reduction, all parameters have a default value and are optional.
       * The actual reduction for a call to apply() is calculated as r = max(reduction,min_defect/start_defect),
       * where start defect is the norm of the residual of x.
       */
      StationaryLinearProblemSolver(const GO& go, LS& ls, const ParameterTree& params)
        : _go(go)
        , _ls(ls)
        , _x()
        , _reduction(params.get<Real>("reduction"))
        , _min_defect(params.get<Real>("min_defect",1e-99))
        , _hanging_node_modifications(params.get<bool>("hanging_node_modifications",false))
        , _keep_matrix(params.get<bool>("keep_matrix",true))
        , _verbose(params.get<int>("verbosity",1))
      {}

      //! Set whether the solver should apply the necessary transformations for calculations on hanging nodes.
      void setHangingNodeModifications(bool b)
      {
        _hanging_node_modifications=b;
      }

      //! Return whether the solver performs the necessary transformations for calculations on hanging nodes.
      bool hangingNodeModifications() const
      {
        return _hanging_node_modifications;
      }

      //! Set whether the jacobian matrix should be kept across calls to apply().
      void setKeepMatrix(bool b)
      {
        _keep_matrix = b;
      }

      //! Return whether the jacobian matrix is kept across calls to apply().
      bool keepMatrix() const
      {
        return _keep_matrix;
      }

      const Result& result() const
      {
        return _res;
      }

      void apply(V& x, bool reuse_matrix = false) {
        _x = &x;
        apply(reuse_matrix);
      }

      void apply (bool reuse_matrix = false)
      {
        Dune::Timer watch;
        double timing,assembler_time=0;

        // assemble matrix; optional: assemble only on demand!
        watch.reset();

        if (!_jacobian)
          {
            _jacobian = std::make_shared<M>(_go);
            timing = watch.elapsed();
            if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
              std::cout << "=== matrix setup (max) " << timing << " s" << std::endl;
            watch.reset();
            assembler_time += timing;
          }
        else if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
          std::cout << "=== matrix setup skipped (matrix already allocated)" << std::endl;

        if (_hanging_node_modifications)
          {
            Dune::PDELab::set_shifted_dofs(_go.localAssembler().trialConstraints(),0.0,*_x); // set hanging node DOFs to zero
            _go.localAssembler().backtransform(*_x); // interpolate hanging nodes adjacent to Dirichlet nodes
          }

        if (!reuse_matrix)
          {
            (*_jacobian) = Real(0.0);
            _go.jacobian(*_x,*_jacobian);
          }

        timing = watch.elapsed();
        // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
        if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
          {
            if (reuse_matrix)
              std::cout << "=== matrix assembly SKIPPED" << std::endl;
            else
              std::cout << "=== matrix assembly (max) " << timing << " s" << std::endl;
          }

        assembler_time += timing;

        // assemble residual
        watch.reset();

        W r(_go.testGridFunctionSpace(),0.0);
        _go.residual(*_x,r);  // residual is additive

        timing = watch.elapsed();
        // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
        if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
          std::cout << "=== residual assembly (max) " << timing << " s" << std::endl;
        assembler_time += timing;
        _res.assembler_time = assembler_time;

        auto defect = _ls.norm(r);

        // compute correction
        watch.reset();
        V z(_go.trialGridFunctionSpace(),0.0);
        auto red = std::max(_reduction,_min_defect/defect);
        if (_go.trialGridFunctionSpace().gridView().comm().rank()==0)
        {
          std::cout << "=== solving (reduction: " << red << ") ";
          if (_verbose>=1)
            std::cout << std::flush;
          else
            std::cout << std::endl;
        }
        _ls.apply(*_jacobian,z,r,red); // solver makes right hand side consistent
        _linear_solver_result = _ls.result();
        timing = watch.elapsed();
        // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
        if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
          std::cout << timing << " s" << std::endl;
        _res.linear_solver_time = timing;

        _res.converged = _linear_solver_result.converged;
        _res.iterations = _linear_solver_result.iterations;
        _res.elapsed = _linear_solver_result.elapsed;
        _res.reduction = _linear_solver_result.reduction;
        _res.conv_rate = _linear_solver_result.conv_rate;
        _res.first_defect = static_cast<double>(defect);
        _res.defect = static_cast<double>(defect)*_linear_solver_result.reduction;
        _res.linear_solver_iterations = _linear_solver_result.iterations;

        // and update
        if (_hanging_node_modifications)
          Dune::PDELab::set_shifted_dofs(_go.localAssembler().trialConstraints(),0.0,*_x); // set hanging node DOFs to zero
        *_x -= z;
        if (_hanging_node_modifications)
          _go.localAssembler().backtransform(*_x); // interpolate hanging nodes adjacent to Dirichlet nodes

        if (!_keep_matrix)
          _jacobian.reset();
      }

      //! Discard the stored Jacobian matrix.
      void discardMatrix()
      {
        if(_jacobian)
          _jacobian.reset();
      }

      const Dune::PDELab::LinearSolverResult<double>& ls_result() const{
        return _linear_solver_result;
      }

      Real reduction() const
      {
        return _reduction;
      }

      void setReduction(Real reduction)
      {
        _reduction = reduction;
      }


    private:
      const GO& _go;
      LS& _ls;
      V* _x;
      shared_ptr<M> _jacobian;
      Real _reduction;
      Real _min_defect;
      Dune::PDELab::LinearSolverResult<double> _linear_solver_result;
      Result _res;
      bool _hanging_node_modifications;
      bool _keep_matrix;
      int _verbose;
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_STATIONARY_LINEARPROBLEM_HH
