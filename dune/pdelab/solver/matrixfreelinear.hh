// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_SOLVER_MATRIXFREELINEAR_HH
#define DUNE_PDELAB_SOLVER_MATRIXFREELINEAR_HH

//
// Note: This is basically a copy of the StationaryLinearProblemSolverResult. A
// smarter implementation that avoids code duplication needs a change in the
// solver backends.
//

#include <iostream>

#include <dune/common/timer.hh>
#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/backend/solver.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

namespace Dune {
  namespace PDELab {

    /** \brief Solve linear problems using matrix-free solvers
     */
    template<typename GO, typename LS, typename V>
    class MatrixFreeStationaryLinearProblemSolver
    {
      /** \brief Field value type */
      typedef typename Dune::template FieldTraits<typename V::ElementType >::real_type Real;
      /** \brief Typo of Jacobian matrix */
      typedef typename GO::Traits::Jacobian M;
      /** \brief Trial grid function space */
      typedef typename GO::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
      /** \brief Vector backend for trial grid function space */
      using W = Dune::PDELab::Backend::Vector<TrialGridFunctionSpace,typename V::ElementType>;
      /** \brief Grid operator to solve for */
      typedef GO GridOperator;

    public:
      /** \brief Class holding results */
      typedef StationaryLinearProblemSolverResult<double> Result;

      /** \brief Construct new solver instance with initial guess of solution
       *
       * \param[in] go Grid operator representing the linear problem
       *            to be solved
       * \param[in] ls Linear solver backend
       * \param[inout] x Solution vector, start with initial guess
       * \param[in] reduction Tolerance, target relative residual reduction
       * \param[in] min_defect Minimal absolute residual
       * \param[in] verbose Verbosity level
       */
      MatrixFreeStationaryLinearProblemSolver(const GO& go, LS& ls, V& x, Real reduction, Real min_defect = 1e-99, int verbose=1)
        : _go(go)
        , _ls(ls)
        , _x(&x)
        , _reduction(reduction)
        , _min_defect(min_defect)
        , _hanging_node_modifications(false)
        , _keep_matrix(true)
        , _verbose(verbose)
      {}

      /** \brief Construct new solver instance
       *
       * \param[in] go Grid operator representing the linear problem
       *            to be solved
       * \param[in] idgo Grid operator implementing the block-diagonal
       *            inverse of go, used for preconditioning
       * \param[in] ls Linear solver backend
       * \param[in] reduction Tolerance, target relative residual reduction
       * \param[in] min_defect Minimal absolute residual
       * \param[in] verbose Verbosity level
       */
      MatrixFreeStationaryLinearProblemSolver (const GO& go, LS& ls, Real reduction, Real min_defect = 1e-99, int verbose=1)
        : _go(go)
        , _ls(ls)
        , _x()
        , _reduction(reduction)
        , _min_defect(min_defect)
        , _hanging_node_modifications(false)
        , _keep_matrix(true)
        , _verbose(verbose)
      {}

      //! Construct a MatrixFreeStationaryLinearProblemSolver for the given objects and read parameters from a ParameterTree.
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
      MatrixFreeStationaryLinearProblemSolver(const GO& go, LS& ls, V& x, const ParameterTree& params)
        : _go(go)
        , _ls(ls)
        , _x(&x)
        , _reduction(params.get<Real>("reduction"))
        , _min_defect(params.get<Real>("min_defect",1e-99))
        , _hanging_node_modifications(params.get<bool>("hanging_node_modifications",false))
        , _keep_matrix(params.get<bool>("keep_matrix",true))
        , _verbose(params.get<int>("verbosity",1))
      {}

      //! Construct a MatrixFreeStationaryLinearProblemSolver for the given objects and read parameters from a ParameterTree.
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
      MatrixFreeStationaryLinearProblemSolver(const GO& go, LS& ls, const ParameterTree& params)
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

      //! Return result object
      const Result& result() const
      {
        return _res;
      }

      //! Provide initial guess and solve linear problem
      void apply(V& x, bool reuse_matrix = false) {
        _x = &x;
        apply(reuse_matrix);
      }

      //! Solve linear problem
      void apply (bool reuse_matrix = false)
      {
        Dune::Timer watch;
        double timing,assembler_time=0;
        int rank=_go.trialGridFunctionSpace().gridView().comm().rank();

        // Assemble matrix if necessary
        watch.reset();
        if (rank==0 && _verbose>=1){
          std::cout << "=== matrix setup not required for matrix free solvers" << std::endl;
        }
        if (_hanging_node_modifications)
          {
            Dune::PDELab::set_shifted_dofs(_go.localAssembler().trialConstraints(),0.0,*_x); // set hanging node DOFs to zero
            _go.localAssembler().backtransform(*_x); // interpolate hanging nodes adjacent to Dirichlet nodes
          }

        timing = watch.elapsed();
        assembler_time += timing;

        // assemble residual
        watch.reset();

        W r(_go.testGridFunctionSpace(),0.0);
        _go.residual(*_x,r);  // residual is additive

        timing = watch.elapsed();

        if (rank==0 && _verbose>=1)
          std::cout << "=== residual assembly (max) " << timing << " s" << std::endl;
        assembler_time += timing;
        _res.assembler_time = assembler_time;

        auto defect = _ls.norm(r);

        // compute correction
        watch.reset();
        V z(_go.trialGridFunctionSpace(),0.0);
        auto red = std::max(_reduction,_min_defect/defect);
        if (rank==0 && _verbose>=1)
        {
          std::cout << "=== solving (reduction: " << red << ") ";
          if (_verbose>=1)
            std::cout << std::flush;
          else
            std::cout << std::endl;
        }

        // TODO: For nonlinear problems we need to set the linearization
        // point. So far this is not supported by the matrix free solver
        // backends.
        //
        // _ls.setLinearizationPoint(*_x);

        _ls.apply(z,r,red);
        _linear_solver_result = _ls.result();
        timing = watch.elapsed();
        // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
        if (rank==0 && _verbose>=1)
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

      //! Return tolerance, i.e. target residual reduction
      Real reduction() const
      {
        return _reduction;
      }

      //! Set tolerance, i.e. target residual reduction
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
