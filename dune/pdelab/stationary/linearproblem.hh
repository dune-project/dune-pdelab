#ifndef DUNE_PDELAB_STATIONARYLINEARPROBLEM_HH
#define DUNE_PDELAB_STATIONARYLINEARPROBLEM_HH

#include <iostream>

#include <dune/common/timer.hh>
#include <dune/common/deprecated.hh>
#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/backendselector.hh>
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
      double linear_solver_time; // Cumulative time for linear sovler
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
      typedef typename V::ElementType Real;
      typedef typename GO::Traits::Jacobian M;
      typedef typename GO::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
      typedef typename Dune::PDELab::BackendVectorSelector<TrialGridFunctionSpace,Real>::Type W;
      typedef GO GridOperator;

    public:
      typedef StationaryLinearProblemSolverResult<double> Result;

      StationaryLinearProblemSolver(const GO& go, V& x, LS& ls, typename V::ElementType reduction, typename V::ElementType min_defect = 1e-99, int verbose=1) DUNE_DEPRECATED_MSG("Use StationaryLinearProblemSolver(const GO&, LS&, V&, ...) instead.")
        : _go(go)
        , _ls(ls)
        , _x(&x)
        , _reduction(reduction)
        , _min_defect(min_defect)
        , _hanging_node_modifications(true)
        , _keep_matrix(true)
        , _verbose(verbose)
      {}

      StationaryLinearProblemSolver(const GO& go, LS& ls, V& x, typename V::ElementType reduction, typename V::ElementType min_defect = 1e-99, int verbose=1)
        : _go(go)
        , _ls(ls)
        , _x(&x)
        , _reduction(reduction)
        , _min_defect(min_defect)
        , _hanging_node_modifications(true)
        , _keep_matrix(true)
        , _verbose(verbose)
      {}

      StationaryLinearProblemSolver (const GO& go, LS& ls, typename V::ElementType reduction, typename V::ElementType min_defect = 1e-99, int verbose=1)
        : _go(go)
        , _ls(ls)
        , _x(nullptr)
        , _reduction(reduction)
        , _min_defect(min_defect)
        , _hanging_node_modifications(true)
        , _keep_matrix(true)
        , _verbose(verbose)
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

      void apply(V& x) {
        _x = &x;
        apply();
      }

      void apply ()
      {
        Dune::Timer watch;
        double timing,assembler_time=0;

        // assemble matrix; optional: assemble only on demand!
        watch.reset();

        if (!_jacobian)
          {
            _jacobian = make_shared<M>(_go);
            timing = watch.elapsed();
            if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
              std::cout << "=== matrix setup (max) " << timing << " s" << std::endl;
            watch.reset();
            assembler_time += timing;
          }
        else if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
          std::cout << "=== matrix setup skipped (matrix already allocated)" << std::endl;

        (*_jacobian) = 0.0;
        if (_hanging_node_modifications)
          {
            Dune::PDELab::set_shifted_dofs(_go.localAssembler().trialConstraints(),0.0,*_x); // set hanging node DOFs to zero
            _go.localAssembler().backtransform(*_x); // interpolate hanging nodes adjacent to Dirichlet nodes
          }
        _go.jacobian(*_x,*_jacobian);

        timing = watch.elapsed();
        // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
        if (_go.trialGridFunctionSpace().gridView().comm().rank()==0 && _verbose>=1)
          std::cout << "=== matrix assembly (max) " << timing << " s" << std::endl;
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

        typename V::ElementType defect = _ls.norm(r);

        // compute correction
        watch.reset();
        V z(_go.trialGridFunctionSpace(),0.0);
        typename V::ElementType red = std::min(_reduction,defect/_min_defect);
        if (_go.trialGridFunctionSpace().gridView().comm().rank()==0)
          std::cout << "=== solving (reduction: " << red << ") ";
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
        _res.first_defect = defect;
        _res.defect = defect*_reduction;
        _res.linear_solver_iterations = _linear_solver_result.iterations;

        // and update
        if (_hanging_node_modifications)
          Dune::PDELab::set_shifted_dofs(_go.localAssembler().trialConstraints(),0.0,*_x); // set hanging node DOFs to zero
        *_x -= z;
        if (_hanging_node_modifications)
          _go.localAssembler().backtransform(*_x); // interpolate hanging nodes adjacent to Dirichlet nodes

        if (!_keep_matrix)
          _jacobian = nullptr;
      }

      //! Discard the stored Jacobian matrix.
      void discardMatrix()
      {
        if(_jacobian)
          _jacobian = nullptr;
      }

      const Dune::PDELab::LinearSolverResult<double>& ls_result() const{
        return _linear_solver_result;
      }

    private:
      const GO& _go;
      LS& _ls;
      V* _x;
      shared_ptr<M> _jacobian;
      typename V::ElementType _reduction;
      typename V::ElementType _min_defect;
      Dune::PDELab::LinearSolverResult<double> _linear_solver_result;
      Result _res;
      bool _hanging_node_modifications;
      bool _keep_matrix;
      int _verbose;
    };

  } // namespace PDELab
} // namespace Dune

#endif
