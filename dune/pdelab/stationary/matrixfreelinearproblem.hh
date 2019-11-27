#ifndef DUNE_PDELAB_MATRIXFREESTATIONARYLINEARPROBLEM_HH
#define DUNE_PDELAB_MATRIXFREESTATIONARYLINEARPROBLEM_HH

#include <iostream>

#include <dune/common/timer.hh>
#include <dune/common/deprecated.hh>
#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/backend/solver.hh>
#include <dune/pdelab/stationary/linearproblembase.hh>

namespace Dune {
  namespace PDELab {

    //===============================================================
    // A class for solving linear stationary problems with matrix-free.
    // methods.
    // It computes the right hand side and solves the problem.
    //===============================================================

    /** \brief Matrix free iterative solver for linear problems in residual
     *         form.
     *
     * Solves problem of the from \f$r(u,v)=0\f$ where \f$r\f$ is a
     * bilinear form.
     */
    template<typename GO, typename LS, typename V>
    class MatrixFreeStationaryLinearProblemSolver
    {
      /** \brief Field value type */
      typedef typename Dune::template FieldTraits<typename V::ElementType >::real_type Real;
      /** \brief Trial grid function space */
      typedef typename GO::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
      /** \brief Vector backend for trial grid function space */
      using W = Dune::PDELab::Backend::Vector<TrialGridFunctionSpace,typename V::ElementType>;
      /** \brief Grid operator to solve for */
      typedef GO GridOperator;

    public:
      /** \brief Class holding results */
      typedef StationaryLinearProblemSolverResult<double> Result;

      /** \brief Construct new solver instance, based on initial guess
       *         of solution
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
        , _verbose(verbose)
        , _hanging_node_modifications(false)
      {}

      /** \brief Set flag that controls hanging node transformations
       *
       * Set whether the solver should apply the necessary
       * transformations for calculations on hanging nodes.
       *
       * \param[in] b Apply hanging node modifications?
       */
      void setHangingNodeModifications(bool b)
      {
        _hanging_node_modifications=b;
      }

      /** \brief Return flag that control hanging node transformations
       *
       * Return whether the solver performs the necessary transformations
       * for calculations on hanging nodes.
       */
      bool hangingNodeModifications() const
      {
        return _hanging_node_modifications;
      }

      /** \brief Return result object */
      const Result& result() const
      {
        return _res;
      }

      void apply(V& x)
      {
        _x = &x;
        apply();
      }

      /** \brief Solve linear problem */
      void apply ()
      {
        Dune::Timer watch;
        double timing,assembler_time=0;
        int rank=_go.trialGridFunctionSpace().gridView().comm().rank();
        // assemble matrix; optional: assemble only on demand!
        watch.reset();
        if (rank==0 && _verbose>=1) {
          std::cout << "=== matrix setup NOT REQUIRED " << std::endl;
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
        if (rank==0)
          std::cout << "=== solving (reduction: " << red << ") " << std::endl;

        // TODO: For nonlinear problems we need to set the linearization
        // point. So far this is not supported by the matrix free solver
        // backends.
        //
        // _ls.setLinearizationPoint(*_x);

        _ls.apply(z,r,red); // solver makes right hand side consistent
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

      /** \brief Return linear solver result object */
      const Dune::PDELab::LinearSolverResult<double>& ls_result() const
      {
        return _linear_solver_result;
      }

      /** \brief Return tolerance, i.e. target residual reduction */
      Real reduction() const
      {
        return _reduction;
      }

      /** \brief Set tolerance, i.e. target residual reduction
       *
       * \param[in] reduction New target residual reduction
       */
      void setReduction(Real reduction)
      {
        _reduction = reduction;
      }


    private:
      /** \brief Grid operator representing linear problem */
      const GO& _go;
      /** \brief Linear solver backend */
      LS& _ls;
      /** \brief Solution grid function */
      V* _x;
      /** \brief Tolerance, i.e. target residual reduction */
      Real _reduction;
      /** \brief Minimal absolute residual
       *
       * abort solve once this tolerance is achieved
       */
      Real _min_defect;
      /** \brief Linear solver result object */
      Dune::PDELab::LinearSolverResult<double> _linear_solver_result;
      /** \brief Achieved residual */
      Result _res;
      /** \brief Verbosity level */
      int _verbose;
      /** \brief Flag controlling hanging node modifications */
      const bool _hanging_node_modifications;
    };

  } // namespace PDELab
} // namespace Dune

#endif
