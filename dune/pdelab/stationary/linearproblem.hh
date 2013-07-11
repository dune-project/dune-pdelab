#ifndef DUNE_PDELAB_STATIONARYLINEARPROBLEM_HH
#define DUNE_PDELAB_STATIONARYLINEARPROBLEM_HH

#include<dune/common/timer.hh>
#include<dune/pdelab/backend/backendselector.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<iostream>

#include "../backend/solver.hh"

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

      StationaryLinearProblemSolverResult() :
        first_defect(0.0), defect(0.0), assembler_time(0.0), linear_solver_time(0.0),
        linear_solver_iterations(0) {}
    };

    template<class GOS, class LS, class V> 
    class StationaryLinearProblemSolver
    {
      typedef typename V::ElementType Real;
      typedef typename GOS::Traits::Jacobian M;
      typedef typename GOS::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
      typedef typename Dune::PDELab::BackendVectorSelector<TrialGridFunctionSpace,Real>::Type W;
      
    public:
      typedef StationaryLinearProblemSolverResult<double> Result;

      StationaryLinearProblemSolver (const GOS& gos_, V& x_, LS& ls_, typename V::ElementType reduction_, typename V::ElementType mindefect_ = 1e-99, int verb_=1)
        : gos(gos_), ls(ls_), x(&x_), reduction(reduction_), mindefect(mindefect_), hangingNodeModifications(true), verbose(verb_)
      {
      }

      StationaryLinearProblemSolver (const GOS& gos_, LS& ls_, V& x_, typename V::ElementType reduction_, typename V::ElementType mindefect_ = 1e-99, int verb_=1)
        : gos(gos_), ls(ls_), x(&x_), reduction(reduction_), mindefect(mindefect_), hangingNodeModifications(true), verbose(verb_)
      {
      }

      StationaryLinearProblemSolver (const GOS& gos_, LS& ls_, typename V::ElementType reduction_, typename V::ElementType mindefect_ = 1e-99, int verb_=1)
          : gos(gos_), ls(ls_), x(0), reduction(reduction_), mindefect(mindefect_), hangingNodeModifications(true), verbose(verb_)
      {
      }

      void setHangingNodeModifications (bool b)
      {
        hangingNodeModifications=b;
      }

      const Result& result() const
      {
        return res;
      }

      void apply (V& x_) {
        x = &x_;
        apply();
      }

      void apply ()
      {
        Dune::Timer watch;
        double timing,assembler_time=0;

        // assemble matrix; optional: assemble only on demand!
        watch.reset();

        M m(gos); 

        timing = watch.elapsed();
        // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
        if (gos.trialGridFunctionSpace().gridView().comm().rank()==0 && verbose>=1)
          std::cout << "=== matrix setup (max) " << timing << " s" << std::endl;
        watch.reset();
        assembler_time += timing;

        m = 0.0;
        if (hangingNodeModifications)
          {
            Dune::PDELab::set_shifted_dofs(gos.localAssembler().trialConstraints(),0.0,*x); // set hanging node DOFs to zero
            gos.localAssembler().backtransform(*x); // interpolate hanging nodes adjacent to Dirichlet nodes
          }
        gos.jacobian(*x,m);

        timing = watch.elapsed();
        // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
        if (gos.trialGridFunctionSpace().gridView().comm().rank()==0 && verbose>=1)
          std::cout << "=== matrix assembly (max) " << timing << " s" << std::endl;
        assembler_time += timing;

        // assemble residual
        watch.reset();

        W r(gos.testGridFunctionSpace(),0.0);
        gos.residual(*x,r);  // residual is additive

        timing = watch.elapsed();
        // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
        if (gos.trialGridFunctionSpace().gridView().comm().rank()==0 && verbose>=1)
          std::cout << "=== residual assembly (max) " << timing << " s" << std::endl;
        assembler_time += timing;
        res.assembler_time = assembler_time;

        typename V::ElementType defect = ls.norm(r);

        // compute correction
        watch.reset();
        V z(gos.trialGridFunctionSpace(),0.0);
        typename V::ElementType red = std::min(reduction,defect/mindefect);
        if (gos.trialGridFunctionSpace().gridView().comm().rank()==0)
          std::cout << "=== solving (reduction: " << red << ") ";
        ls.apply(m,z,r,red); // solver makes right hand side consistent
        linearsolverresult = ls.result();
        timing = watch.elapsed();
        // timing = gos.trialGridFunctionSpace().gridView().comm().max(timing);
        if (gos.trialGridFunctionSpace().gridView().comm().rank()==0 && verbose>=1)
          std::cout << timing << " s" << std::endl;
        res.linear_solver_time = timing;

        res.converged = linearsolverresult.converged;
        res.iterations = linearsolverresult.iterations;
        res.elapsed = linearsolverresult.elapsed;
        res.reduction = linearsolverresult.reduction;
        res.conv_rate = linearsolverresult.conv_rate;
        res.first_defect = defect;
        res.defect = defect*reduction;
        res.linear_solver_iterations = linearsolverresult.iterations;

        // and update
        if (hangingNodeModifications)
          Dune::PDELab::set_shifted_dofs(gos.localAssembler().trialConstraints(),0.0,*x); // set hanging node DOFs to zero
        *x -= z;
        if (hangingNodeModifications)
          gos.localAssembler().backtransform(*x); // interpolate hanging nodes adjacent to Dirichlet nodes
      }

      const Dune::PDELab::LinearSolverResult<double>& ls_result() const{
        return linearsolverresult;
      }

    private:
      const GOS& gos;
      LS& ls;
      V* x;
      typename V::ElementType reduction;
      typename V::ElementType mindefect;
      Dune::PDELab::LinearSolverResult<double> linearsolverresult;
      Result res;
      bool hangingNodeModifications;
      int verbose;
    };

  } // namespace PDELab
} // namespace Dune

#endif
