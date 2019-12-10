#ifndef DUNE_PDELAB_SOLVER_UTILITY_HH
#define DUNE_PDELAB_SOLVER_UTILITY_HH

namespace Dune {
  namespace PDELab {

    template <typename RFType>
    struct PDESolverResult : LinearSolverResult<RFType>
    {
      RFType first_defect;       // the first defect
      RFType defect;             // the final defect
      double assembler_time;     // Cumulative time for matrix assembly
      double linear_solver_time; // Cumulative time for linear solver
      int linear_solver_iterations; // Total number of linear iterations

      PDESolverResult()
      {
        clear();
      }

      void clear()
      {
        LinearSolverResult<RFType>::clear();
        first_defect = 0.0;
        defect = 0.0;
        assembler_time = 0.0;
        linear_solver_time = 0.0;
        linear_solver_iterations = 0;
      }
    };

  } // namespace PDELab
} // namespace Dune

#endif
