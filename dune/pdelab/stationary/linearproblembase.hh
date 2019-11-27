// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_STATIONARY_LINEARPROBLEMBASE_HH
#define DUNE_PDELAB_STATIONARY_LINEARPROBLEMBASE_HH

#include <dune/pdelab/backend/solver.hh>

namespace Dune {
  namespace PDELab {

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
  } // end namespace PDELab
} // end namespace Dune
#endif
