// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_NEWTON_NEWTONBASE_HH
#define DUNE_PDELAB_NEWTON_NEWTONBASE_HH

#include <dune/common/exceptions.hh>
#include <dune/pdelab/backend/solver.hh>

namespace Dune {
  namespace PDELab {

    // Exception classes used in NewtonSolver
    class NewtonError : public Exception {};
    class NewtonDefectError : public NewtonError {};
    class NewtonLinearSolverError : public NewtonError {};
    class NewtonLineSearchError : public NewtonError {};
    class NewtonNotConverged : public NewtonError {};

    // Status information of Newton's method
    template<class RFType>
    struct NewtonResult : LinearSolverResult<RFType>
    {
      RFType first_defect;       // the first defect
      RFType defect;             // the final defect
      double assembler_time;     // Cumulative time for matrix assembly
      double linear_solver_time; // Cumulative time for linear sovler
      int linear_solver_iterations; // Total number of linear iterations

      NewtonResult() :
        first_defect(0.0), defect(0.0), assembler_time(0.0), linear_solver_time(0.0),
        linear_solver_iterations(0) {}
    };

    enum struct LineSearchStrategy {

      /** \brief don't do any linesearch or damping */
      noLineSearch,

      /** \brief perform a linear search for the optimal damping parameter with multiples of damping

          the strategy was described in <a href="http://dx.doi.org/10.1007/BF01406516">[Hackbusch and Reusken, 1989]</a> */
     hackbuschReusken,

      /** \brief same as hackbuschReusken, but doesn't fail if the best update is still not good enough */
      hackbuschReuskenAcceptBest
      };

  } // end namespace PDELab
} // end namespace Dune
#endif
