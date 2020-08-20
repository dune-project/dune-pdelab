// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_NEWTONERRORS_HH
#define DUNE_PDELAB_SOLVER_NEWTONERRORS_HH

#ifndef DOXYGEN

namespace Dune::PDELab
{
  // Exception classes used in NewtonSolver
  class NewtonError : public Exception {};
  class NewtonDefectError : public NewtonError {};
  class NewtonLinearSolverError : public NewtonError {};
}

#endif // DOXYGEN

#endif
