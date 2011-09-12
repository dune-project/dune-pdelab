// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PETSCUTILITY_HH
#define DUNE_PETSCUTILITY_HH

#if HAVE_PETSC

#include <petscsys.h>
#include <dune/common/exceptions.hh>

namespace Dune {
  namespace PDELab {

    class PetscException
      : public Dune::Exception
    {};

  } // namespace PDELab
} // namespace Dune


#define PETSC_GUARD_START {

#define PETSC_GUARD_END(err_var) \
  if (err_var != 0) DUNE_THROW(Dune::PDELab::PetscException,"Petsc problem"); \
}

#define PETSC_CALL(x) \
PETSC_GUARD_START \
PetscErrorCode __petsc_err = (x); \
PETSC_GUARD_END(__petsc_err)

#endif // HAVE_PETSC

#endif // DUNE_PETSCUTILITY_HH
