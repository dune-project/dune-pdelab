include(UseEigen)
include(UseOpenMP)

if(DUNE_PDELAB_ENABLE_TRACING)
  include(UsePerfetto)
endif()

# Trying to run a sequential UG in a parallel PDELab application
# will result in very subtle errors. We therefore issue a warning
# and set a preprocessor variable to give meaningful error messages.
if(MPI_FOUND AND UG_FOUND)
  if(NOT UG_PARALLEL)
    message(WARNING "You are using a sequential UG in a parallel PDELab environment!")
    set(PDELAB_SEQUENTIAL_UG 1)
  endif()
endif()
