
dune_add_test(SOURCES recipe-geometry-grid.cc)
dune_add_test(SOURCES recipe-integrating-grid-functions.cc)
dune_add_test(SOURCES recipe-grid-function-operations.cc)
dune_add_test(SOURCES recipe-blocking.cc)
dune_add_test(SOURCES recipe-linear-system-assembly.cc)
dune_add_test(SOURCES recipe-linear-system-solution-istl.cc)
dune_add_test(SOURCES recipe-linear-system-solution-pdelab.cc)
dune_add_test(SOURCES recipe-communication.cc)
dune_add_test(SOURCES recipe-operator-splitting.cc
              MPI_RANKS 1 2 4
              TIMEOUT 30)
