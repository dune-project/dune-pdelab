dune_add_test(
  NAME test-operator-forward-runge-kutta-scalar
  SOURCES test-runge-kutta-scalar.cc
  LINK_LIBRARIES ${GTEST_MAIN_LIBRARIES})

dune_add_test(
  NAME test-operator-forward-runge-kutta-pde
  SOURCES test-runge-kutta-pde.cc)
