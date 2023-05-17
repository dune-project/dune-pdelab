dune_add_test(
  NAME test-operator-forward-runge-kutta-scalar
  SOURCES test-runge-kutta-scalar.cc
  LINK_LIBRARIES GTest::gtest_main)

dune_add_test(
  NAME test-operator-forward-runge-kutta-pde
  SOURCES test-runge-kutta-pde.cc)
