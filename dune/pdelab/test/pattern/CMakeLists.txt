dune_add_test(
  SOURCES test-sparsity-pattern.cc
  LINK_LIBRARIES GTest::gtest_main)

dune_add_test(
  NAME test-pattern-basis-to-matrix
  SOURCES test-basis-to-matrix.cc)
