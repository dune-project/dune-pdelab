#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "fixture.hh"
#include "ordering.hh"

// declare a dummy class that forwards it's template argument as the fixture to test
template <typename Fixture> class ForwardFixture : public Fixture {};

// declare the templated fixture
TYPED_TEST_SUITE_P(ForwardFixture);

// declare test to make on each space
TYPED_TEST_P(ForwardFixture, TestOrdering) {
  test_ordering(this->makeFixtureBasis(this->_grid->leafGridView()));
}

// register the test
REGISTER_TYPED_TEST_SUITE_P(ForwardFixture, TestOrdering);

// instantiate the test for each of the spaces in the BasisFixtures
INSTANTIATE_TYPED_TEST_SUITE_P(Ordering, ForwardFixture, BasisFixtures);

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  Dune::MPIHelper::instance(argc, argv);
  return RUN_ALL_TESTS();
}
