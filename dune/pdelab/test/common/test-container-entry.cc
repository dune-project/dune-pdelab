#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/common/container_entry.hh>

#include <dune/typetree/treepath.hh>

#include <dune/common/tuplevector.hh>

#include <gtest/gtest.h>

#include <vector>
#include <random>

TEST(TestContainerEntry, Vector) {
  std::size_t size = 10'000;
  std::uniform_int_distribution unif(-100, 100);
  std::default_random_engine re;
  re.seed(0);

  std::vector<int> container(size);

  for (int& val : container)
    val = unif(re);

  using Dune::PDELab::containerEntry;
  using Dune::TypeTree::treePath;

  for (std::size_t i = 0; i != size; ++i)
    EXPECT_EQ(container[i], containerEntry(container, treePath(i)));
}

TEST(TestContainerEntry, BlockedVector) {
  std::size_t size = 1'000;
  std::uniform_int_distribution unif(-100, 100);
  std::default_random_engine re;
  re.seed(0);

  std::vector<std::vector<int>> container(size);

  for (auto& val0 : container) {
    val0.resize(size);
    for (int& val1 : val0)
      val1 = unif(re);
  }

  using Dune::PDELab::containerEntry;
  using Dune::TypeTree::treePath;

  for (std::size_t i = 0; i != size; ++i)
    for (std::size_t j = 0; j != size; ++j)
    EXPECT_EQ(container[i][j], containerEntry(container, treePath(i,j)));
}


TEST(TestContainerEntry, TupleVector) {
  using Dune::makeTupleVector;
  auto container = makeTupleVector(makeTupleVector(1, 2.), makeTupleVector(3., std::vector<int>{3, 4}));

  using Dune::PDELab::containerEntry;
  using Dune::TypeTree::treePath;

  using namespace Dune::Indices;
  EXPECT_EQ(container[_0][_0],    containerEntry(container, treePath(_0, _0)));
  EXPECT_NEAR(container[_0][_1],  containerEntry(container, treePath(_0, _1)), 1e-12);
  EXPECT_NEAR(container[_1][_0],  containerEntry(container, treePath(_1, _0)), 1e-12);
  EXPECT_EQ(container[_1][_1][0], containerEntry(container, treePath(_1, _1, 0)));
  EXPECT_EQ(container[_1][_1][1], containerEntry(container, treePath(_1, _1, 1)));
}

struct CustomContainer {
  friend constexpr int containerEntry(CustomContainer, Dune::PDELab::Concept::MultiIndex auto) {
    return 1;
  }
};

TEST(TestContainerEntry, CustomContainer) {
  using namespace Dune::Indices;
  CustomContainer container;

  using Dune::PDELab::containerEntry;
  using Dune::TypeTree::treePath;

  static_assert(1 == containerEntry(container, treePath()));
  static_assert(1 == containerEntry(container, treePath(0,1,2,3,4)));
  static_assert(1 == containerEntry(container, treePath(_0, _0)));
}
