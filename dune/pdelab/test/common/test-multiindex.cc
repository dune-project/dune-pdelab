#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/concepts/multiindex.hh>
#include <dune/pdelab/common/multiindex.hh>

#include <dune/typetree/treepath.hh>

#include <unordered_map>

#include <gtest/gtest.h>

TEST(TestMultiIndex, TestReservedMultiIndex) {
  using namespace Dune::PDELab;
  using namespace Dune::Indices;

  Concept::MultiIndex auto mi = MultiIndex<std::size_t,5>();
  mi.resize(4); mi[0] = _1; mi[1] = 3; mi[2] = _2; mi[3] = 5;

  EXPECT_EQ(mi[_0], 1);
  EXPECT_EQ(mi[3], 5);

  EXPECT_EQ(back(mi), 5);
  EXPECT_EQ(back(push_back(mi, 3)), 3);
  EXPECT_EQ(back(pop_back(mi)), 2);
  EXPECT_EQ(back(pop_back(pop_back(mi))), 3);

  EXPECT_EQ(front(push_front(mi,0)), 0);
  EXPECT_EQ(front(pop_front(mi)), 3);
  EXPECT_EQ(front(pop_front(pop_front(mi))), 2);

  using MI = std::decay_t<decltype(mi)>;
  std::unordered_map<MI, int> map;

  map[mi] = 2;

  Concept::MultiIndex auto mi0 = mi;
  mi0.resize(2); mi0[0] = _2; mi0[1] = 4;
  map[mi0] = 3;

  EXPECT_EQ(map[mi], 2);
  EXPECT_EQ(map[mi0], 3);
}

TEST(TestMultiIndex, TestHybridMultiIndex) {
  using Dune::TypeTree::treePath;
  using namespace Dune::PDELab;
  using namespace Dune::Indices;

  constexpr Concept::MultiIndex auto mi = treePath(_1,3,_2,5);

  EXPECT_EQ(mi[_0], 1);
  EXPECT_EQ(mi[3], 5);

  EXPECT_EQ(back(mi), 5);
  static_assert(back(push_back(mi, _3)) == 3);
  EXPECT_EQ(back(push_back(mi, 3)), 3);
  static_assert(back(pop_back(mi)) == 2);
  EXPECT_EQ(back(pop_back(mi)), 2);
  EXPECT_EQ(back(pop_back(pop_back(mi))), 3);

  static_assert(front(mi) == 1);
  static_assert(front(push_front(mi,_0)) == 0);
  EXPECT_EQ(front(push_front(mi,0)), 0);
  EXPECT_EQ(front(pop_front(mi)), 3);
  static_assert(front(pop_front(pop_front(mi))) == 2);
  EXPECT_EQ(front(pop_front(pop_front(mi))), 2);

  // using MI = std::decay_t<decltype(mi)>;
  // std::unordered_map<MI, int> map;

  // map[mi] = 2;
  // map[treePath(_1,4,_2,6)] = 3;

  // EXPECT_EQ(map[mi], 2);
  // EXPECT_EQ(map[treePath(_1,4,_2,6)], 3);
}
