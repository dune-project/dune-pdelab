#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/common/for_each.hh>

#include <dune/common/tuplevector.hh>

#include <gtest/gtest.h>

#include <atomic>
#include <vector>
#include <map>
#include <random>

std::uniform_int_distribution<std::size_t> unif(0, 100);
std::default_random_engine re;

TEST(TestForEach, Vector) {
  re.seed(100);

  std::size_t size = 10'000;
  std::vector<std::size_t> container(size);

  for (auto& val : container)
    val = unif(re);

  using Dune::PDELab::forEach;
  using std::as_const;

  std::size_t count = 0;
  forEach(as_const(container), [&container,&count](const auto& entry){
    EXPECT_EQ(container[count++], entry);
  });

  forEach(container, [](auto& entry){
    entry = unif(re);
  });

  std::size_t sum = 0;
  forEach(as_const(container), [&container,&sum](const auto& entry, int i){
    EXPECT_EQ(container[i], entry);
    sum += entry;
  });

  std::atomic<std::size_t> asum{0};
  forEach(std::execution::par_unseq, as_const(container), [&container,&asum](const auto& entry, int i){
    EXPECT_EQ(container[i], entry);
    asum.fetch_add(entry, std::memory_order::relaxed);
  });

  EXPECT_EQ(sum, asum.load());
}

TEST(TestForEach, Map) {
  re.seed(200);

  std::size_t size = 100;
  std::map<std::size_t, std::size_t> container;

  for (std::size_t i = 0; i != size; ++i)
    container[unif(re)] = unif(re);

  using Dune::PDELab::forEach;
  using std::as_const;

  std::size_t sum = 0;
  forEach(as_const(container), [&container,&sum](const auto& entry){
    EXPECT_EQ(container[entry.first], entry.second);
    sum += entry.second;
  });

  std::atomic<std::size_t> asum{0};
  forEach(std::execution::par_unseq, as_const(container), [&container,&asum](const auto& entry){
    EXPECT_EQ(container[entry.first], entry.second);
    asum.fetch_add(entry.second, std::memory_order::relaxed);
  });

  EXPECT_EQ(sum, asum.load());
}

TEST(TestForEach, TupleVector) {
  using Dune::makeTupleVector;
  auto container = makeTupleVector(makeTupleVector(-1, 0u, 1.), std::vector{3.,4.,5.});

  using Dune::PDELab::forEach;

  double sum = 0.;
  forEach(as_const(container), [&sum](const auto& entry1){
    forEach(entry1, [&sum](const auto& entry2){
      sum += entry2;
    });
  });

  EXPECT_NEAR(sum,  12, 1e-12);

  std::atomic<std::size_t> asum{0};
  forEach(std::execution::par_unseq, as_const(container), [&asum](const auto& entry1){
    forEach(std::execution::par_unseq, entry1, [&asum](const auto& entry2){
      asum.fetch_add(entry2, std::memory_order::relaxed);
    });
  });

  EXPECT_NEAR(asum.load(),  12, 1e-12);
}
