#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/pattern/sparsity_pattern.hh>
#include <dune/pdelab/common/multiindex.hh>

#include <dune/pdelab/concepts/size_provider.hh>

#include <dune/typetree/treepath.hh>

#include <gtest/gtest.h>

#include <random>

// ************ Leaf cases ********************

struct FlatSizeProvider {
  using size_type = std::size_t;
  using SizePrefix = Dune::PDELab::MultiIndex<size_type, 1>;

  size_type size(SizePrefix prefix) const {
    return prefix.empty() ? _size : 0;
  }

  size_type _size;
};

static_assert(Dune::PDELab::Concept::SizeProvider<FlatSizeProvider>);

TEST(TestPatterns, TestLowerTriangularLeafPattern) {

  FlatSizeProvider sp{5};
  Dune::PDELab::LeafSparsePattern<FlatSizeProvider,FlatSizeProvider> pattern{sp, {}, sp, {}, 2};

  using Dune::TypeTree::treePath;
  for (std::size_t row = 0; row != sp.size({}); ++row)
    for (std::size_t col = 0; col != row; ++col)
      pattern.addLink(treePath(row), treePath(col));

  pattern.sort();

  auto row_nnzs = pattern.rowNonZeros();
  EXPECT_EQ(row_nnzs, (std::vector<std::size_t>{0,1,2,3,4}) );
  for (std::size_t row = 0; row != row_nnzs.size(); ++row) {
    std::size_t cols = 0;
    for(auto it = pattern.columnBegin(row); it != pattern.columnEnd(row); ++it) {
      EXPECT_EQ(cols, *it);
      ++cols;
    }
    EXPECT_EQ(cols, row_nnzs[row]);
  }

  EXPECT_FALSE(pattern < pattern);
}

TEST(TestPatterns, TestReversedLowerTriangularLeafPattern) {

  FlatSizeProvider sp{5};
  Dune::PDELab::LeafSparsePattern<FlatSizeProvider,FlatSizeProvider> pattern{sp, {}, sp, {}, 2};

  using Dune::TypeTree::treePath;
  for (std::size_t row = 0; row != sp.size({}); ++row)
    for (std::size_t col = 0; col != row; ++col)
      pattern.addLink(treePath(row), treePath(row - col - 1));

  pattern.sort();

  auto row_nnzs = pattern.rowNonZeros();
  EXPECT_EQ(row_nnzs, (std::vector<std::size_t>{0,1,2,3,4}) );
  for (std::size_t row = 0; row != row_nnzs.size(); ++row) {
    std::size_t cols = 0;
    for(auto it = pattern.columnBegin(row); it != pattern.columnEnd(row); ++it) {
      EXPECT_EQ(cols, *it);
      ++cols;
    }
    EXPECT_EQ(cols, row_nnzs[row]);
  }

  EXPECT_FALSE(pattern < pattern);
}

TEST(TestPatterns, TestUpperTriangularLeafPattern) {

  FlatSizeProvider sp{5};
  Dune::PDELab::LeafSparsePattern<FlatSizeProvider,FlatSizeProvider> pattern{sp, {}, sp, {}, 2};

  using Dune::TypeTree::treePath;
  for (std::size_t row = 0; row != sp.size({}); ++row)
    for (std::size_t col = row; col != sp.size({}); ++col)
      pattern.addLink(treePath(row), treePath(col));

  pattern.sort();

  auto row_nnzs = pattern.rowNonZeros();
  EXPECT_EQ(row_nnzs, (std::vector<std::size_t>{5,4,3,2,1}) );
  for (std::size_t row = 0; row != row_nnzs.size(); ++row) {
    std::size_t cols = 0;
    for(auto it = pattern.columnBegin(row); it != pattern.columnEnd(row); ++it) {
      EXPECT_EQ(cols + row, *it);
      ++cols;
    }
    EXPECT_EQ(cols, row_nnzs[row]);
  }

  EXPECT_FALSE(pattern < pattern);
}

TEST(TestPatterns, TestRandomLeafPattern) {

// test different fill and buffer rations to go to all code paths of the pattern
#ifndef NDEBUG
std::size_t max_size_shift = 6;
#else
std::size_t max_size_shift = 16;
#endif

  for (std::size_t size = 2<<0; size != 2u<<max_size_shift; size <<= 1) {
    for (auto fill_ratio : {.25, .5, .75, 1.}) {
      for (auto buffer_ratio : { 0.5, 1.0, 1.5}) {
        std::size_t max_entries = static_cast<std::size_t>(fill_ratio*std::min<std::size_t>(size,10));
        std::size_t entries_per_row = static_cast<std::size_t>(buffer_ratio*max_entries);

        FlatSizeProvider sp{size};
        Dune::PDELab::LeafSparsePattern<FlatSizeProvider,FlatSizeProvider> pattern{sp, {}, sp, {}, entries_per_row};


        std::uniform_int_distribution<std::size_t> unif(0, sp.size({})-1);
        std::default_random_engine re;
        re.seed(0);
        std::vector<std::set<std::size_t>> pattern_set(sp.size({}));

        auto start_t = std::chrono::system_clock::now();

        using Dune::TypeTree::treePath;
        for (std::size_t row = 0; row != sp.size({}); ++row)
          for (std::size_t i = 0; i != max_entries; ++i) {
            std::size_t col = unif(re);
            pattern.addLink(treePath(row), treePath(col));
            pattern_set[row].insert(col);
          }

        std::chrono::duration<double, std::micro> duration = std::chrono::system_clock::now() - start_t;
        // std::cout << size << " \t" << max_entries << " \t" << entries_per_row << " \t" << duration.count() << std::endl;

        pattern.sort();

        auto row_nnzs = pattern.rowNonZeros();
        for (std::size_t row = 0; row != row_nnzs.size(); ++row) {
          EXPECT_EQ(row_nnzs[row], pattern_set[row].size());
          std::size_t cols = 0;
          for(auto it = pattern.columnBegin(row); it != pattern.columnEnd(row); ++it) {
            EXPECT_TRUE(pattern_set[row].contains(*it)) << "Expected value " << *it;
            pattern_set[row].erase(*it);
            ++cols;
          }
          EXPECT_TRUE(pattern_set[row].empty());
          EXPECT_EQ(cols, row_nnzs[row]);
        }

        EXPECT_FALSE(pattern < pattern);
      }
    }
  }
}

TEST(TestPatterns, TestCompareLeafPattern) {

  FlatSizeProvider sp{5};

  using Dune::TypeTree::treePath;
  for (std::size_t i = 0; i != sp.size({})-1; ++i) {
    Dune::PDELab::LeafSparsePattern<FlatSizeProvider,FlatSizeProvider> pattern_a{sp, {}, sp, {}, 2};
    Dune::PDELab::LeafSparsePattern<FlatSizeProvider,FlatSizeProvider> pattern_b{sp, {}, sp, {}, 2};

    for (std::size_t j = 0; j != i; ++j) {
      pattern_a.addLink(treePath(0), treePath(j));
      pattern_b.addLink(treePath(0), treePath(j));
    }

    pattern_a.sort();

    for (std::size_t j = i; j != sp.size({}); ++j) {
      pattern_b.addLink(treePath(0), treePath(j));
      pattern_b.sort();

      EXPECT_LT(pattern_a.rowNonZeros()[0], pattern_b.rowNonZeros()[0]);

      // check that there is no side effects
      EXPECT_FALSE(pattern_a < pattern_a);
      EXPECT_FALSE(pattern_b < pattern_b);

      // b must be bigger than a since it has more entries (and everything else is equal)
      EXPECT_TRUE(pattern_a < pattern_b);
    }
  }
}


// // ************ Blocked cases ********************

struct BlockedSizeProvider {
  using size_type = std::size_t;
  using SizePrefix = Dune::PDELab::MultiIndex<size_type, 2>;

  size_type size(SizePrefix prefix) const {
    return prefix.size() < 2 ? _size : 0;
  }

  size_type _size;
};

static_assert(Dune::PDELab::Concept::SizeProvider<BlockedSizeProvider>);

TEST(TestPatterns, TestLowerTriangularBlockedPattern) {

  BlockedSizeProvider sp{5};
  using LeafPattern = Dune::PDELab::LeafSparsePattern<BlockedSizeProvider,BlockedSizeProvider>;
  Dune::PDELab::BlockedSparsePattern<LeafPattern> pattern{sp, {}, sp, {}, [](auto,auto){return 2;}};

  using Dune::TypeTree::treePath;
  for (std::size_t row = 0; row != sp.size({}); ++row)
    for (std::size_t col = 0; col != row; ++col)
      pattern.addLink(treePath(row,0), treePath(col,0));

  pattern.sort();

  auto row_nnzs = pattern.rowNonZeros();
  EXPECT_EQ(row_nnzs, (std::vector<std::size_t>{0,1,2,3,4}) );
  for (std::size_t row = 0; row != row_nnzs.size(); ++row) {
    std::size_t cols = 0;
    for(auto it = pattern.columnBegin(row); it != pattern.columnEnd(row); ++it) {
      EXPECT_EQ(cols, *it);
      ++cols;
    }
    EXPECT_EQ(cols, row_nnzs[row]);
  }

  EXPECT_FALSE(pattern < pattern);
}


TEST(TestPatterns, TestReversedLowerTriangularBlockedPattern) {

  BlockedSizeProvider sp{5};
  using LeafPattern = Dune::PDELab::LeafSparsePattern<BlockedSizeProvider,BlockedSizeProvider>;
  Dune::PDELab::BlockedSparsePattern<LeafPattern> pattern{sp, {}, sp, {}, [](auto,auto){return 2;}};

  using Dune::TypeTree::treePath;
  for (std::size_t row = 0; row != sp.size({}); ++row)
    for (std::size_t col = 0; col != row; ++col)
      pattern.addLink(treePath(row,0), treePath(row - col - 1,0));

  pattern.sort();

  auto row_nnzs = pattern.rowNonZeros();
  EXPECT_EQ(row_nnzs, (std::vector<std::size_t>{0,1,2,3,4}) );
  for (std::size_t row = 0; row != row_nnzs.size(); ++row) {
    std::size_t cols = 0;
    for(auto it = pattern.columnBegin(row); it != pattern.columnEnd(row); ++it) {
      EXPECT_EQ(cols, *it);
      ++cols;
    }
    EXPECT_EQ(cols, row_nnzs[row]);
  }

  EXPECT_FALSE(pattern < pattern);
}

TEST(TestPatterns, TestUpperTriangularBlockedPattern) {

  BlockedSizeProvider sp{5};
  using LeafPattern = Dune::PDELab::LeafSparsePattern<BlockedSizeProvider,BlockedSizeProvider>;
  Dune::PDELab::BlockedSparsePattern<LeafPattern> pattern{sp, {}, sp, {}, [](auto,auto){return 2;}};

  using Dune::TypeTree::treePath;
  for (std::size_t row = 0; row != sp.size({}); ++row)
    for (std::size_t col = row; col != sp.size({}); ++col)
      pattern.addLink(treePath(row,0), treePath(col,0));

  pattern.sort();

  auto row_nnzs = pattern.rowNonZeros();
  EXPECT_EQ(row_nnzs, (std::vector<std::size_t>{5,4,3,2,1}) );
  for (std::size_t row = 0; row != row_nnzs.size(); ++row) {
    std::size_t cols = 0;
    for(auto it = pattern.columnBegin(row); it != pattern.columnEnd(row); ++it) {
      EXPECT_EQ(cols + row, *it);
      ++cols;
    }
    EXPECT_EQ(cols, row_nnzs[row]);
  }

  EXPECT_FALSE(pattern < pattern);
}



TEST(TestPatterns, TestRandomBlockedPattern) {
  // in this case we only check the outer block
  // that is, multi indices for the inner block, are always (0,0)

// test different fill and buffer rations to go to all code paths of the pattern
#ifndef NDEBUG
std::size_t max_size_shift = 6;
#else
std::size_t max_size_shift = 11;
#endif

  for (std::size_t size = 2<<0; size != 2u<<max_size_shift; size <<= 1) {
    for (auto fill_ratio : {.25, .5, .75, 1.}) {
      for (auto buffer_ratio : {0.5, 1.0, 1.5}) {
        std::size_t max_entries = static_cast<std::size_t>(fill_ratio*std::min<std::size_t>(size,10));
        std::size_t entries_per_row = static_cast<std::size_t>(buffer_ratio*max_entries);

        BlockedSizeProvider sp{size};
        using LeafPattern = Dune::PDELab::LeafSparsePattern<BlockedSizeProvider,BlockedSizeProvider>;

        Dune::PDELab::BlockedSparsePattern<LeafPattern> pattern{sp, {}, sp, {},
        [=](auto row_prefx,auto col_prefix) -> std::size_t {
          return row_prefx.empty() ? entries_per_row : 1;
        }};


        std::uniform_int_distribution<std::size_t> unif(0, sp.size({})-1);
        std::default_random_engine re;
        re.seed(0);
        std::vector<std::set<std::size_t>> pattern_set(sp.size({}));

        auto start_t = std::chrono::system_clock::now();

        using Dune::TypeTree::treePath;
        for (std::size_t row = 0; row != sp.size({}); ++row)
          for (std::size_t i = 0; i != max_entries; ++i) {
            std::size_t col = unif(re);
            pattern.addLink(treePath(row,0), treePath(col,0));
            pattern_set[row].insert(col);
          }

        std::chrono::duration<double, std::micro> duration = std::chrono::system_clock::now() - start_t;
        // std::cout << size << " \t" << max_entries << " \t" << entries_per_row << " \t" << duration.count() << std::endl;

        pattern.sort();

        auto row_nnzs = pattern.rowNonZeros();
        for (std::size_t row = 0; row != row_nnzs.size(); ++row) {
          EXPECT_EQ(row_nnzs[row], pattern_set[row].size());
          std::size_t cols = 0;
          for(auto it = pattern.columnBegin(row); it != pattern.columnEnd(row); ++it) {
            EXPECT_TRUE(pattern_set[row].contains(*it)) << "Expected value " << *it;
            pattern_set[row].erase(*it);
            ++cols;
          }
          EXPECT_TRUE(pattern_set[row].empty());
          EXPECT_EQ(cols, row_nnzs[row]);
        }

        EXPECT_FALSE(pattern < pattern);
      }
    }
  }
}

TEST(TestPatterns, TestRandomBlockedBlockedPattern) {
  // in this case we check the outer and inner block
  // that is, multi indices for the inner block, are also random


// test different fill and buffer rations to go to all code paths of the pattern
#ifndef NDEBUG
std::size_t max_size_shift = 4;
#else
std::size_t max_size_shift = 6;
#endif

  for (std::size_t size = 2<<1; size != 2u<<max_size_shift; size <<= 1) {
    for (auto fill_ratio : {.25, .5, .75, 1.}) {
      for (auto buffer_ratio : {0.5, 1.0, 1.5}) {
        std::size_t max_entries = static_cast<std::size_t>(fill_ratio*std::min<std::size_t>(size,10));
        std::size_t entries_per_row = static_cast<std::size_t>(buffer_ratio*max_entries);

        BlockedSizeProvider sp{size};
        using LeafPattern = Dune::PDELab::LeafSparsePattern<BlockedSizeProvider,BlockedSizeProvider>;

        Dune::PDELab::BlockedSparsePattern<LeafPattern> pattern{sp, {}, sp, {},[=](auto,auto){ return entries_per_row; }};


        std::uniform_int_distribution<std::size_t> unif(0, sp.size({})-1);
        std::default_random_engine re;
        re.seed(0);
        std::set<std::array<std::size_t,4>> pattern_set;

        auto start_t = std::chrono::system_clock::now();

        using Dune::TypeTree::treePath;
        for (std::size_t row_out = 0; row_out != sp.size({}); ++row_out) {
          for (std::size_t i_out = 0; i_out != max_entries; ++i_out) {
            std::size_t col_out = unif(re);
            for (std::size_t row_in = 0; row_in != sp.size({}); ++row_in) {
              for (std::size_t i_in = 0; i_in != max_entries; ++i_in) {
                std::size_t col_in = unif(re);
                auto row = treePath(row_out,row_in);
                auto col = treePath(col_out,col_in);
                pattern.addLink(row, col);
                pattern_set.insert({row_out, row_in, col_out, col_in});
              }
            }
          }
        }

        std::chrono::duration<double, std::micro> duration = std::chrono::system_clock::now() - start_t;
        // std::cout << size << " \t" << max_entries << " \t" << entries_per_row << " \t" << duration.count() << std::endl;

        pattern.sort();

        auto row_out_nnzs = pattern.rowNonZeros();
        for (std::size_t row_out = 0; row_out != row_out_nnzs.size(); ++row_out) {
          EXPECT_LE(row_out_nnzs[row_out], max_entries);
          std::size_t cols_out = 0;
          for(auto it_out = pattern.columnBegin(row_out); it_out != pattern.columnEnd(row_out); ++it_out) {
            auto col_out = *it_out;
            auto& pattern_in = it_out.pattern();

            auto row_in_nnzs = pattern_in.rowNonZeros();
            for (std::size_t row_in = 0; row_in != row_in_nnzs.size(); ++row_in) {
              // EXPECT_LE(row_in_nnzs[row_in], max_entries); // since col_out may be randomly repeated, we cannot guarantee this
              std::size_t cols_in = 0;
              for(auto it_in = pattern_in.columnBegin(row_in); it_in != pattern_in.columnEnd(row_in); ++it_in) {
                auto col_in = *it_in;
                auto entry = std::array<std::size_t,4>{row_out, row_in, col_out, col_in};
                EXPECT_TRUE(pattern_set.contains(entry)) << "Wrong value in pattern " << row_out << " " << row_in << " " << col_out << " " << col_in << std::endl;
                pattern_set.erase(entry);
                ++cols_in;
              }
              EXPECT_EQ(cols_in, row_in_nnzs[row_in]);
            }
            ++cols_out;
          }
          EXPECT_EQ(cols_out, row_out_nnzs[row_out]);
        }

        EXPECT_TRUE(pattern_set.empty());
        for (auto [row_out, row_in, col_out ,col_in] : pattern_set)
          ADD_FAILURE() << "Missing pattern set " << row_out << " " << row_in << " " << col_out << " " << col_in << std::endl;
        EXPECT_FALSE(pattern < pattern);
      }
    }
  }
}
