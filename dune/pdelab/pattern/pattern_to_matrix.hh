#ifndef DUNE_PDELAB_PATTERN_PATTERN_TO_MATRIX_HH
#define DUNE_PDELAB_PATTERN_PATTERN_TO_MATRIX_HH

#include <numeric>
#include <map>

namespace Dune::PDELab::inline Experimental {

template<class Pattern, class Mat>
void patternToMatrix(
  const Pattern& pattern,
  Mat& matrix,
  typename Pattern::RowSizeProvider::MultiIndex row_prefix = {},
  typename Pattern::ColSizeProvider::MultiIndex col_prefix = {})
{
  using size_type = typename Pattern::size_type;
  auto row_size = pattern.rowSizeProvider().size(row_prefix);
  auto col_size = pattern.colSizeProvider().size(col_prefix);
  if constexpr (requires {matrix.endindices();}) { // sparse matrices
    auto row_nnz = pattern.rowNonZeros();
    size_type nnz = std::accumulate(begin(row_nnz), end(row_nnz), size_type{0});
    size_type longest_row = 0;

    matrix.setSize(row_size, col_size, 0 /*nnz*/);
    matrix.setBuildMode(Mat::random);

    for (size_type row = 0; row != pattern.rows(); ++row) {
      nnz -= row_nnz[row];
      longest_row = std::max(longest_row, row_nnz[row]);
      matrix.setrowsize(row, row_nnz[row]);
    }
    assert(nnz == 0);
    matrix.endrowsizes();

    for (size_type row = 0; row != pattern.rows(); ++row)
      matrix.setIndices(row,pattern.columnBegin(row),pattern.columnEnd(row));
  } else {
    if constexpr (requires {matrix.setSize(row_size, col_size);})
      matrix.setSize(row_size, col_size);
  }

  if constexpr (not Pattern::isLeaf) {
    auto comp_val = [](auto a, auto b){return *a < *b;};
    using SubPattern = typename Pattern::SubPattern;
    using Block = typename Mat::block_type;
    std::map<SubPattern const*,Block const *, decltype(comp_val)> reuse_map(comp_val);
    // allocate data array
    if constexpr (requires {matrix.endindices();})  // sparse matrices
      matrix.endindices();
    // add sub patterns on sub matrices
    typename Pattern::RowSizePrefix row_sub_prefix = row_prefix;
    row_sub_prefix.push_back(0);
    typename Pattern::RowSizePrefix col_sub_prefix = col_prefix;
    col_sub_prefix.push_back(0);

    for (auto row_it = matrix.begin(); row_it != matrix.end(); ++row_it) {
      auto row = row_it.index();
      row_sub_prefix.back() = row;
      auto pattern_it = pattern.columnBegin(row);
      for(auto col_it = row_it->begin(); col_it != row_it->end(); ++col_it, ++pattern_it) {
        assert(col_it.index() == *pattern_it);
        col_sub_prefix.back() = col_it.index();
        auto& sub_container = *col_it;

        // we gain ownership of pattern storage
        const auto& sub_pattern = pattern_it.pattern();
        // if pattern storage is *not equivalent* to another pattern
        // we insert it in the map, otherwise it is discarded at the end of this scope
        auto [reuse_it, block_is_new] = reuse_map.insert({&sub_pattern,nullptr});
        if (block_is_new) {
          // not equivalent pattern in map, we need to allocate matrix from stracth
          patternToMatrix(sub_pattern, sub_container, row_sub_prefix, col_sub_prefix);
          reuse_it->second = &sub_container;
        } else {
          // copy construct bcrs matrix and thus share underlying pattern!
          assert(reuse_it->second);
          sub_container = *reuse_it->second;
        }
      }
      assert(pattern_it == pattern.columnEnd(row));
    }
    // pattern.clear();
  } else {
    // free temporary index storage in pattern before allocating data array in matrix
    // if (row_prefix.size() == 0 and col_prefix.size() == 0)
    //   pattern.clear();
    // allocate data array
    if constexpr (requires {matrix.endindices();})  // sparse matrices
      matrix.endindices();
  }
}


} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_PATTERN_PATTERN_TO_MATRIX_HH
