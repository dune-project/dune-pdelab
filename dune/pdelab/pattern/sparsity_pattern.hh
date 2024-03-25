#ifndef DUNE_PDELAB_PATTERN_SPARSITY_PATTERN_HH
#define DUNE_PDELAB_PATTERN_SPARSITY_PATTERN_HH

#include <dune/pdelab/concepts/multiindex.hh>
#include <dune/pdelab/concepts/size_provider.hh>

#include <dune/common/iteratorfacades.hh>
#include <dune/common/exceptions.hh>

#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <vector>
#include <memory>
#include <limits>
#include <cassert>
#include <mutex>
#include <functional>

namespace Dune::PDELab::inline Experimental {


template<Concept::SizeProvider RowSizeProvider_, Concept::SizeProvider ColSizeProvider_, class SubPattern_ = void*>
class SparsePatternBase
{

public:
  //! SparsePattern cannot contain nested subpatterns. This entry is only
  //! required for TMP purposes.
  using SubPattern = SubPattern_;
  static const bool isLeaf = std::is_same<SubPattern, void*>{};

  using RowSizeProvider = RowSizeProvider_;
  using ColSizeProvider = ColSizeProvider_;

  using RowSizePrefix = typename RowSizeProvider::SizePrefix;
  using ColSizePrefix = typename ColSizeProvider::SizePrefix;

  //! size type used by SparsePattern.
  using size_type = std::size_t;

  using Indices = std::conditional_t<
    isLeaf,
    std::vector<size_type>,
    std::vector<std::pair<size_type, SubPattern>>>;

  using Overflow = std::conditional_t<
    isLeaf,
    std::set<std::pair<size_type, size_type>>,
    std::map<std::pair<size_type, size_type>, std::shared_ptr<SubPattern>>>;

public:

  SparsePatternBase() = default;
  SparsePatternBase(SparsePatternBase&&) = default;
  SparsePatternBase(const SparsePatternBase&) = default;

  SparsePatternBase& operator=(SparsePatternBase&&) = default;
  SparsePatternBase& operator=(const SparsePatternBase&) = default;

  //! Constructs a SparsePattern for the given pair of orderings and reserves
  //! space for the provided average number of entries per row.
  /**
   * \param row_size_provider    Ordering describing the row structure
   * \param col_size_provider    Ordering describing the column structure
   * \param entries_per_row An estimate of the average number of entries per
   * row.
   */
  SparsePatternBase(const RowSizeProvider& row_size_provider,
                      typename RowSizeProvider::SizePrefix row_prefix,
                      const ColSizeProvider& col_size_provider,
                      typename ColSizeProvider::SizePrefix col_prefix,
                      std::function<size_type(RowSizePrefix, ColSizePrefix)> fentries_per_row)
    : _row_size_provider(&row_size_provider)
    , _col_size_provider(&col_size_provider)
    , _fentries_per_row{fentries_per_row}
    , _row_prefix(row_prefix)
    , _col_prefix(col_prefix)
    , _rows(_row_size_provider->size(_row_prefix))
    , _cols(_col_size_provider->size(_col_prefix))
    // The +1 is to ensure that requested entries per row have O(er) complexity
    , _entries_per_row(std::min(_cols, _fentries_per_row(_row_prefix, _col_prefix) + 1))
    , _sorted(false)
  {
    if constexpr (isLeaf)
      _indices.resize(_rows * _entries_per_row, empty_col);
    else
      _indices.resize(_rows * _entries_per_row, { empty_col, SubPattern{} });
  }

protected:
  static const size_type& getIndexCol(const typename Indices::value_type& t)
  {
    if constexpr (isLeaf)
      return t;
    else
      return t.first;
  }

  static size_type& getIndexCol(typename Indices::value_type& t)
  {
    if constexpr (isLeaf)
      return t;
    else
      return t.first;
  }

  static const size_type& getOverflowCol(
    const typename Overflow::const_iterator& t)
  {
    if constexpr (isLeaf)
      return t->second;
    else
      return t->first.second;
  }

  static const size_type& getOverflowRow(
    const typename Overflow::const_iterator& t)
  {
    if constexpr (isLeaf)
      return t->first;
    else
      return t->first.first;
  }

public:
  //! Returns a vector with the size of all rows in the pattern.
  std::vector<size_type> rowNonZeros() const
  {
    std::vector<size_type> row_sizes(_rows, 0);
    auto idx_it = std::begin(_indices);
    auto idx_end = std::begin(_indices) + _entries_per_row;
    auto ovf_it = std::begin(_overflow);
    auto ovf_end = std::end(_overflow);
    for (size_type row = 0; row < _rows; ++row) {
      size_type row_size = 0;
      // count non-empty column entries
      while (idx_it != idx_end and getIndexCol(*idx_it) != empty_col)
        ++idx_it, ++row_size;
      // add overflow entries
      while (ovf_it != ovf_end and getOverflowRow(ovf_it) == row)
        ++ovf_it, ++row_size;
      // write row size int vector and advance to next row
      row_sizes[row] = row_size;
      // advance index range to next row
      idx_it = idx_end;
      std::advance(idx_end, _entries_per_row);
    }

    return row_sizes;
  }

  std::size_t rows() const {
    return _rows;
  }

  std::size_t cols() const {
    return _cols;
  }

  //! Iterator over all column indices for a given row
  struct const_column_iterator
    : public ForwardIteratorFacade<const_column_iterator, const size_type>
  {

    const size_type& dereference() const
    {
      if (_in_overflow)
        return getOverflowCol(_ovf_it);
      else
        return getIndexCol(*_idx_it);
    }

    const SubPattern& pattern() const
    {
      static_assert(not isLeaf);
      if (_in_overflow)
        return *_ovf_it->second;
      else
        return _idx_it->second;
    }

    void increment()
    {
      if (_in_overflow) {
        if (++_ovf_it == _ovf_end or getOverflowRow(_ovf_it) != _row) {
          // we have exhausted the row, invalidate iterator
          _at_end = true;
        }
      } else {
        if (_idx_it != _idx_end)
          ++_idx_it;
        if (_idx_it == _idx_end || getIndexCol(*_idx_it) == empty_col) {
          _in_overflow = true;
          // we have exhausted the row, invalidate iterator
          if (_ovf_it == _ovf_end || getOverflowRow(_ovf_it) > _row)
            _at_end = true;
        }
      }
    }

    bool equals(const const_column_iterator& other) const
    {
      if (_row != other._row)
        return false;
      if (_at_end || other._at_end)
        return _at_end && other._at_end;
      if (_in_overflow)
        return _ovf_it == other._ovf_it;
      else
        return _idx_it == other._idx_it;
    }

    const_column_iterator(const SparsePatternBase& pattern, size_type row, bool at_end)
      : _row(row)
      , _in_overflow(false)
      , _at_end(at_end)
      , _idx_it(pattern._indices.begin() + row * pattern._entries_per_row)
      , _idx_end(pattern._indices.begin() + (row + 1) * pattern._entries_per_row)
      , _ovf_it(pattern._overflow.lower_bound(std::make_pair(row, 0)))
      , _ovf_end(pattern._overflow.end())
    {
      assert(row < pattern._rows);
      // catch corner case with completely empty row
      if ((!_at_end) &&
          (_idx_it == _idx_end || getIndexCol(*_idx_it) == empty_col)) {
        _in_overflow = true;
        _at_end = _ovf_it == _ovf_end || getOverflowRow(_ovf_it) != _row;
      }
    }

    size_type _row;
    bool _in_overflow;
    bool _at_end;
    typename Indices::const_iterator _idx_it;
    typename Indices::const_iterator _idx_end;
    typename Overflow::const_iterator _ovf_it;
    const typename Overflow::const_iterator _ovf_end;
  };

  //! Returns an iterator to the first column index of row i.
  const_column_iterator columnBegin(size_type row) const
  {
    return const_column_iterator(*this, row, false);
  }

  //! Returns an iterator past the last column index of row i.
  const_column_iterator columnEnd(size_type row) const
  {
    return const_column_iterator(*this, row, true);
  }


  //! Iterator over all column indices for a given row, unique but in arbitrary
  //! order.
  struct column_iterator
    : public ForwardIteratorFacade<column_iterator, const size_type>
  {

    const size_type& dereference() const
    {
      if (_in_overflow)
        return getOverflowCol(_ovf_it);
      else
        return getIndexCol(*_idx_it);
    }

    SubPattern& pattern() const
    {
      static_assert(not isLeaf);
      if (_in_overflow)
        return *_ovf_it->second;
      else
        return _idx_it->second;
    }

    void increment()
    {
      if (_in_overflow) {
        if (++_ovf_it == _ovf_end or getOverflowRow(_ovf_it) != _row) {
          // we have exhausted the row, invalidate iterator
          _at_end = true;
        }
      } else {
        if (_idx_it != _idx_end)
          ++_idx_it;
        if (_idx_it == _idx_end || getIndexCol(*_idx_it) == empty_col) {
          _in_overflow = true;
          // we have exhausted the row, invalidate iterator
          if (_ovf_it == _ovf_end || getOverflowRow(_ovf_it) > _row)
            _at_end = true;
        }
      }
    }

    bool equals(const column_iterator& other) const
    {
      if (_row != other._row)
        return false;
      if (_at_end || other._at_end)
        return _at_end && other._at_end;
      if (_in_overflow)
        return _ovf_it == other._ovf_it;
      else
        return _idx_it == other._idx_it;
    }

    column_iterator(SparsePatternBase& pattern, size_type row, bool at_end)
      : _row(row)
      , _in_overflow(false)
      , _at_end(at_end)
      , _idx_it(pattern._indices.begin() + row * pattern._entries_per_row)
      , _idx_end(pattern._indices.begin() + (row + 1) * pattern._entries_per_row)
      , _ovf_it(pattern._overflow.lower_bound(std::make_pair(row, 0)))
      , _ovf_end(pattern._overflow.end())
    {
      assert(row < pattern._rows);
      // catch corner case with completely empty row
      if ((!_at_end) &&
          (_idx_it == _idx_end || getIndexCol(*_idx_it) == empty_col)) {
        _in_overflow = true;
        _at_end = _ovf_it == _ovf_end || getOverflowRow(_ovf_it) != _row;
      }
    }

    size_type _row;
    bool _in_overflow;
    bool _at_end;
    typename Indices::iterator _idx_it;
    typename Indices::iterator _idx_end;
    typename Overflow::iterator _ovf_it;
    const typename Overflow::iterator _ovf_end;
  };

  //! Returns an iterator to the first column index of row i.
  column_iterator columnBegin(size_type row)
  {
    return column_iterator(*this, row, false);
  }

  //! Returns an iterator past the last column index of row i.
  column_iterator columnEnd(size_type row)
  {
    return column_iterator(*this, row, true);
  }



protected:
  bool is_sorted(size_type row) const {
    // get iterators for current row
    auto begin = std::begin(_indices) + _entries_per_row * row;
    auto end = begin + _entries_per_row;
#ifndef NDEBUG
    // we always keep range sorted when the it's filled
    auto last = std::prev(end);
    return (getIndexCol(*last) != empty_col);
#else
    return std::is_sorted(begin, end, [](auto& a, auto& b) {
      return getIndexCol(a) < getIndexCol(b);
    });
#endif
  }

  void sort(size_type row)
  {
    // get iterators for current row
    auto begin = std::begin(_indices) + _entries_per_row * row;
    auto end = begin + _entries_per_row;
    std::sort(begin, end, [](auto& a, auto& b) {
      return getIndexCol(a) < getIndexCol(b);
    });
  }

public:

  void sort()
  {
    if (not _sorted) {
      for (size_type row = 0; row < _rows; ++row) {
        if (not is_sorted(row))
          sort(row);
        if constexpr (not isLeaf)
          for (auto it = columnBegin(row); it != columnEnd(row); ++it)
            it.pattern().sort();
      }
    }
    _sorted = true;
  }


  //! Remove all links stored in the pattern
  void clear()
  {
    if constexpr (isLeaf)
      _indices.assign(_indices.size(), empty_col);
    else
      _indices.assign(_indices.size(), { empty_col, SubPattern{} });
    _overflow.clear();
  }

  size_type entriesPerRow() const
  {
    return _entries_per_row;
  }

  size_type overflowCount() const
  {
    return _overflow.size();
  }

  const RowSizeProvider& rowSizeProvider() const
  {
    assert(_row_size_provider);
    return *_row_size_provider;
  }

  const ColSizeProvider& colSizeProvider() const
  {
    assert(_col_size_provider);
    return *_col_size_provider;
  }

  bool operator<(const SparsePatternBase& other) const
  {
    if (not (_sorted and other._sorted))
      DUNE_THROW(InvalidStateException, "Patterns must be sorted before before comparing them");
    assert(_sorted and other._sorted);

    if (this == &other) // irreflexivity
      return false;

    bool r = _rows < other._rows or _cols < other._cols or _entries_per_row < other._entries_per_row or _overflow < other._overflow;
    /* raw lexicopgraphic compare on indices does not work because entries
    * filled with empty_col are also counted as indices (not as a size mismacth).
    * Now, since empty_col is bigger than any other number, it counts as
    * "major than", countrary to the behavior of vectors with different lenght
    * which is what we want.
    */
    r |= std::lexicographical_compare(
              _indices.begin(), _indices.end(),
              other._indices.begin(), other._indices.end(),
              [&](auto& _a, auto& _b){
                size_type a = getIndexCol(_a), b = getIndexCol(_b);
                if constexpr (std::is_unsigned_v<size_type> and std::is_integral_v<size_type>) {
                  // In this case we leverage on the unsigned integral types
                  // where the max value is folded to zero by simply adding one.
                  return a+1 < b+1;
                } else {
                  // In this case we explicitely check for (a,b) to be empty_col
                  // notice that this lambda shall give weak ordered guarantees
                  if (a == empty_col) // special case
                    return (b != a); // notice irreflexivity: (b == a) => !(a < b)
                  else
                    return a < b;
                }
              }
            );
    if (r)
      return true;
    if constexpr (not isLeaf) {
      for (std::size_t row = 0; row < _rows; ++row) {
        auto it = this->columnBegin(row);
        auto other_it = other.columnBegin(row);
        while (it != this->columnEnd(row) and other_it != other.columnEnd(row)) {
          if (it.pattern() < other_it.pattern())
            return true;
          ++it, ++other_it;
        }
        if (other_it != other.columnEnd(row))
          return true;
      }
    }
    return false;
  }

  void addLink(Concept::MultiIndex auto row_suffix, Concept::MultiIndex auto col_suffix)
  {
    // in case of empty indices, no link is needed  FIXME???
    // if (row_suffix.size() == 0 or col_suffix.size() == 0)
    //   return;

    // extract block indices for current level
    size_type row = front(row_suffix);
    assert(row < this->rowSizeProvider().size(_row_prefix));
    size_type col = front(col_suffix);
    assert(col < this->colSizeProvider().size(_col_prefix));

    // get iterators for current row
    auto start = _indices.begin();
    auto begin = start + _entries_per_row * row;
    auto end = begin + _entries_per_row;
    auto last = std::prev(end);

    // is_sorted <=> idx range is filled left to right, so we know when is full
    bool idx_full = (getIndexCol(*last) != empty_col);

    SubPattern* sub_pattern_ptr = nullptr;
    auto make_sub_pattern = [&]() {
      if constexpr (isLeaf) {
        return nullptr;
      } else {
        auto sub_row_prefix = push_back(this->_row_prefix, row);
        auto sub_col_prefix = push_back(this->_col_prefix, col);
        return SubPattern{this->rowSizeProvider(), sub_row_prefix,
                          this->colSizeProvider(), sub_col_prefix,
                          _fentries_per_row};
      }
    };

    if (not idx_full) {
      // if idx range is not full, cols are not sorted -> linear search
      // looking up a position to set the column index. stop if found.
      auto it = find_if(begin, end, [col](auto& k) {
        return getIndexCol(k) == col or getIndexCol(k) == empty_col;
      });
      assert(it != end);

      if constexpr (not isLeaf) {
        if (getIndexCol(*it) == empty_col)
          it->second = make_sub_pattern();
        sub_pattern_ptr = &(it->second);
      }
      getIndexCol(*it) = col;

      if (it == last) {
        this->sort(row); // we just filled idx range -> sort it
        if constexpr (not isLeaf) {
          // because positions in _indices changes, we need to find sub_pattern_ptr again
          auto it = std::lower_bound(begin, end, col, [](auto& a, auto& b) { return getIndexCol(a) < b; });
          assert(getIndexCol(*it) == col);
          sub_pattern_ptr = &(it->second);
        }
      }
      else
        _sorted = false; // we put value in first possible spot, cannot guarantee ordering
    } else {
      //! @invariant: [begin,end) is sorted and values are partially ordered wrt
      // overflow. Thus, if idx range is full, cols are sorted -> binary search
      assert(is_sorted(row));
      auto it = std::lower_bound(
        begin, end, col, [](auto& a, auto& b) { return getIndexCol(a) < b; });

      if (it != end) {
        if (getIndexCol(*it) == col) {
          // column is in idx range, get sub pattern
          if constexpr (not isLeaf)
            sub_pattern_ptr = &(it->second);
        } else {
          // not found but col order is in idx range
          // store last value in range to store it later in the overflow
          auto overflow_content = std::move(*prev(end));

          // shift values in idx range to the right
          std::move_backward(it, prev(end), end);

          // send last value of the shift to overflow
          // and set row value in the newly created space
          if constexpr (isLeaf) {
            _overflow.insert({ row, overflow_content });
            *it = col;
          } else {
            _overflow.insert({ { row, overflow_content.first }, std::make_shared<SubPattern>(std::move(overflow_content.second)) });
            *it = { col, make_sub_pattern() };
            sub_pattern_ptr = &it->second;
          }
        }
      } else {
        // row belongs to overflow partition
        if constexpr (isLeaf) {
          _overflow.insert({ row, col });
        } else {
          auto [ovf_it, _] = _overflow.insert({ { row, col }, std::make_shared<SubPattern>(make_sub_pattern()) });
          sub_pattern_ptr = ovf_it->second.get();
        }
      }
    }
    if constexpr (not isLeaf) {
      assert(sub_pattern_ptr);
      sub_pattern_ptr->addLink(pop_front(row_suffix), pop_front(col_suffix));
    }
  }

protected:
  //! Marker value indicating an empty array entry.
  inline static constexpr size_type empty_col = std::numeric_limits<size_type>::max();

  RowSizeProvider const * _row_size_provider = nullptr;
  ColSizeProvider const * _col_size_provider = nullptr;

  Indices _indices;
  Overflow _overflow;

  std::function<size_type(RowSizePrefix, ColSizePrefix)> _fentries_per_row;

  RowSizePrefix _row_prefix;
  ColSizePrefix _col_prefix;

  size_type _rows, _cols, _entries_per_row;

  bool _sorted;
};

template<Concept::SizeProvider RowSizeProvider, Concept::SizeProvider ColSizeProvider>
class LeafSparsePattern
  : public SparsePatternBase<RowSizeProvider, ColSizeProvider>
{
  using Base = SparsePatternBase<RowSizeProvider, ColSizeProvider>;

public:

  using size_type = typename Base::size_type;

  using Base::Base;

  LeafSparsePattern(const RowSizeProvider& row_size_provider,
                      typename RowSizeProvider::SizePrefix row_prefix,
                      const ColSizeProvider& col_size_provider,
                      typename ColSizeProvider::SizePrefix col_prefix,
                      size_type entries_per_row)
  : Base{row_size_provider, row_prefix, col_size_provider, col_prefix, [=](auto,auto){return entries_per_row;}}
  {}
};

template<class SubPattern>
class BlockedSparsePattern
  : public SparsePatternBase<typename SubPattern::RowSizeProvider,
                             typename SubPattern::ColSizeProvider,
                             SubPattern>
{
  using Base = SparsePatternBase<typename SubPattern::RowSizeProvider,
                                 typename SubPattern::ColSizeProvider,
                                SubPattern>;
public:
  using Base::Base;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_PATTERN_SPARSITY_PATTERN_HH
