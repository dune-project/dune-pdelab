// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_BCRSPATTERN_HH
#define DUNE_PDELAB_BACKEND_ISTL_BCRSPATTERN_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <utility>
#include <vector>
#include <algorithm>
#include <set>

#include <dune/common/iteratorfacades.hh>

#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>

namespace Dune {
  namespace PDELab {
    namespace ISTL {

      //! Pattern builder for generic BCRS-like sparse matrices.
      /**
       * BCRSPattern is a pattern builder for unstructured sparse matrices
       * for operators mapping from a vector that conforms to RowOrdering to a vector
       * that conforms to ColOrdering.
       *
       * BCRSPattern has much better runtime performance and requires far less memory
       * than the older pattern constructon method in PDELab. By letting the user specify
       * the average number of nonzeroes per row, it is possible to use a more efficient
       * array-based storage scheme for the majority of the pattern entries, only using
       * expensive map-like lookups for those entries that exceed that average.
       *
       * BCRSPattern requires a recent version of the BCRSMatrix with support for row-wise
       * setting of column indices and split allocation of column index and data arrays.
       *
       * Note that unlike the implicit construction mode of the BCRSMatrix itself, this
       * pattern builder will neither throw an exception if the number of nonzeroes was set
       * too low nor retain excess memory if it was set too high after the pattern construction
       * is complete. Performance will degrade if the user-provided estimate is too far away
       * from the real value.
       */
      template<typename RowOrdering, typename ColOrdering>
      class BCRSPattern
      {

      public:

        //! size type used by BCRSPattern.
        typedef typename RowOrdering::Traits::size_type size_type;

        //! BCRSPattern cannot contain nested subpatterns. This entry is only required for TMP purposes.
        typedef void SubPattern;

      private:

        //! Marker value indicating an empty array entry.
        static const size_type empty = ~static_cast<size_type>(0);

        //! Functor for looking up a column index within a row.
        /**ss
         * Looking up column indices requires a special comparison iterator,
         * as we want to either return the position of the actual index if it
         * has already been inserted or of the first empty matrix entry that we
         * can use to store the index.
         *
         * \warning Only use for sequential search algorithms, binary searches
         *          will not work correctly!
         */
        struct PaddedColumnCriterion
        {

          bool operator()(size_type k) const
          {
            return k == _j || k == empty;
          }

          PaddedColumnCriterion(size_type j)
            : _j(j)
          {}

          const size_type _j;

        };


        typedef typename std::vector<size_type>::iterator IndicesIterator;
        typedef typename std::set<std::pair<size_type,size_type> >::iterator OverflowIterator;

        typedef typename std::vector<size_type>::const_iterator ConstIndicesIterator;
        typedef typename std::set<std::pair<size_type,size_type> >::const_iterator ConstOverflowIterator;

      public:

        //! Add a link between the row indicated by ri and the column indicated by ci.
        template<typename RI, typename CI>
        void add_link(const RI& ri, const CI& ci)
        {
          // extract block indices for current level
          size_type i = ri.back();
          size_type j = ci.back();

          IndicesIterator start = _indices.begin();
          IndicesIterator begin = start + _entries_per_row*i;
          IndicesIterator end = start + _entries_per_row*(i+1);

          // Does the entry (i,j) already exist?
          IndicesIterator it = std::find_if(begin,end,PaddedColumnCriterion(j));
          if (it != end)
            {
              // Yes, just write out j. This does the right thing regardless of whether
              // it points at the correct column value or at an empty entry.
              *it = j;
            }
          else
            {
              // The row is already full -> spill into map
              _overflow.insert(std::make_pair(i,j));
            }
        }

#ifndef DOXYGEN

        // implementation detail for nested patterns
        template<typename RI, typename CI>
        void recursive_add_entry(const RI& ri, const CI& ci)
        {
          add_link(ri,ci);
        }

#endif //DOXYGEN

        //! Stream the sizes of all rows into the output iterator rit.
        template<typename I>
        void sizes(I rit) const
        {
          ConstIndicesIterator it = _indices.begin();
          ConstIndicesIterator end = _indices.begin() + _entries_per_row;
          ConstOverflowIterator oit = _overflow.begin();
          ConstOverflowIterator oend = _overflow.end();
          for (size_type i = 0; i < _row_ordering.blockCount(); ++i, ++rit, end+=_entries_per_row)
            {
              size_type s = 0;
              // count non-empty column entries, break when first empty one is found.
              for (; it != end; ++it)
                {
                  if (*it == empty)
                    break;
                  ++s;
                }
              it = end;
              // add overflow entries
              for (; oit != oend && oit->first == i; ++oit)
                ++s;
              *rit = s;
            }
        }

        //! Returns a vector with the size of all rows in the pattern.
        std::vector<size_type> sizes() const
        {
          std::vector<size_type> r(_row_ordering.blockCount());
          sizes(r.begin());
          return std::move(r);
        }

        //! Iterator over all column indices for a given row, unique but in arbitrary order.
        struct iterator
          : public ForwardIteratorFacade<iterator, const size_type>
        {

#ifndef DOXYGEN

          const size_type& dereference() const
          {
            if (_in_overflow)
              return _oit->second;
            else
              return *_it;
          }

          void increment()
          {
            if (_in_overflow)
              {
                if (_oit != _oend)
                  ++_oit;
                if (_oit->first == _row)
                  return;
                // we have exhausted the row, invalidate iterator
                _at_end = true;
              }
            else
              {
                if (_it != _end)
                  ++_it;
                if (_it == _end || *_it == empty)
                  {
                    _in_overflow = true;
                    // we have exhausted the row, invalidate iterator
                    if (_oit == _oend || _oit->first > _row)
                      _at_end = true;
                  }
              }
          }

          bool equals(const iterator& other) const
          {
            if (_row != other._row)
              return false;
            if (_at_end || other._at_end)
              return _at_end && other._at_end;
            if (_in_overflow)
              return _oit == other._oit;
            else
              return _it == other._it;
          }

          iterator(const BCRSPattern& p, size_type row, bool at_end)
            : _row(row)
            , _in_overflow(false)
            , _at_end(at_end)
            , _it(p._indices.begin() + row * p._entries_per_row)
            , _end(p._indices.begin() + (row+1) * p._entries_per_row)
            , _oit(p._overflow.lower_bound(std::make_pair(row,0)))
            , _oend(p._overflow.end())
          {
            // catch corner case with completely empty row
            if ((!_at_end) && (_it == _end || *_it == empty))
              {
                _in_overflow = true;
                _at_end = _oit == _oend || _oit->first != _row;
              }
          }

          size_type _row;
          bool _in_overflow;
          bool _at_end;
          typename std::vector<size_type>::const_iterator _it;
          typename std::vector<size_type>::const_iterator _end;
          typename std::set<std::pair<size_type,size_type> >::const_iterator _oit;
          const typename std::set<std::pair<size_type,size_type> >::const_iterator _oend;

#endif // DOXYGEN

        };

        //! Returns an iterator to the first column index of row i.
        iterator begin(size_type i) const
        {
          return iterator(*this,i,false);
        }

        //! Returns an iterator past the last column index of row i.
        iterator end(size_type i) const
        {
          return iterator(*this,i,true);
        }

        //! Constructs a BCRSPattern for the given pair of orderings and reserves space for the provided average number of entries per row.
        /**
         * \param row_ordering    Ordering describing the row structure
         * \param col_ordering    Ordering describing the column structure
         * \param entries_per_row An estimate of the average number of entries per row.
         */
        BCRSPattern(const RowOrdering& row_ordering, const ColOrdering& col_ordering, size_type entries_per_row)
          : _row_ordering(row_ordering)
          , _col_ordering(col_ordering)
          , _entries_per_row(entries_per_row)
          , _indices(row_ordering.blockCount()*entries_per_row,size_type(empty))
        {}

        const RowOrdering& rowOrdering() const
        {
          return _row_ordering;
        }

        const ColOrdering& colOrdering() const
        {
          return _row_ordering;
        }

        //! Discard all internal data.
        /**
         * The purpose of this method is to release all internal memory before calling
         * BCRSMatrix::endindices(). That way, the matrix creation process never consumes
         * substantially more memory as required by the matrix after construction, as the
         * second copy of the column indices is about as large as the data array.
         */
        void clear()
        {
          _indices = std::vector<size_type>();
          _overflow = std::set<std::pair<size_type,size_type> >();
        }

        size_type entriesPerRow() const
        {
          return _entries_per_row;
        }

        size_type overflowCount() const
        {
          return _overflow.size();
        }

      private:

        const RowOrdering& _row_ordering;
        const ColOrdering& _col_ordering;
        const size_type _entries_per_row;

        std::vector<size_type> _indices;
        std::set<std::pair<size_type,size_type> > _overflow;

      };


      //! Pattern builder for nested hierarchies of generic BCRS-like sparse matrices.
      /**
       * NestedPattern contains a dense set of subpatterns for each matrix block. Those
       * blocks can be nested (i.e. be NestedPatterns again) or BCRSPattern instances.
       */

      template<typename RowOrdering, typename ColOrdering, typename SubPattern_>
      class NestedPattern
      {

      public:

        //! The pattern type used for each block.
        typedef SubPattern_ SubPattern;

        //! size type used by NestedPattern.
        typedef typename SubPattern::size_type size_type;

        //! Add a link between the row indicated by ri and the column indicated by ci.
        /**
         * This method just forwards the call to the relevant block as indicated by the
         * tail members of ri and ci.
         */
        template<typename RI, typename CI>
        void add_link(const RI& ri, const CI& ci)
        {
          recursive_add_entry(ri.view(),ci.view());
        }

#ifndef DOXYGEN

        template<typename RI, typename CI>
        void recursive_add_entry(const RI& ri, const CI& ci)
        {
          _sub_patterns[ri.back() * _col_ordering.blockCount() + ci.back()].recursive_add_entry(ri.back_popped(),ci.back_popped());
        }

#endif // DOXYGEN

        template<typename EntriesPerRow>
        NestedPattern(const RowOrdering& row_ordering, const ColOrdering& col_ordering, const EntriesPerRow& entries_per_row)
          : _row_ordering(row_ordering)
          , _col_ordering(col_ordering)
        {
          size_type rows = row_ordering.blockCount();
          size_type cols = col_ordering.blockCount();
          for (size_type i = 0; i < rows; ++i)
            for (size_type j = 0; j < cols; ++j)
              _sub_patterns.push_back(
                SubPattern(
                  _row_ordering.childOrdering(i),
                  _col_ordering.childOrdering(j),
                  entries_per_row[i][j]
                  )
                );
        }

        NestedPattern(const RowOrdering& row_ordering, const ColOrdering& col_ordering, size_type entries_per_row)
          : _row_ordering(row_ordering)
          , _col_ordering(col_ordering)
        {
          size_type rows = row_ordering.blockCount();
          size_type cols = col_ordering.blockCount();
          for (size_type i = 0; i < rows; ++i)
            for (size_type j = 0; j < cols; ++j)
              _sub_patterns.push_back(
                SubPattern(
                  _row_ordering.childOrdering(i),
                  _col_ordering.childOrdering(j),
                  entries_per_row
                  )
                );
        }

        //! Returns the subpattern associated with block (i,j).
        SubPattern& subPattern(size_type i, size_type j)
        {
          return _sub_patterns[i * _col_ordering.blockCount() + j];
        }

      private:

        const RowOrdering& _row_ordering;
        const ColOrdering& _col_ordering;
        std::vector<SubPattern> _sub_patterns;

      };


    } // namespace ISTL
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_BCRSPATTERN_HH
