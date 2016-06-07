// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_FLAT_PATTERN_HH
#define DUNE_PDELAB_BACKEND_ISTL_FLAT_PATTERN_HH

#include <utility>
#include <vector>
#include <algorithm>
#include <set>

#include <dune/common/iteratoradapters.hh>

#include <dune/istl/ellmatrix/host.hh>

#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>

namespace Dune {
  namespace PDELab {
    namespace istl {
      namespace flat {

      template<typename RowOrdering, typename ColOrdering>
      class Pattern
      {

        static const std::size_t empty = ~static_cast<std::size_t>(0);

      public:

        template<typename RI, typename CI>
        void add_link(const RI& ri, const CI& ci)
        {
          auto i = ri.back();
          auto j = ci.back();
          auto start = _indices.begin();
          auto begin = start + _entries_per_row*i;
          auto end = start + _entries_per_row*(i+1);
          auto it = std::find_if(begin,end,[=](const std::size_t k){ return k == j || k == (~0ul);});
          if (it != end)
            {
              *it = j;
            }
          else
            {
              _overflow.insert(std::make_pair(i,j));
            }
        }

        template<typename I>
        void sizes(I rit) const
        {
          auto it = _indices.begin();
          auto end = _indices.begin() + _entries_per_row;
          auto oit = _overflow.begin();
          auto oend = _overflow.end();
          for (std::size_t i = 0; i < _row_ordering.blockCount(); ++i, ++rit, end+=_entries_per_row)
            {
              std::size_t s = 0;
              for (; it != end; ++it)
                {
                  if (*it == empty)
                    break;
                  ++s;
                }
              it = end;
              for (; oit != oend && oit->first == i; ++oit)
                ++s;
              *rit = s;
            }
        }

        std::vector<std::size_t> sizes() const
        {
          std::vector<std::size_t> r(_row_ordering.blockCount());
          sizes(r);
          return std::move(r);
        }

        struct RowStreamer
        {

          template<typename T>
          T streamRow(T out_it)
          {
            T out_begin = out_it;
            for (; _it != _end; ++_it, ++out_it)
              {
                if (*_it == empty)
                  break;
                *out_it = *_it;
              }
            _it = _end;
            for (; _oit != _oend && _oit->first == _row; ++_oit, ++out_it)
              *out_it = _oit->second;
            ++_row;
            _end += _entries_per_row;
            std::sort(out_begin,out_it);
            return out_it;
          }

          template<typename T>
          T streamRow(std::size_t i, T out_it)
          {
            setRow(i);
            return streamRow(out_it);
          }

          std::size_t nextRowIndex() const
          {
            return _row;
          }

          bool operator==(const RowStreamer& r) = delete;
          bool operator!=(const RowStreamer& r) = delete;

          void setRow(std::size_t i)
          {
            if (i == _row)
              return;
            _row = i;
            _it = _begin + i*_entries_per_row;
            _end = _begin + (i+1)*_entries_per_row;
            _oit = _overflow.lower_bound(std::make_pair(i,0));
          }

        private:

          friend class Pattern;

          explicit RowStreamer(const Pattern& p, std::size_t row = 0)
            : _row(row)
            , _entries_per_row(p._entries_per_row)
            , _begin(p._indices.begin())
            , _it(p._indices.begin())
            , _end(p._indices.begin() + p._entries_per_row)
            , _oit(p._overflow.begin())
            , _oend(p._overflow.end())
            , _overflow(p._overflow)
          {}

          std::size_t _row;
          const std::size_t _entries_per_row;
          const std::vector<std::size_t>::const_iterator _begin;
          std::vector<std::size_t>::const_iterator _it;
          std::vector<std::size_t>::const_iterator _end;
          std::set<std::pair<std::size_t,std::size_t> >::const_iterator _oit;
          const std::set<std::pair<std::size_t,std::size_t> >::const_iterator _oend;
          const std::set<std::pair<std::size_t,std::size_t> >& _overflow;

        };

        RowStreamer indexStreamer() const
        {
          return RowStreamer(*this);
        }

        Pattern(const RowOrdering& row_ordering, const ColOrdering& col_ordering, std::size_t entries_per_row)
          : _row_ordering(row_ordering)
          , _col_ordering(col_ordering)
          , _entries_per_row(entries_per_row)
          , _indices(row_ordering.blockCount()*entries_per_row,~0ul)
        {}

        const RowOrdering& rowOrdering() const
        {
          return _row_ordering;
        }

        const ColOrdering& colOrdering() const
        {
          return _row_ordering;
        }

      private:

        friend struct RowStreamer;

        const RowOrdering& _row_ordering;
        const ColOrdering& _col_ordering;
        const std::size_t _entries_per_row;

        std::vector<std::size_t> _indices;
        std::set<std::pair<std::size_t,std::size_t> > _overflow;

      };

      } // namespace flat
    } // namespace istl
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_FLAT_PATTERN_HH
