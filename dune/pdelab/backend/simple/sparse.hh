// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_SIMPLE_SPARSE_HH
#define DUNE_PDELAB_BACKEND_SIMPLE_SPARSE_HH

#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <memory>
#include <unordered_set>

#include <dune/common/typetraits.hh>
#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/simple/descriptors.hh>

namespace Dune {
  namespace PDELab {
    namespace Simple {

      class SparseMatrixPattern
        : public std::vector< std::unordered_set<std::size_t> >
      {

      public:

        typedef std::unordered_set<std::size_t> col_type;

        template<typename RI, typename CI>
        void add_link(const RI& ri, const CI& ci)
        {
          (*this)[ri.back()].insert(ci.back());
        }

        SparseMatrixPattern(std::size_t rows)
        : std::vector< std::unordered_set<std::size_t> >(rows)
        {}

      };

      template<template<typename> class C, typename ET, typename I>
      struct SparseMatrixData
      {
        typedef ET ElementType;
        typedef I  index_type;
        typedef std::size_t size_type;
        std::size_t _rows;
        std::size_t _cols;
        std::size_t _non_zeros;
        C<ElementType> _data;
        C<index_type>  _colindex;
        C<index_type>  _rowoffset;
      };

      /**
         \brief Simple backend for CSR matrices

         Matrix stored as a:
         nnz 		Number of nonzero elements
         data 		CSR format data array of the matrix
         indices 	CSR format index array of the matrix
         rowptr 	CSR format index pointer array of the matrix

         \example
         Consider the following 3x3 matrix
            [1, 0, 2]
            [0, 0, 3]
            [4, 5, 6]
         the data would be stored as
         nnz=6
         data=[1 2 3 4 5 6]
         indices=[0 2 2 0 1 2]
         rowptr=[0 2 3 6]
       */
      template<typename GFSV, typename GFSU, template<typename> class C, typename ET, typename I>
      class SparseMatrixContainer
        : public Backend::impl::Wrapper<SparseMatrixData<C,ET,I> >
      {

      public:

        typedef SparseMatrixData<C,ET,I> Container;

      private:

        friend Backend::impl::Wrapper<Container>;

      public:

        typedef ET ElementType;

        typedef ElementType field_type;
        typedef typename Container::size_type size_type;
        typedef I index_type;

        typedef GFSU TrialGridFunctionSpace;
        typedef GFSV TestGridFunctionSpace;

        typedef typename GFSV::Ordering::Traits::ContainerIndex RowIndex;
        typedef typename GFSU::Ordering::Traits::ContainerIndex ColIndex;

        template<typename RowCache, typename ColCache>
        using LocalView = UncachedMatrixView<SparseMatrixContainer,RowCache,ColCache>;

        template<typename RowCache, typename ColCache>
        using ConstLocalView = ConstUncachedMatrixView<const SparseMatrixContainer,RowCache,ColCache>;

        typedef SparseMatrixPattern Pattern;

        template<typename GO>
        SparseMatrixContainer(const GO& go)
          : _container(std::make_shared<Container>())
        {
          allocate_matrix(_container, go, ElementType(0));
        }

        template<typename GO>
        SparseMatrixContainer(const GO& go, const ElementType& e)
          : _container(std::make_shared<Container>())
        {
          allocate_matrix(_container, go, e);
        }

        //! Creates an SparseMatrixContainer without allocating an underlying ISTL matrix.
        explicit SparseMatrixContainer(Backend::unattached_container = Backend::unattached_container())
        {}

        //! Creates an SparseMatrixContainer with an empty underlying ISTL matrix.
        explicit SparseMatrixContainer(Backend::attached_container)
        : _container(std::make_shared<Container>())
        {}

        SparseMatrixContainer(const SparseMatrixContainer& rhs)
          : _container(std::make_shared<Container>(*(rhs._container)))
        {}

        SparseMatrixContainer& operator=(const SparseMatrixContainer& rhs)
        {
          if (this == &rhs)
            return *this;
          if (attached())
          {
            (*_container) = (*(rhs._container));
          }
          else
          {
            _container = std::make_shared<Container>(*(rhs._container));
          }
          return *this;
        }

        void detach()
        {
          _container.reset();
        }

        void attach(std::shared_ptr<Container> container)
        {
          _container = container;
        }

        bool attached() const
        {
          return bool(_container);
        }

        const std::shared_ptr<Container>& storage() const
        {
          return _container;
        }

        size_type N() const
        {
          return _container->_rows;
        }

        size_type M() const
        {
          return _container->_cols;
        }

        SparseMatrixContainer& operator=(const ElementType& e)
        {
          std::fill(_container->_data.begin(),_container->_data.end(),e);
          return *this;
        }

        SparseMatrixContainer& operator*=(const ElementType& e)
        {
          using namespace std::placeholders;
          std::transform(_container->_data.begin(),_container->_data.end(),_container->_data.begin(),std::bind(std::multiplies<ET>(),e,_1));
          return *this;
        }

        template<typename V>
        void mv(const V& x, V& y) const
        {
          assert(y.N() == N());
          assert(x.N() == M());
          for (std::size_t r = 0; r < N(); ++r)
          {
            y.base()[r] = sparse_inner_product(r,x);
          }
        }

        template<typename V>
        void usmv(const ElementType alpha, const V& x, V& y) const
        {
          assert(y.N() == N());
          assert(x.N() == M());
          for (std::size_t r = 0; r < N(); ++r)
          {
            y.base()[r] += alpha * sparse_inner_product(r,x);
          }
        }

        ElementType& operator()(const RowIndex& ri, const ColIndex& ci)
        {
          // entries are in ascending order
          auto begin = _container->_colindex.begin() + _container->_rowoffset[ri[0]];
          auto end = _container->_colindex.begin() + _container->_rowoffset[ri[0]+1];
          auto it = std::lower_bound(begin,end,ci[0]);
          assert (it != end);
          return _container->_data[it - _container->_colindex.begin()];
        }

        const ElementType& operator()(const RowIndex& ri, const ColIndex& ci) const
        {
          // entries are in ascending order
          auto begin = _container->_colindex.begin() + _container->_rowoffset[ri[0]];
          auto end = _container->_colindex.begin() + _container->_rowoffset[ri[0]+1];
          auto it = std::lower_bound(begin,end,ci[0]);
          assert(it != end);
          return _container->_data[it - _container->_colindex.begin()];
        }

        const Container& base() const
        {
          return *_container;
        }

        Container& base()
        {
          return *_container;
        }

      private:

        const Container& native() const
        {
          return *_container;
        }

        Container& native()
        {
          return *_container;
        }

      public:

        void flush()
        {}

        void finalize()
        {}

        void clear_row(const RowIndex& ri, const ElementType& diagonal_entry)
        {
          std::fill(
            _container->_data.begin() + _container->_rowoffset[ri[0]],
            _container->_data.begin() + _container->_rowoffset[ri[0]+1], ElementType(0));
          (*this)(ri,ri) = diagonal_entry;
        }

        void clear_row_block(const RowIndex& ri, const ElementType& diagonal_entry)
        {
          std::fill(
            _container->_data.begin() + _container->_rowoffset[ri[0]],
            _container->_data.begin() + _container->_rowoffset[ri[0]+1], ElementType(0));
          (*this)(ri,ri) = diagonal_entry;
        }

      protected:
        template<typename GO>
        static void allocate_matrix(std::shared_ptr<Container> & c, const GO & go, const ElementType& e)
        {
          typedef typename Pattern::col_type col_type;
          Pattern pattern(go.testGridFunctionSpace().ordering().blockCount());
          go.fill_pattern(pattern);

          c->_rows = go.testGridFunctionSpace().size();
          c->_cols = go.trialGridFunctionSpace().size();
          // compute row offsets
          c->_rowoffset.resize(c->_rows+1);
          size_type offset = 0;
          auto calc_offset = [=](const col_type & entry) mutable -> size_t { offset += entry.size(); return offset; };
          std::transform(pattern.begin(), pattern.end(),
            c->_rowoffset.begin()+1,
            calc_offset);
          // compute non-zeros
          c->_non_zeros = c->_rowoffset.back();
          // allocate col/data vectors
          c->_data.resize(c->_non_zeros, e);
          c->_colindex.resize(c->_non_zeros);
          // copy pattern
          auto colit = c->_colindex.begin();
          c->_rowoffset[0] = 0;
          for (auto & row : pattern)
          {
            auto last = std::copy(row.begin(),row.end(),colit);
            std::sort(colit, last);
            colit = last;
          }
        }

        template<typename V>
        ElementType sparse_inner_product (std::size_t row, const V & x) const {
            std::size_t begin = _container->_rowoffset[row];
            std::size_t end = _container->_rowoffset[row+1];
            auto accu = [](const ElementType & a, const ElementType & b) -> ElementType { return a+b; };
            auto mult = [=](const ElementType & a, const index_type & i) -> ElementType { return a * x.base()[i]; };
            return std::inner_product(
              &_container->_data[begin],
              &_container->_data[end],
              &_container->_colindex[begin],
              ElementType(0),
              accu, mult
              );
        }

        std::shared_ptr< Container > _container;
      };

    } // namespace Simple
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_SIMPLE_SPARSE_HH
