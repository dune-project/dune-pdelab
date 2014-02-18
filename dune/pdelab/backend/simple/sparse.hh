// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_SIMPLE_SPARSE_HH
#define DUNE_PDELAB_BACKEND_SIMPLE_SPARSE_HH

#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>

#include <dune/common/typetraits.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/pdelab/backend/tags.hh>
#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/simple/descriptors.hh>

namespace Dune {
  namespace PDELab {
    namespace simple {

      template<typename RowOrdering, typename ColOrdering>
      class SparseMatrixPattern
        : public std::vector< unordered_map<std::size_t> >
      {

      public:

        typedef std::unordered_map<std::size_t> col_type;

        template<typename RI, typename CI>
        void add_link(const RI& ri, const CI& ci)
        {
          this->resize(_row_ordering.blockCount());
          (*this)[ri.back()].insert(ci.back());
        }

        SparseMatrixPattern(const RowOrdering& row_ordering, const ColOrdering& col_ordering)
          : _row_ordering(row_ordering)
          , _col_ordering(col_ordering)
        {}

      private:

        const RowOrdering& _row_ordering;
        const ColOrdering& _col_ordering;

      };

      template<template<typename> class C, typename ET, typename I>
      struct SparseMatrixData
      {
        typedef ET ElementType;
        typedef ET index_type;
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
         indptr 	CSR format index pointer array of the matrix

         \example
         Consider the following 3x3 matrix
            [1, 0, 2]
            [0, 0, 3]
            [4, 5, 6]
         the data would be stored as
         nnz=6
         data=[1 2 3 4 5 6]
         indices=[0 2 2 0 1 2]
         indptr=[0 2 3 6]
       */
      template<typename GFSV, typename GFSU, template<typename> class C, typename ET, typename I=std::size_t >
      class SparseMatrixContainer
      {

      public:

        typedef SparseMatrixData<C,ET,I> Container;
        typedef ET ElementType;

        typedef ElementType field_type;
        typedef typename Container::size_type size_type;
        typedef I index_type;

        typedef GFSU TrialGridFunctionSpace;
        typedef GFSV TestGridFunctionSpace;

        typedef typename GFSV::Ordering::Traits::ContainerIndex RowIndex;
        typedef typename GFSU::Ordering::Traits::ContainerIndex ColIndex;

        template<typename RowCache, typename ColCache>
        using LocalView = UncachedMatrixView<MatrixContainer,RowCache,ColCache>;

        template<typename RowCache, typename ColCache>
        using ConstLocalView = ConstUncachedMatrixView<const MatrixContainer,RowCache,ColCache>;

        typedef OrderingBase<
          typename GFSV::Ordering::Traits::DOFIndex,
          typename GFSV::Ordering::Traits::ContainerIndex
          > RowOrdering;

        typedef OrderingBase<
          typename GFSU::Ordering::Traits::DOFIndex,
          typename GFSU::Ordering::Traits::ContainerIndex
          > ColOrdering;

        typedef Pattern<RowOrdering,ColOrdering> Pattern;

        template<typename GO>
        SparseMatrixContainer(const GO& go)
          : _container(make_shared<Container>())
        {
          allocate_matrix(_container, go, E(0));
        }

        template<typename GO>
        SparseMatrixContainer(const GO& go, const E& e)
          : _container(make_shared<Container>())
        {
          allocate_matrix(_container, go, e);
        }

        //! Creates an ISTLMatrixContainer without allocating an underlying ISTL matrix.
        explicit MatrixContainer(tags::unattached_container = tags::unattached_container())
        {}

        //! Creates an ISTLMatrixContainer with an empty underlying ISTL matrix.
        explicit MatrixContainer(tags::attached_container)
        : _container(make_shared<Container>())
        {}

        MatrixContainer(const MatrixContainer& rhs)
          : _container(make_shared<Container>(*(rhs._container)))
        {
          _container->_rows = 0;
          _container->_cols = 0;
          _container->_non_zeros = 0;
        }

        MatrixContainer& operator=(const MatrixContainer& rhs)
        {
          if (this == &rhs)
            return *this;
          if (attached())
          {
            (*_container) = (*(rhs._container));
          }
          else
          {
            _container = make_shared<Container>(*(rhs._container));
          }
        }

        void detach()
        {
          _container.reset();
        }

        void attach(shared_ptr<Container> container)
        {
          _container = container;
        }

        bool attached() const
        {
          return bool(_data);
        }

        const shared_ptr<Container>& storage() const
        {
          return _container;
        }

        size_type N() const
        {
          return _container->_cols;
        }

        size_type M() const
        {
          return _container->_cols;
        }

        MatrixContainer& operator=(const E& e)
        {
          std::fill(_container->begin(),_container->end(),e);
          return *this;
        }

        MatrixContainer& operator*=(const E& e)
        {
          using namespace std::placeholders;
          std::transform(_container->begin(),_container->end(),_container->begin(),std::bind(std::multiplies<E>(),e,_1));
          return *this;
        }

        template<typename V>
        void mv(const V& x, V& y) const
        {
          auto rowit = _container->begin();
          for (auto& v : y)
          {
            v = std::inner_product(rowit,rowit + _cols,x.begin(),E(0));
            rowit += _cols;
          }
        }

        template<typename V>
        void usmv(const E alpha, const V& x, V& y) const
        {
          auto rowit = _container->begin();
          for (auto& v : y)
          {
            v += alpha * std::inner_product(rowit,rowit + _cols,x.begin(),E(0));
            rowit += _cols;
          }
        }

        E& operator()(const RowIndex& ri, const ColIndex& ci)
        {
          // entries are in ascending order
          auto begin = _container->_colindex.begin() + _container->_rowoffset[ri[0]];
          auto end = _container->_colindex.begin() + _container->_rowoffset[ri[0]+1];
          auto it = std::lower_bound(begin,end,ci[0]);
          assert(it != _container->_colindex.end());
          return _container->data[it - _container->_colindex.begin()];
        }

        const E& operator()(const RowIndex& ri, const ColIndex& ci) const
        {
          // entries are in ascending order
          auto begin = _container->_colindex.begin() + _container->_rowoffset[ri[0]];
          auto end = _container->_colindex.begin() + _container->_rowoffset[ri[0]+1];
          auto it = std::lower_bound(begin,end,ci[0]);
          assert(it != _container->_colindex.end());
          return _container->data[it - _container->_colindex.begin()];
        }

        const Container& base() const
        {
          return *_container;
        }

        Container& base()
        {
          return _container;
        }

        void flush()
        {}

        void finalize()
        {}

        void clear_row(const RowIndex& ri, const E& diagonal_entry)
        {
          std::fill(_container->_data.begin() + _container->_rowoffset[ri[0]], _container->begin() + _container->_rowoffset[ri[0]+1], E(0));
          (*this)(ri,ri) = diagonal_entry;
        }

      protected:
        template<typename GO>
        void allocate_matrix(shared_ptr<Container> c, const GO & go, const E & e)
        {
          typedef typename Pattern::col_type col_type;
          Pattern pattern(go.testGridFunctionSpace().ordering(),go.trialGridFunctionSpace().ordering());
          go.fill_pattern(pattern);

          c->_rows = go.testGridFunctionSpace().size();
          c->_cols = go.trialGridFunctionSpace().size();
          // compute row offsets
          c->_rowoffset.resize(c->_rows+1);
          std::partial_sum(
            pattern.begin(), pattern.end(), c->_rowoffset.begin()+1,
            [](size_t x, const col_type & entry) -> size_t { return x+entry.size();}
            );
          // compute non-zeros
          c->_non_zeros = c->_rowoffsets.back();
          // allocate col/data vectors
          c->_data.resize(c->_non_zeros, e);
          c->_colindex.resize(c->_non_zeros);
          // copy pattern
          auto colit = c->_colindex.begin();
          c->_rowoffset[0] = 0;
          for (size_t r = 0; r < pattern.size(); r++)
          {
            std::copy(pattern[r].begin(),pattern[r].end(),colit);
          }
        }

        shared_ptr< Container > _container;
      };

    } // namespace simple

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_SIMPLE_SPARSE_HH
