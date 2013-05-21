// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_BLOCKMATRIXDIAGONAL_HH
#define DUNE_PDELAB_BACKEND_ISTL_BLOCKMATRIXDIAGONAL_HH

#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istl/utility.hh>
#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>

namespace Dune {
  namespace PDELab {
    namespace istl {

#ifndef DOXYGEN

      // implementation namespace for matrix diagonal

      namespace diagonal {

        // TMP for determining the type of vector to use for the matrix diagonal
        template<typename M>
        struct matrix_element_vector;

        // At FieldMatrix level, we keep the whole matrix
        template<typename E, int n, int m>
        struct matrix_element_vector<
          FieldMatrix<E,n,m>
          >
        {
          typedef FieldMatrix<E,n,m> type;
        };

        // At BCRSMatrix level, we use a BlockVector and recursively apply the
        // TMP to the block type.
        template<typename Block, typename Allocator>
        struct matrix_element_vector<
          BCRSMatrix<Block,Allocator>
          >
        {
          typedef BlockVector<
            typename matrix_element_vector<Block>::type,
            Allocator
            > type;
        };


        // Function for extracting the diagonal from a matrix.
        // For the FieldMatrix, we just copy the complete matrix.
        template<typename FieldMatrix>
        void matrix_element_vector_from_matrix(tags::field_matrix, FieldMatrix& c, const FieldMatrix& matrix)
        {
          c = matrix;
        }

        // For the BCRSMatrix, we recursively copy diagonal blocks.
        template<typename BlockVector, typename BCRSMatrix>
        void matrix_element_vector_from_matrix(tags::block_vector, BlockVector& c, const BCRSMatrix& m)
        {
          const std::size_t rows = m.N();
          c.resize(rows,false);
          for (std::size_t i = 0; i < rows; ++i)
            matrix_element_vector_from_matrix(container_tag(c[i]),c[i],m[i][i]);
        }


        // Function for inverting the diagonal.
        // The FieldMatrix supports direct inverson.
        template<typename FieldMatrix>
        void invert_blocks(tags::field_matrix, FieldMatrix& c)
        {
          c.invert();
        }

        // For a BCRSMatrix, we recursively invert all diagonal entries.
        template<typename BlockVector>
        void invert_blocks(tags::block_vector, BlockVector& c)
        {
          const std::size_t rows = c.size();
          for (std::size_t i = 0; i < rows; ++i)
            invert_blocks(container_tag(c[i]),c[i]);
        }


        // Matrix-vector product between matrix consisting of only the
        // diagonal and a matching vector.
        // For the FieldMatrix, simply call its matrix-vector product.
        template<typename FieldMatrix, typename X, typename Y>
        void mv(tags::field_matrix, const FieldMatrix& c, const X& x, Y& y)
        {
          c.mv(x,y);
        }

        // For the BCRSMatrix, recursively apply this function to the
        // individual blocks.
        template<typename BlockVector, typename X, typename Y>
        void mv(tags::block_vector, const BlockVector& c, const X& x, Y& y)
        {
          const std::size_t rows = c.size();
          for (std::size_t i = 0; i < rows; ++i)
            mv(container_tag(c[i]),c[i],x[i],y[i]);
        }


        template<typename FieldMatrix, typename CI>
        std::size_t row_size(tags::field_matrix, const FieldMatrix& c, const CI& ci, int i)
        {
          return FieldMatrix::cols;
        }

        template<typename FieldMatrix>
        std::size_t row_size(tags::field_matrix, const FieldMatrix& c)
        {
          return FieldMatrix::cols;
        }

        template<typename BlockVector, typename CI>
        std::size_t row_size(tags::block_vector, const BlockVector& c, const CI& ci, int i)
        {
          return row_size(container_tag(c[0]),c[0]);
        }

        template<typename BlockVector>
        std::size_t row_size(tags::block_vector, const BlockVector& c)
        {
          return row_size(container_tag(c[0]),c[0]);
        }

        template<typename FieldMatrix, typename CI>
        typename FieldMatrix::field_type* row_begin(tags::field_matrix_1_any, FieldMatrix& c, const CI& ci, int i)
        {
          assert(i == -1);
          return &(*c[0].begin());
        }

        template<typename FieldMatrix, typename CI>
        typename FieldMatrix::field_type* row_begin(tags::field_matrix_n_any, FieldMatrix& c, const CI& ci, int i)
        {
          assert(i == 0);
          return &(*c[ci[0]].begin());
        }

        template<typename FieldMatrix, typename CI>
        typename FieldMatrix::field_type* row_end(tags::field_matrix_1_any, FieldMatrix& c, const CI& ci, int i)
        {
          assert(i == -1);
          return &(*c[0].end());
        }

        template<typename FieldMatrix, typename CI>
        typename FieldMatrix::field_type* row_end(tags::field_matrix_n_any, FieldMatrix& c, const CI& ci, int i)
        {
          assert(i == 0);
          return &(*c[ci[0]].end());
        }


        template<typename BlockVector, typename CI>
        typename BlockVector::field_type* row_begin(tags::block_vector, BlockVector& c, const CI& ci, std::size_t i)
        {
          return row_begin(container_tag(c[ci[i]]),c[ci[i]],ci,i-1);
        }

        template<typename BlockVector, typename CI>
        typename BlockVector::field_type* row_end(tags::block_vector, BlockVector& c, const CI& ci, std::size_t i)
        {
          return row_end(container_tag(c[ci[i]]),c[ci[i]],ci,i-1);
        }


      } // namespace diagonal

#endif // DOXYGEN


      template<typename M>
      struct BlockMatrixDiagonal
      {

        typedef typename raw_type<M>::type Matrix;

        struct MatrixElementVector
        {

          typedef typename diagonal::matrix_element_vector<Matrix>::type Container;
          typedef typename Container::field_type field_type;
          typedef field_type* iterator;

          Container _container;

          MatrixElementVector(const M& m)
          {
            diagonal::matrix_element_vector_from_matrix(container_tag(_container),_container,raw(m));
          }

          void invert()
          {
            diagonal::invert_blocks(container_tag(_container),_container);
          }

          template<typename X, typename Y>
          void mv(const X& x, Y& y) const
          {
            diagonal::mv(container_tag(_container),_container,raw(x),raw(y));
          }

          template<typename ContainerIndex>
          std::size_t row_size(const ContainerIndex& ci) const
          {
            return diagonal::row_size(container_tag(_container),_container,ci,ci.size()-1);
          }

          template<typename ContainerIndex>
          iterator row_begin(const ContainerIndex& ci)
          {
            return diagonal::row_begin(container_tag(_container),_container,ci,ci.size()-1);
          }

          template<typename ContainerIndex>
          iterator row_end(const ContainerIndex& ci)
          {
            return diagonal::row_end(container_tag(_container),_container,ci,ci.size()-1);
          }

        };


        template<typename GFS>
        class AddMatrixElementVectorDataHandle
          : public Dune::CommDataHandleIF<AddMatrixElementVectorDataHandle<GFS>,typename Matrix::field_type>
        {

        public:

          typedef typename Matrix::field_type DataType;
          typedef typename GFS::Traits::SizeType size_type;

          AddMatrixElementVectorDataHandle(const GFS& gfs, MatrixElementVector& v)
            : _gfs(gfs)
            , _index_cache(gfs)
            , _v(v)
          {}

          //! returns true if data for this codim should be communicated
          bool contains(int dim, int codim) const
          {
            return _gfs.dataHandleContains(codim);
          }

          //!  \brief returns true if size per entity of given dim and codim is a constant
          bool fixedsize(int dim, int codim) const
          {
            return _gfs.dataHandleFixedSize(codim);
          }

          /*!  \brief how many objects of type DataType have to be sent for a given entity

            Note: Only the sender side needs to know this size.
          */
          template<typename Entity>
          size_type size(Entity& e) const
          {
            _index_cache.update(e);

            size_type s = 0;
            for (size_type i = 0; i < _index_cache.size(); ++i)
              s += _v.row_size(_index_cache.containerIndex(i));
            return s;
          }

          //! \brief pack data from user to message buffer
          template<typename MessageBuffer, typename Entity>
          void gather(MessageBuffer& buff, const Entity& e) const
          {
            _index_cache.update(e);
            for (size_type i = 0; i < _index_cache.size(); ++i)
              {
                const CI& ci = _index_cache.containerIndex(i);
                for (RowIterator it = _v.row_begin(ci),
                       end_it = _v.row_end(ci);
                     it != end_it;
                     ++it)
                  buff.write(*it);
              }
          }

          /*! \brief unpack data from message buffer to user

            n is the number of objects sent by the sender
          */
          template<typename MessageBuffer, typename Entity>
          void scatter(MessageBuffer& buff, const Entity& e, size_type n)
          {
            _index_cache.update(e);
            for (size_type i = 0; i < _index_cache.size(); ++i)
              {
                const CI& ci = _index_cache.containerIndex(i);
                for (RowIterator it = _v.row_begin(ci),
                       end_it = _v.row_end(ci);
                     it != end_it;
                     ++it)
                  {
                    DataType x;
                    buff.read(x);
                    *it += x;
                  }
              }
          }

        private:

          typedef EntityIndexCache<GFS> IndexCache;
          typedef typename IndexCache::ContainerIndex CI;
          typedef typename MatrixElementVector::iterator RowIterator;

          const GFS& _gfs;
          mutable IndexCache _index_cache;
          MatrixElementVector& _v;

        };

      };

    } // namespace istl
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_BLOCKMATRIXDIAGONAL_HH
