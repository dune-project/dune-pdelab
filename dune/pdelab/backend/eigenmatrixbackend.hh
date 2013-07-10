#ifndef DUNE_PDELAB_EIGENMATRIXBACKEND_HH
#define DUNE_PDELAB_EIGENMATRIXBACKEND_HH

#include <utility>
#include <vector>
#include <set>

#if HAVE_EIGEN

#include <Eigen/Sparse>

namespace Dune
{
  namespace PDELab
  {

    //! Eigen backend for LinearOperatorSpace
    class SparseEigenMatrixBackend
    {
    public:

      //! container construction
      template< typename E = double >
        class Matrix: public Eigen::SparseMatrix< E, Eigen::RowMajor >
      {
      public:
        typedef Eigen::SparseMatrix< E, Eigen::RowMajor > EigenMatrixType;
        typedef typename EigenMatrixType::Index size_type;
        typedef E ElementType;
        typedef E field_type;
        typedef EigenMatrixType BaseT;
        typedef SparseEigenMatrixBackend Backend;

        //! construct container
        template<typename T>
        Matrix (const T& t)
          : BaseT(t.globalSizeV(), t.globalSizeU())
        {
          Pattern pattern(t.globalSizeV(), t.globalSizeU());
          t.fill_pattern(pattern);

          //      this->reserve(Eigen::VectorXi::Constant(this->outerSize(), pattern.max_nnzs_));

          for (size_type row = 0; (unsigned) row < this->outerSize(); ++row)
          {
            this->startVec(row);
            for (typename std::set< size_type >::iterator it = pattern[row].begin();
                 it != pattern[row].end(); ++it)
            {
              size_type column = *it;
              //                    std::cout << "insert row=" << row << " column=" << column << std::endl;
              this->insertBackByOuterInner(row, column);
              //                    this->insert(row, column);
            }
          }
          this->finalize();
          this->makeCompressed();
        }

        //! set from element
        Matrix& operator= (const E& x)
        {
          if (x == 0.0)
            this->setZero();
          else
          {
            for (int k = 0; k < this->outerSize(); ++k)
            {
              for (typename BaseT::InnerIterator it((*this), k); it; ++it)
              {
                BaseT & mat = static_cast< BaseT & >(*this);
                mat.coeffRef(it.row(), it.col()) = x;
              }
            }
          }
          return *this;
        }

        //! for debugging and RB Access
        BaseT& base ()
        {
          return *this;
        }

        //! for debugging and RB access
        const BaseT& base () const
        {
          return *this;
        }
      };

      //! The size type
      typedef Eigen::SparseMatrix<double, Eigen::RowMajor>::Index size_type;

      //! extract type of container element
      template< class C >
        struct Value
        {
          typedef typename C::ElementType Type;
        };

      //! type to store sparsity pattern of block indices
      class Pattern: public std::vector< std::set< size_type > >
      {
        typedef std::vector< std::set< size_type > > BaseT;
      public:
        Pattern (size_type m_, size_type n_)
          : max_nnzs_(5)
        {
          this->resize(m_);
        }

        void add_link (size_type i, size_type j)
        {
          (*this)[i].insert(j);
          size_t tmp_nnzs = (*this)[i].size();
          if (tmp_nnzs > max_nnzs_) max_nnzs_ = tmp_nnzs;
        }
      public:
        size_t max_nnzs_;
      };

      //==========================
      // individual element access
      //==========================

      template<typename LFSV, typename LFSU, typename E>
       struct Accessor
       {
         typedef E ElementType;
         typedef SparseEigenMatrixBackend Backend;
         typedef typename Backend::size_type size_type;
         typedef typename Backend::template Matrix<ElementType> Matrix;

         Accessor(Matrix& matrix, const LFSV& lfsv, const LFSU& lfsu)
           : _matrix(matrix)
           , _lfsv(lfsv)
           , _lfsu(lfsu)
         {}

         void set(size_type i, size_type j, const typename Matrix::ElementType& v)
         {
           Backend::access(_matrix,_lfsv.globalIndex(i),_lfsu.globalIndex(j)) = v;
         }

         void add(size_type i, size_type j, const typename Matrix::ElementType& v)
         {
           Backend::access(_matrix,_lfsv.globalIndex(i),_lfsu.globalIndex(j)) += v;
         }

         typename Matrix::ElementType get(size_type i, size_type j) const
         {
           return Backend::access(_matrix,_lfsv.globalIndex(i),_lfsu.globalIndex(j));
         }


         void setGlobal(size_type gi, size_type gj, const typename Matrix::ElementType& v)
         {
           Backend::access(_matrix,gi,gj) = v;
         }

         void addGlobal(size_type gi, size_type gj, const typename Matrix::ElementType& v)
         {
           Backend::access(_matrix,gi,gj) += v;
         }

         typename Matrix::field_type getGlobal(size_type gi, size_type gj) const
         {
           return Backend::access(_matrix,gi,gj);
         }

         Matrix& _matrix;
         const LFSV& _lfsv;
         const LFSU& _lfsu;

       };

      //! get const_reference to container element
      /**
       * \note this method does not depend on T!
       */
      template< typename C >
        static const typename C::ElementType& access (const C& c, size_type row, size_type column)
        {
          return c.coeff(row, column);
        }

      //! get non const_reference to container element
      /**
       * \note this method does not depend on T!
       */
      template< typename C >
        static typename C::ElementType& access (C& c, size_type row, size_type column)
        {
          return c.coeffRef(row, column);
        }

      //! clear one row of the matrix
      template< typename C, typename RI >
        static void clear_row (RI row, C& c, const typename C::ElementType& diag_val)
        {
          typename C::BaseT::InnerIterator it(c, row);
          for (; it; ++it)
          {
            c.coeffRef(it.row(), it.col()) = 0.0;
          }
          access(c, row, row) = diag_val;
        }

      template<typename C>
      static void flush(C& matrix)
      {}

      template<typename C>
      static void finalize(C& matrix)
      {}

    };

  } // namespace PDELab
} // namespace Dune

#endif /* HAVE_EIGEN */

#endif /* DUNE_PDELAB_EIGENMATRIXBACKEND_HH */
// -*- tab-width: 4; indent-tabs-mode: nil -*-
// vi: set et ts=4 sw=2 sts=2:
