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
    namespace EIGEN
    {

      template<typename M>
      struct MatrixPatternInserter
      {
        MatrixPatternInserter(M & mat)
          : _matrix(mat)
        {}

        template<typename RI, typename CI>
        void add_link(const RI& ri, const CI& ci)
        {
          _matrix.coeffRef(ri.back(),ci.back());
        }

        M & _matrix;
      };

      template<typename GFSV, typename GFSU, typename ET, int _Options = Eigen::RowMajor>
      class SparseMatrixContainer
      {

      public:

        typedef Eigen::SparseMatrix<ET,_Options> Container;
        typedef ET ElementType;

        typedef ElementType field_type;
        typedef typename Container::Index size_type;
        typedef typename Container::Index index_type;

        typedef GFSU TrialGridFunctionSpace;
        typedef GFSV TestGridFunctionSpace;

        typedef typename GFSV::Ordering::Traits::ContainerIndex RowIndex;
        typedef typename GFSU::Ordering::Traits::ContainerIndex ColIndex;

        template<typename RowCache, typename ColCache>
        using LocalView = UncachedMatrixView<SparseMatrixContainer,RowCache,ColCache>;

        template<typename RowCache, typename ColCache>
        using ConstLocalView = ConstUncachedMatrixView<const SparseMatrixContainer,RowCache,ColCache>;

        typedef OrderingBase<
          typename GFSV::Ordering::Traits::DOFIndex,
          typename GFSV::Ordering::Traits::ContainerIndex
          > RowOrdering;

        typedef OrderingBase<
          typename GFSU::Ordering::Traits::DOFIndex,
          typename GFSU::Ordering::Traits::ContainerIndex
          > ColOrdering;

        typedef MatrixPatternInserter<Container> Pattern;

        template<typename GO>
        SparseMatrixContainer(const GO& go, int avg_per_row)
          : _container(make_shared<Container>())
        {
          allocate_matrix(_container, go, avg_per_row, ElementType(0));
        }

        template<typename GO>
        SparseMatrixContainer(const GO& go, int avg_per_row, const ElementType& e)
          : _container(make_shared<Container>())
        {
          allocate_matrix(_container, go, avg_per_row, e);
        }

        //! Creates an SparseMatrixContainer without allocating an underlying ISTL matrix.
        explicit SparseMatrixContainer(tags::unattached_container = tags::unattached_container())
        {}

        //! Creates an SparseMatrixContainer with an empty underlying ISTL matrix.
        explicit SparseMatrixContainer(tags::attached_container)
        : _container(make_shared<Container>())
        {}

        SparseMatrixContainer(const SparseMatrixContainer& rhs)
          : _container(make_shared<Container>(*(rhs._container)))
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
            _container = make_shared<Container>(*(rhs._container));
          }
          return *this;
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
          return bool(_container);
        }

        const shared_ptr<Container>& storage() const
        {
          return _container;
        }

        size_type N() const
        {
          return _container->rows();
        }

        size_type M() const
        {
          return _container->cols();
        }

        SparseMatrixContainer& operator=(const ElementType& e)
        {
          if(!_container->isCompressed()) _container->makeCompressed();
          for (int i=0; i<_container->nonZeros(); i++)
            _container->valuePtr()[i] = 1.0;
          // we require a sufficiently new Eigen version to use setConstant (newer than in debian testing)
          // _container->setConstant(e);
          return *this;
        }

        SparseMatrixContainer& operator*=(const ElementType& e)
        {
          (*_container) *= e;
          return *this;
        }

        template<typename V>
        void mv(const V& x, V& y) const
        {
          y = (*_container) * x;
        }

        template<typename V>
        void usmv(const ElementType alpha, const V& x, V& y) const
        {
          y += alpha * (*_container) * x;
        }

        ElementType& operator()(const RowIndex& ri, const ColIndex& ci)
        {
          return _container->coeffRef(ri[0],ci[0]);
        }

        const ElementType operator()(const RowIndex& ri, const ColIndex& ci) const
        {
          return _container->coeff(ri[0],ci[0]);
        }

        const Container& base() const
        {
          return *_container;
        }

        Container& base()
        {
          return *_container;
        }

        void flush()
        {}

        void finalize()
        {}

        void clear_row(const RowIndex& ri, const ElementType& diagonal_entry)
        {
          _container->middleRows(ri[0],1) *= 0.0;
          (*this)(ri,ri) = diagonal_entry;
        }

      protected:
        template<typename GO>
        static void allocate_matrix(shared_ptr<Container> & c, const GO & go, int avg_per_row, const ElementType& e)
        {
          // guess size
          int rows = go.testGridFunctionSpace().ordering().blockCount();
          int cols = go.trialGridFunctionSpace().ordering().blockCount();
          c->resize(rows,cols);
          c->reserve(Eigen::VectorXi::Constant(rows,avg_per_row));
          // setup pattern
          Pattern pattern(*c);
          go.fill_pattern(pattern);
          // compress matrix
          c->makeCompressed();
        }

        shared_ptr< Container > _container;
      };

    } // end namespace EIGEN
  } // namespace PDELab
} // namespace Dune

#endif /* HAVE_EIGEN */

#endif /* DUNE_PDELAB_EIGENMATRIXBACKEND_HH */
// -*- tab-width: 4; indent-tabs-mode: nil -*-
// vi: set et ts=4 sw=2 sts=2:
