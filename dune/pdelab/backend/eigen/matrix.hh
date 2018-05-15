// -*- tab-width: 4; indent-tabs-mode: nil -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_EIGEN_MATRIX_HH
#define DUNE_PDELAB_BACKEND_EIGEN_MATRIX_HH

#include <utility>
#include <vector>
#include <set>

#if HAVE_EIGEN

#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/ordering/orderingbase.hh>
#include <Eigen/Sparse>

namespace Dune
{
  namespace PDELab
  {
    namespace Eigen
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
          _matrix.coeffRef(ri.back(),ci.back()) = 0.0;
        }

        M & _matrix;
      };

      /** \brief Wrapper class for a sparse matrix from the Eigen library
       *
       * \tparam GFSV
       * \tparam GFSU
       * \tparam ET Number type used for the matrix entries
       */
      template<typename GFSV, typename GFSU, typename ET, int _Options>
      class MatrixContainer
        : public Backend::impl::Wrapper<::Eigen::SparseMatrix<ET,_Options>>
      {

      public:

        typedef ::Eigen::SparseMatrix<ET,_Options> Container;

      private:

        friend Backend::impl::Wrapper<Container>;

      public:

        typedef ET ElementType;

        typedef ElementType field_type;
        typedef typename Container::Index size_type;
        typedef typename Container::Index index_type;

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

        typedef MatrixPatternInserter<Container> Pattern;

        template<typename GO>
        MatrixContainer(const GO& go)
          : _container(std::make_shared<Container>())
        {
          allocate_matrix(_container, go, ElementType(0));
        }

        template<typename GO>
        MatrixContainer(const GO& go, const ElementType& e)
          : _container(std::make_shared<Container>())
        {
          allocate_matrix(_container, go, e);
        }

        //! Creates an MatrixContainer without allocating an underlying Eigen matrix.
        explicit MatrixContainer(Backend::unattached_container = Backend::unattached_container())
        {}

        //! Creates an MatrixContainer with an empty underlying Eigen matrix.
        explicit MatrixContainer(Backend::attached_container)
        : _container(std::make_shared<Container>())
        {}

        MatrixContainer(const MatrixContainer& rhs)
          : _container(std::make_shared<Container>(*(rhs._container)))
        {}

        MatrixContainer& operator=(const MatrixContainer& rhs)
        {
          if (this == &rhs)
            return *this;
          if (attached())
          {
            base() = rhs.base();
          }
          else
          {
            _container = std::make_shared<Container>(rhs.base());
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
          return _container->rows();
        }

        size_type M() const
        {
          return _container->cols();
        }

        MatrixContainer& operator=(const ElementType& e)
        {
          if(!_container->isCompressed()) _container->makeCompressed();
          for (int i=0; i<_container->nonZeros(); i++)
            _container->valuePtr()[i] = e;
          // we require a sufficiently new Eigen version to use setConstant (newer than in debian testing)
          // _container->setConstant(e);
          return *this;
        }

        MatrixContainer& operator*=(const ElementType& e)
        {
          base() *= e;
          return *this;
        }

        template<typename V>
        void mv(const V& x, V& y) const
        {
          y.base() = base() * x.base();
        }

        template<typename V>
        void usmv(const ElementType alpha, const V& x, V& y) const
        {
          y.base() += alpha * (base() * x.base());
        }

        ElementType& operator()(const RowIndex& ri, const ColIndex& ci)
        {
          return _container->coeffRef(ri[0],ci[0]);
        }

        const ElementType operator()(const RowIndex& ri, const ColIndex& ci) const
        {
          return _container->coeffRef(ri[0],ci[0]);
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
          _container->middleRows(ri[0],1) *= 0.0;
          _container->coeffRef(ri[0],ri[0]) = diagonal_entry;
        }

      protected:
        template<typename GO>
        static void allocate_matrix(std::shared_ptr<Container> & c, const GO & go, const ElementType& e)
        {
          // guess size
          int rows = go.testGridFunctionSpace().ordering().blockCount();
          int cols = go.trialGridFunctionSpace().ordering().blockCount();
          c->resize(rows,cols);
          size_type nz = go.matrixBackend().avg_nz_per_row;
          if (nz)
            c->reserve(::Eigen::VectorXi::Constant(rows,nz));
          // setup pattern
          Pattern pattern(*c);
          go.fill_pattern(pattern);
          // compress matrix
          c->makeCompressed();
        }

        std::shared_ptr< Container > _container;
      };

    } // end namespace EIGEN
  } // namespace PDELab
} // namespace Dune

#endif /* HAVE_EIGEN */

#endif // DUNE_PDELAB_BACKEND_EIGEN_MATRIX_HH
