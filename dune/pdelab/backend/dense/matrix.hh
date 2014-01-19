// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_DENSE_MATRIX_HH
#define DUNE_PDELAB_BACKEND_DENSE_MATRIX_HH

#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>

#include <dune/common/typetraits.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/pdelab/backend/tags.hh>
#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/dense/descriptors.hh>

namespace Dune {
  namespace PDELab {
    namespace dense {

      template<typename GFSV, typename GFSU, typename C>
      class MatrixContainer
      {

      public:

        typedef C Container;
        typedef typename Container::value_type ElementType;
        typedef ElementType E;

        typedef ElementType field_type;
        typedef typename Container::size_type size_type;

        typedef GFSU TrialGridFunctionSpace;
        typedef GFSV TestGridFunctionSpace;

        typedef typename GFSV::Ordering::Traits::ContainerIndex RowIndex;
        typedef typename GFSU::Ordering::Traits::ContainerIndex ColIndex;

        template<typename RowCache, typename ColCache>
        using LocalView = UncachedMatrixView<MatrixContainer,RowCache,ColCache>;

        template<typename RowCache, typename ColCache>
        using ConstLocalView = ConstUncachedMatrixView<const MatrixContainer,RowCache,ColCache>;


        template<typename GO>
        MatrixContainer(const GO& go)
          : _rows(go.testGridFunctionSpace().size())
          , _cols(go.trialGridFunctionSpace().size())
          , _container(make_shared<Container>(_rows*_cols,E(0)))
        {}

        template<typename GO>
        MatrixContainer(const GO& go, const E& e)
          : _rows(go.testGridFunctionSpace().size())
          , _cols(go.trialGridFunctionSpace().size())
          , _container(make_shared<Container>(_rows*_cols,e))
        {}

        //! Creates an ISTLMatrixContainer without allocating an underlying ISTL matrix.
        explicit MatrixContainer(tags::unattached_container = tags::unattached_container())
          : _rows(0)
          , _cols(0)
        {}

        //! Creates an ISTLMatrixContainer with an empty underlying ISTL matrix.
        explicit MatrixContainer(tags::attached_container)
          : _rows(0)
          , _cols(0)
          , _container(make_shared<Container>())
        {}

        MatrixContainer(const MatrixContainer& rhs)
          : _rows(rhs._rows)
          , _cols(rhs._cols)
          , _container(make_shared<Container>(*(rhs._container)))
        {}

        MatrixContainer& operator=(const MatrixContainer& rhs)
        {
          if (this == &rhs)
            return *this;
          if (_rows == 0 && _cols == 0)
            {
              _rows = rhs._rows;
              _cols = rhs._cols;
            }
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
          return _rows;
        }

        size_type M() const
        {
          return _cols;
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
          return (*_container)[ri[0]*_cols + ci[0]];
        }

        const E& operator()(const RowIndex& ri, const ColIndex& ci) const
        {
          return (*_container)[ri[0]*_cols + ci[0]];
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

        void clear_row(const RowIndex& ri, const E& diagonal_entry)
        {
          std::fill(_container->begin() + ri[0]*_cols,_container->begin() + (ri[0]+1)*_cols,E(0));
          (*this)(ri,ri) = diagonal_entry;
        }

      private:

        std::size_t _rows;
        std::size_t _cols;
        shared_ptr<Container> _container;

      };

    } // namespace dense

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_DENSE_MATRIX_HH
