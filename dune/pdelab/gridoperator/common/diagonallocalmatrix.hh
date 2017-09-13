// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDOPERATOR_COMMON_DIAGONALLOCALMATRIX_HH
#define DUNE_PDELAB_GRIDOPERATOR_COMMON_DIAGONALLOCALMATRIX_HH

#include <dune/pdelab/gridoperator/common/localmatrix.hh>

namespace Dune {
  namespace PDELab {

    /**
     * \addtogroup PDELab
     * \{
     */


    //! A dense matrix for storing data associated with the degrees of freedom of a pair of LocalFunctionSpaces.
    /**
     * This container represents a dense matrix based on a std::vector-like storage container and supports accessing
     * its entries indexed by pairs of (LocalFunctionSpace,DOF of LocalFunctionSpace). If the LocalFunctionSpaces
     * contain tags indicating whether they are trial or test spaces, the access methods will also assert that the
     * first space is a test space and the second space is a trial space.
     *
     * \tparam T            The type of values to store in the matrix.
     * \tparam W            The type of weight applied in a WeightedAccumulationView.
     */	template<typename T, typename W = T>
     class DiagonalLocalMatrix
        {
        public:

          //! The type of the underlying storage container.
          /**
           * \warning This is not a matrix-like container anymore, but a std::vector-like one!
           */
          typedef std::vector<T> BaseContainer;

          //! The value type of this container.
          typedef typename BaseContainer::value_type  value_type;

          //! The size type of this container.
          typedef typename BaseContainer::size_type    size_type;

          //! The reference type of this container.
          typedef typename BaseContainer::reference    reference;

          //! The const reference type of this container.
          typedef typename BaseContainer::const_reference const_reference;

          //! The weight type of this container.
          /**
           * A value of this type will be used to assign a weight to contributions in
           * a WeightedAccumulationView.
           */
          typedef W weight_type;

          //! An accumulate-only view of this container that automatically applies a weight to all contributions.
          typedef WeightedMatrixAccumulationView<DiagonalLocalMatrix> WeightedAccumulationView;

          struct iterator
            : public Dune::BidirectionalIteratorFacade<iterator,value_type>
          {

            iterator()
              : _m(nullptr)
              , _i(0)
            {}

            iterator(DiagonalLocalMatrix& m, size_type i)
              : _m(&m)
              , _i(i)
            {}

            bool equals(const iterator& other) const
            {
              return _m == other._m && _i == other._i;
            }

            value_type& dereference() const
            {
              return _m->getEntry(_i,_i);
            }

            void increment()
            {
              ++_i;
            }

            void decrement()
            {
              --_i;
            }

            size_type row() const
            {
              return _i;
            }

            size_type col() const
            {
              return _i;
            }

            DiagonalLocalMatrix* _m;
            size_type _i;

          };

          iterator begin()
          {
            return iterator(*this,0);
          }

          iterator end()
          {
            return iterator(*this,nrows());
          }

          //! Default constructor
          DiagonalLocalMatrix () {}

          //! Construct a LocalMatrix with r rows and c columns.
          DiagonalLocalMatrix (size_type r, size_type c)
            : _container(r)
            , _rows(r)
            , _cols(c)
          {
            assert(r == c);
          }

          //! Construct a LocalMatrix with r rows and c columns and initialize its entries with t.
          DiagonalLocalMatrix (size_type r, size_type c, const T& t)
            : _container(r,t)
            , _rows(r)
            , _cols(c)
          {
            assert(r == c);
          }

          //! Resize the matrix.
          void resize (size_type r, size_type c)
          {
            assert(r == c);
            _container.resize(r);
            _rows = r;
            _cols = c;
          }

          //! Assign t to all entries of the matrix.
          DiagonalLocalMatrix& operator=(const T& t)
          {
            std::fill(_container.begin(),_container.end(),t);
            return *this;
          }

          //! Resize the matrix and assign t to all entries.
          void assign (size_type r, size_type c, const T& t)
          {
            assert(r == c);
            _container.assign(r,t);
            _rows = r;
            _cols = c;
          }

          //! Access the value associated with the i-th DOF of lfsv and the j-th DOF of lfsu.
          /**
           * \param lfsv The LocalFunctionSpace whose DOF will determine the row index of the entry.
           * \param i    The index of the DOF of lfsv that will determine the row index of the entry.
           * \param lfsu The LocalFunctionSpace whose DOF will determine the column index of the entry.
           * \param j    The index of the DOF of lfsu that will determine the column index of the entry.
           */
          template<typename LFSU, typename LFSV>
          T& operator() (const LFSV& lfsv, size_type i, const LFSU& lfsu, size_type j)
          {
            return getEntry(lfsv.localIndex(i),lfsu.localIndex(j));
          }

          //! Access the value associated with the i-th DOF of lfsv and the j-th DOF of lfsu (const version).
          /**
           * \param lfsv The LocalFunctionSpace whose DOF will determine the row index of the entry.
           * \param i    The index of the DOF of lfsv that will determine the row index of the entry.
           * \param lfsu The LocalFunctionSpace whose DOF will determine the column index of the entry.
           * \param j    The index of the DOF of lfsu that will determine the column index of the entry.
           */
          template<typename LFSU, typename LFSV>
          const T& operator() (const LFSV& lfsv, size_type i, const LFSU& lfsu, size_type j) const
          {
            return getEntry(lfsv.localIndex(i),lfsu.localIndex(j));
          }

          //! Multiplies all entries of the matrix with x.
          DiagonalLocalMatrix& operator *= (const T& x)
          {
            using namespace std::placeholders;
            std::transform(
              _container.begin(),
              _container.end(),
              _container.begin(),
              std::bind(std::multiplies<T>(),x,_1)
              );
            return *this;
          }

          //! Returns the number of rows.
          size_type nrows () const
          {
            return _rows;
          }

          //! Returns the number of columns.
          size_type ncols () const
          {
            return _cols;
          }

          //! y += A x
          template<class X, class R>
          void umv (const X& x, R& y) const
          {
            for (size_type i=0; i<_rows; ++i)
              {
                for (size_type j=0; j<_cols; j++)
                  accessBaseContainer(y)[i] += getEntry(i,j) * accessBaseContainer(x)[j];
              }
          }

          //! y += alpha A x
          template<class X, class R>
          void usmv (const value_type& alpha, const X& x, R& y) const
          {
            for (size_type i=0; i<_rows; ++i)
              {
                for (size_type j=0; j<_cols; j++)
                  accessBaseContainer(y)[i] += alpha * getEntry(i,j) * accessBaseContainer(x)[j];
              }
          }

          //! Returns a weighted accumulate-only view of this matrix with the given weight.
          WeightedAccumulationView weightedAccumulationView(weight_type weight)
          {
            return WeightedAccumulationView(*this,weight);
          }

          //! Returns the underlying storage container.
          /**
           * \warning This is not a matrix-like container anymore, but a std::vector-like one!
           */
          BaseContainer& base()
          {
            return _container;
          }

          //! Returns the underlying storage container (const version).
          /**
           * \warning This is not a matrix-like container anymore, but a std::vector-like one!
           */
          const BaseContainer& base() const
          {
            return _container;
          }

          //! Direct (unmapped) access to the (i,j)-th entry of the matrix.
          /**
           * \note This method is just a convenience method for helper function doing things
           *       like printing the matrix, as the BaseContainer is a linear std::vector-like
           *       object and does not know about the shape of the matrix anymore.
           *
           * \warning This method should normally not be called by normal users. Only use it if you
           *          really know what you are doing!
           */
          value_type& getEntry(size_type i, size_type j)
          {
            assert(i == j);
            return _container[i];
          }

          //! Direct (unmapped) access to the (i,j)-th entry of the matrix (const version).
          /**
           * \note This method is just a convenience method for helper function doing things
           *       like printing the matrix, as the BaseContainer is a linear std::vector-like
           *       object and does not know about the shape of the matrix anymore.
           *
           * \warning This method should normally not be called by normal users. Only use it if you
           *          really know what you are doing!
           */
          const value_type& getEntry(size_type i, size_type j) const
          {
            assert(i == j);
            return _container[i];
          }

        private:

          std::vector<T> _container;
          size_type _rows, _cols;
     };

    /**
     * \} group PDELab
     */

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDOPERATOR_COMMON_DIAGONALLOCALMATRIX_HH
