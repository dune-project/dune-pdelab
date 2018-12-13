// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDOPERATOR_COMMON_LOCALMATRIX_HH
#define DUNE_PDELAB_GRIDOPERATOR_COMMON_LOCALMATRIX_HH

#include <dune/common/iteratorfacades.hh>

#include <dune/pdelab/gridfunctionspace/localvector.hh>

namespace Dune {
  namespace PDELab {

    /**
     * \addtogroup PDELab
     * \{
     */

    //! An accumulate-only view on a local matrix that automatically takes into account an accumulation weight.
    template<typename C>
    class WeightedMatrixAccumulationView
    {

    public:

      //! The type of the underlying container.
      typedef C Container;

      //! The type of the storage container underlying the LocalVector.
      /**
       * \warning This is not a matrix-like container anymore, but a std::vector-like one!
       */
      typedef typename Container::BaseContainer BaseContainer;

      //! The value type of the entries.
      typedef typename C::value_type value_type;

      //! The type of the weight applied when accumulating contributions.
      typedef typename C::weight_type weight_type;

      //! \brief Export this type for uniform handling of the containers
      //!        themselves and their views.
      typedef WeightedMatrixAccumulationView WeightedAccumulationView;

      //! \brief Returns a WeighedAccumulationView with some weight in
      //!        addition to this view's weight
      WeightedAccumulationView weightedAccumulationView(weight_type weight)
      {
        return WeightedAccumulationView(container(),weight*this->weight());
      }

      //! The size_type of the underlying container.
      typedef typename C::size_type size_type;

      //! Returns the weight associated with this view.
      /**
       * \note This can be used together with rawAccumulate() to avoid applying the weight at
       * each loop iteration.
       */
      weight_type weight() const
      {
        return _weight;
      }

      //! Resets the weighting coefficient of the view.
      /**
       * \warning Only call this method when you know what you are doing! It is especially not meant
       *          to be called from inside local operators.
       */
      void setWeight(weight_type weight)
      {
        _weight = weight;
      }

      //! Applies the current weight to v and adds the result to the matrix entry associated with the i-th entry of lfsv and the j-th entry of lfsu.
      template<typename LFSU, typename LFSV>
      void accumulate(const LFSV& lfsv, size_type i,
                      const LFSU& lfsu, size_type j,
                      value_type v)
      {
        _modified = true;
        _container(lfsv,i,lfsu,j) += _weight * v;
      }

      //! Adds v to the (i,j)-th entry of the underlying container without applying the current weight.
      /**
       * \warning When using this method, you must take care of applying the weight yourself, otherwise
       *          your program may exhibit strange behavior or even calculate wrong results!
       */
      template<typename LFSU, typename LFSV>
      void rawAccumulate(const LFSV& lfsv, size_type i,
                         const LFSU& lfsu, size_type j,
                         value_type v)
      {
        _modified = true;
        _container(lfsv,i,lfsu,j) += v;
      }

      //! Returns the number of rows of the underlying container.
      size_type nrows() const
      {
        return _container.nrows();
      }

      //! Returns the number of colums of the underlying container.
      size_type ncols() const
      {
        return _container.ncols();
      }

      //! Returns whether this view has been written to.
      bool modified() const
      {
        return _modified;
      }

      //! Resets the modification state of the view to not modified.
      /**
       * \warning Never call this method from within a local operator, or
       *          your local residual / matrix contributions will be lost!
       */
      void resetModified()
      {
        _modified = false;
      }

      //! Returns the container (of type LocalMatrix) that this view is based on.
      Container& container()
      {
        _modified = true;
        return _container;
      }

      //! Returns the container (of type LocalMatrix) that this view is based on (const version).
      const Container& container() const
      {
        return _container;
      }

      //! Returns the storage container of the underlying LocalMatrix.
      /**
       * \warning This is not a matrix-like container anymore, but a std::vector-like one!
       */
      BaseContainer& base()
      {
        _modified = true;
        return _container.base();
      }

      //! Returns the storage container of the underlying LocalVector (const version).
      /**
       * \warning This is not a matrix-like container anymore, but a std::vector-like one!
       */
      const BaseContainer& base() const
      {
        return _container.base();
      }

      //! Constructor
      WeightedMatrixAccumulationView(Container& container, weight_type weight)
        : _container(container)
        , _weight(weight)
        , _modified(false)
      {}

    private:
      Container& _container;
      weight_type _weight;
      bool _modified;
    };


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
     class LocalMatrix
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
          typedef WeightedMatrixAccumulationView<LocalMatrix> WeightedAccumulationView;

          struct iterator
            : public Dune::BidirectionalIteratorFacade<iterator,value_type>
          {

            iterator()
              : _m(nullptr)
              , _i(0)
              , _j(0)
            {}

            iterator(LocalMatrix& m, size_type i, size_type j)
              : _m(&m)
              , _i(i)
              , _j(j)
            {}

            bool equals(const iterator& other) const
            {
              return _m == other._m && _i == other._i && _j == other._j;
            }

            value_type& dereference() const
            {
              return _m->getEntry(_i,_j);
            }

            void increment()
            {
              if (_j < _m->ncols() - 1)
                ++_j;
              else
                {
                  ++_i;
                  _j = 0;
                }
            }

            void decrement()
            {
              if (_j > 0)
                --_j;
              else
                {
                  --_i;
                  _j = _m->ncols() - 1;
                }
            }

            size_type row() const
            {
              return _i;
            }

            size_type col() const
            {
              return _j;
            }

            LocalMatrix* _m;
            size_type _i;
            size_type _j;

          };

          iterator begin()
          {
            return iterator(*this,0,0);
          }

          iterator end()
          {
            return ncols() > 0 ? iterator(*this, nrows(), 0) : begin();
          }

          //! Default constructor
          LocalMatrix () {}

          //! Construct a LocalMatrix with r rows and c columns.
          LocalMatrix (size_type r, size_type c)
            : _container(r*c)
            , _rows(r)
            , _cols(c)
          {}

          //! Construct a LocalMatrix with r rows and c columns and initialize its entries with t.
          LocalMatrix (size_type r, size_type c, const T& t)
            : _container(r*c,t)
            , _rows(r)
            , _cols(c)
          {}

          //! Resize the matrix.
          void resize (size_type r, size_type c)
          {
            _container.resize(r*c);
            _rows = r;
            _cols = c;
          }

          //! Assign t to all entries of the matrix.
          LocalMatrix& operator=(const T& t)
          {
            std::fill(_container.begin(),_container.end(),t);
            return *this;
          }

          //! Resize the matrix and assign t to all entries.
          void assign (size_type r, size_type c, const T& t)
          {
            _container.assign(r*c,t);
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
          LocalMatrix& operator *= (const T& x)
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
            return _container[j*_rows + i];
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
            return _container[j*_rows + i];
          }

        private:

          std::vector<T> _container;
          size_type _rows, _cols;
     };

    template<class Stream, class T, class W>
    Stream &operator<<(Stream &stream, const LocalMatrix<T,W> &m) {
      for(int r = 0; r < m.nrows(); ++r) {
        if(m.ncols() >= 1)
          stream << m.getEntry(r, 0);
        for(int c = 1; c < m.ncols(); ++c)
          stream << "\t" << m.getEntry(r, c);
        stream << "\n";
      }
      return stream;
    }

    /**
     * \} group PDELab
     */

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDOPERATOR_COMMON_LOCALMATRIX_HH
