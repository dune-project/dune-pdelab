// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALMATRIX_HH
#define DUNE_PDELAB_LOCALMATRIX_HH

#include <dune/pdelab/gridfunctionspace/localvector.hh>

namespace Dune {
  namespace PDELab {

    //! An accumulate-only view on a local matrix that automatically takes into account an accumulation weight.
    template<typename C>
    class WeightedMatrixAccumulationView
    {

    public:

      //! The type of the underlying container.
      typedef C Container;

      //! The value type of the entries.
      typedef typename C::value_type value_type;

      //! The type of the weight applied when accumulating contributions.
      typedef typename C::weight_type weight_type;

      //! The size_type of the underlying container.
      typedef typename C::size_type size_type;

      //! A special wrapper type to enable backwards compatibility with current containers when directly accessing entries.
      typedef WeightedContainerEntryProxy<value_type,weight_type> reference;

      //! Returns the weight associated with this view.
      /**
       * \note This can be used together with rawAccumulate() to avoid applying the weight at
       * each loop iteration.
       */
      weight_type weight()
      {
        _modified = true;
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

      //! Applies the current weight to v and adds the result to the (i,j)-th entry of the container.
      void accumulate(size_type i, size_type j, value_type v)
      {
        _modified = true;
        _container(i,j) += _weight * v;
      }

      //! Adds v to the (i,j)-th entry of the underlying container without applying the current weight.
      /**
       * \warning When using this method, you must take care of applying the weight yourself, otherwise
       *          your program may exhibit strange behavior or even calculate wrong results!
       */
      void rawAccumulate(size_type i, size_type j, value_type v)
      {
        _container(i,j) += v;
      }

      //! Returns a reference proxy to an entry of the underlying container.
      /**
       * \returns A proxy that wraps the entry in the underlying container and
       *          only allows application of operator+=() and operator-=(). In
       *          both cases, the argument will automatically be multiplied by
       *          the weight associated with this view.
       *
       * \deprecated Direct access to the individual container entries is deprecated
       *             and will be removed in a future release. Please use accumulate()
       *             or a combination of weight() and rawAccumulate() to update the
       *             entries of this container.
       */
      reference operator()(size_type i, size_type j) DUNE_DEPRECATED
      {
        _modified = true;
        return reference(_container(i,j),_weight);
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
        return false;
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

      // Constructor
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

	// a simple container that stores a dense matrix in a std::vector
	template<typename T, typename W = T>
	class LocalMatrix
	{
	public:

      typedef T value_type;
      typedef W weight_type;
      typedef int size_type;

      //! An accumulate-only view of the matrix that automatically takes into account an accumulation weight.
      typedef WeightedMatrixAccumulationView<LocalMatrix> WeightedAccumulationView;

	  LocalMatrix () {}

	  LocalMatrix (size_type r, size_type c)
		: m(r*c), rows(r), cols(c)
	  {}

	  LocalMatrix (size_type r, size_type c, const T& t)
		: m(r*c,t), rows(r), cols(c)
	  {}

	  void resize (size_type r, size_type c)
	  {
		m.resize(r*c);
		rows = r;
		cols = c;
	  }

	  void assign (size_type r, size_type c, const T& t)
	  {
		m.assign(r*c,t);
		rows = r;
		cols = c;
	  }

	  const T& operator() (size_type i, size_type j) const
	  {
		return m[j*rows+i];
	  }

	  T& operator() (size_type i, size_type j)
	  {
		return m[j*rows+i];
	  }

      LocalMatrix& operator *= (const T& x)
      {
        for (size_t i=0; i<m.size(); ++i) m[i] *= x;
        return *this;
      }

	  size_type nrows () const
	  {
		return rows;
	  }

	  size_type ncols () const
	  {
		return cols;
	  }

      //! y = A x
      template<class X, class R>
      void umv (const X& x, R& y) const
      {
        for (size_type i=0; i<rows; ++i)
        {
          for (size_type j=0; j<cols; j++)
            y[i] += (*this)(i,j) * x[j];
        }
      }

      //! Returns a weighted accumulate-only view of this matrix with the given weight.
      WeightedAccumulationView weightedAccumulationView(weight_type weight)
      {
        return WeightedAccumulationView(*this,weight);
      }

	private:
	  std::vector<T> m;
	  size_type rows, cols;
	};

    template<class Stream, class T>
    Stream &operator<<(Stream &stream, const LocalMatrix<T> &m) {
      for(int r = 0; r < m.nrows(); ++r) {
        if(m.ncols() >= 1)
          stream << m(r, 0);
        for(int c = 1; c < m.ncols(); ++c)
          stream << "\t" << m(r, c);
        stream << "\n";
      }
      return stream;
    }

  } // namespace PDELab
} // namespace Dune

#endif
