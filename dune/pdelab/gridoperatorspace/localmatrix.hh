// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALMATRIX_HH
#define DUNE_PDELAB_LOCALMATRIX_HH

#include<vector>

namespace Dune {
  namespace PDELab {

	// a simple container that stores a dense matrix in a std::vector
	template<typename T>
	class LocalMatrix
	{
	public:
      typedef T value_type;
      typedef int size_type;
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
