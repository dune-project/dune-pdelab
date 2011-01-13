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
	  LocalMatrix () {}

	  LocalMatrix (int r, int c)
		: m(r*c), rows(r), cols(c)
	  {}

	  LocalMatrix (int r, int c, const T& t)
		: m(r*c,t), rows(r), cols(c)
	  {}

	  void resize (int r, int c)
	  {
		m.resize(r*c);
		rows = r; 
		cols = c;
	  }

	  void assign (int r, int c, const T& t)
	  {
		m.assign(r*c,t);
		rows = r; 
		cols = c;
	  }

	  const T& operator() (int i, int j) const
	  {
		return m[j*rows+i];
	  }

	  T& operator() (int i, int j)
	  {
		return m[j*rows+i];
	  }

      LocalMatrix& operator *= (const T& x)
      {
        for (size_t i=0; i<m.size(); ++i) m[i] *= x;
        return *this;
      }

	  int nrows () const
	  {
		return rows;
	  }

	  int ncols () const
	  {
		return cols;
	  }

      //! y = A x
      template<class X, class R>
      void umv (const X& x, R& y) const
      {
        for (int i=0; i<rows; ++i)
        {
          for (int j=0; j<cols; j++)
            y[i] += (*this)(i,j) * x[j];
        }
      }
      
	private:
	  std::vector<T> m;
	  int rows, cols;
	};


  } // namespace PDELab
} // namespace Dune

#endif
