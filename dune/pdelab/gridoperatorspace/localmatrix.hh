// -*- tab-width: 4; indent-tabs-mode: nil -*-
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

	  const T& operator() (int i, int j) const
	  {
		return m[j*rows+i];
	  }

	  T& operator() (int i, int j)
	  {
		return m[j*rows+i];
	  }

	  int nrows () const
	  {
		return rows;
	  }

	  int ncols () const
	  {
		return cols;
	  }

	private:
	  std::vector<T> m;
	  int rows, cols;
	};


  } // namespace PDELab
} // namespace Dune

#endif
