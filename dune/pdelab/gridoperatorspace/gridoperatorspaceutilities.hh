// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_GRIDOPERATORSPACEUTILITIES_HH
#define DUNE_PDELAB_GRIDOPERATORSPACEUTILITIES_HH

#include<dune/common/exceptions.hh>
#include"../common/geometrywrapper.hh"
#include"localmatrix.hh"

namespace Dune {
  namespace PDELab {


	/**@ingroup FlatOperatorSpaceGroup
	   \brief Entry in sparsity pattern

	   The sparsity pattern of a linear operator is described by by connecting
	   degrees of freedom in one element with degrees of freedom in the
	   same element (intra) or an intersecting element (inter). 

	   This numbering is with respect to the depth-first canonical order of the 
	   degrees of freedom of an entity.

	   \nosubgrouping
	*/
	class SparsityLink : public Dune::tuple<int,int>
	{
	public:
	  //! \brief Standard constructor for uninitialized local index
	  SparsityLink ()
	  {}

	  //! \brief Initialize all components
	  SparsityLink (int i, int j)
		: Dune::tuple<int,int>(i,j)
	  {}

	  //! \brief Return first component
	  inline int i () const
	  {
		return Dune::get<0>(*this);
	  } 

	  //! \brief Return second component
	  inline int j () const
	  {
		return Dune::get<1>(*this);
	  } 

	  //! \brief Set both components
	  void set (int i, int j)
	  {
		Dune::get<0>(*this) = i;
		Dune::get<1>(*this) = j;
	  } 
	};

	/**@ingroup FlatOperatorSpaceGroup
	   \brief Layout description for a sparse linear operator

	   \nosubgrouping
	*/
	class LocalSparsityPattern : public std::vector<SparsityLink>
	{};

	//================================================
	// Default matrix backend
	//================================================

	// Simple Backend for std::vector
	class StdVectorFlatMatrixBackend
	{
	public:
	  // Matrix construction
	  template<typename T, typename E>
	  class Matrix : public std::vector<E>
	  {
		typedef std::vector<E> BaseT;
	  public:
		friend class StdVectorFlatMatrixBackend; // for access to line length

		typedef E ElementType;

		Matrix (const T& t) 
		  : n(t.globalSizeU()), 
			BaseT(t.globalSizeU()*t.globalSizeV()) 
		{}

	  private:
		std::vector<int>::size_type line () const
		{
		  return n;
		}

		std::vector<int>::size_type n;
	  };

	  // extract type of matrix element from a Matrix
	  template<class C>
	  struct Value
	  {
		typedef typename C::value_type Type;
	  };

	  //! The size type
	  typedef std::vector<int>::size_type size_type;

	  // get const_reference to container element
	  template<typename C>
	  static const typename C::value_type& const_access (const C& c, size_type i, size_type j)
	  {
		return c.operator[](i*c.line()+j);
	  }

	  // type to store sparsity pattern
	  class Pattern
	  {
	  public:
		void add_link (size_type i, size_type j)
		{
		}
	  };

	  // get non const_reference to container element 
	  // note: this method does not depend on T!
	  template<typename C>
	  static typename C::value_type& access (C& c, size_type i, size_type j)
	  {
		return c.operator[](i*c.line()+j);
	  }

// 	  // read a submatrix given by global indices
// 	  template<typename C, typename RI, typename CI, typename T>
// 	  static void read (const C& c, 
// 						const RI& row_index, const CI& col_index, 
// 						LocalMatrix<T>& submatrix)
// 	  {
// 		submatrix.resize(row_index.size(),col_index.size());
// 		for (int j=0; j<col_index.size(); j++)
// 		  for (int i=0; i<row_index.size(); i++)
// 			submatrix(i,j) = c.operator[](row_index[i]*c.line()+col_index[j]);
// 	  }

// 	  // write a submatrix given by global indices
// 	  template<typename C, typename RI, typename CI, typename T>
// 	  static void write (const RI& row_index, const CI& col_index, 
// 						 const LocalMatrix<T>& submatrix, C& c)
// 	  {
// 		for (int j=0; j<col_index.size(); j++)
// 		  for (int i=0; i<row_index.size(); i++)
// 			c.operator[](row_index[i]*c.line()+col_index[j]) = submatrix(i,j);
// 	  }

// 	  // write a submatrix given by global indices
// 	  template<typename C, typename RI, typename CI, typename T>
// 	  static void add (const RI& row_index, const CI& col_index, 
// 					   const LocalMatrix<T>& submatrix, C& c)
// 	  {
// 		for (int j=0; j<col_index.size(); j++)
// 		  for (int i=0; i<row_index.size(); i++)
// 			c.operator[](row_index[i]*c.line()+col_index[j]) += submatrix(i,j);
// 	  }

	  // clear one row of the matrix
	  template<typename C, typename RI>
	  static void clear_row (RI row, C& c)
	  {
		for (int j=0; j<c.line(); j++)
		  c.operator[](row*c.line()+j) = 0;
	  }
	};



	// compile time switching of function call
    template<typename LA, bool doIt>
    struct LocalAssemblerCallSwitch
    {
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
	  {
	  }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_volume_apply (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
	  {
	  }
      template<typename EG, typename LFSU, typename LFSV>
      static void pattern_volume (const LA& la, const EG& eg, const LFSU& lfsu, const LFSV& lfsv, LocalSparsityPattern& pattern)
      {
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, LocalMatrix<R>& mat)
      {
      }
    };
    template<typename LA>
    struct LocalAssemblerCallSwitch<LA,true>
    {
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
	  {
		la.alpha_volume(eg,lfsu,x,lfsv,r);
	  }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_volume_apply (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
	  {
		la.jacobian_volume_apply(eg,lfsu,x,lfsv,y);
	  }
      template<typename EG, typename LFSU, typename LFSV>
      static void pattern_volume (const LA& la, const EG& eg, const LFSU& lfsu, const LFSV& lfsv, LocalSparsityPattern& pattern)
      {
        la.pattern_volume(eg,lfsu,lfsv,pattern);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, LocalMatrix<R>& mat)
      {
        la.jacobian_volume(eg,lfsu,x,lfsv,mat);
      }
    };


    // derive from this class to add numerical jacobian for volume
    template<typename Imp>
    class NumericalJacobianVolume
    {
    public:

      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, LocalMatrix<R>& mat) const
      {
        const R epsilon=1E-8; // problem: this depends on data type R!
        const int m=lfsv.size();
        const int n=lfsu.size();

        X u(x);
        std::vector<R> down(m,0.0),up(m);

        mat.resize(m,n);
        asImp().alpha_volume(eg,lfsu,u,lfsv,down);	
        for (int j=0; j<n; j++) // loop over columns
          {
            for (int k=0; k<m; k++) up[k]=0.0;
            R delta = epsilon*(1.0+std::abs(u[j]));
            u[j] += delta;
            asImp().alpha_volume(eg,lfsu,u,lfsv,up);
            for (int i=0; i<m; i++)
              mat(i,j) = (up[i]-down[i])/delta;
            u[j] = x[j];
          }
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };

	// derive from this class to add numerical evaluation of jacobian apply
	template<typename Imp>
	class NumericalJacobianVolumeApply
	{
	public:

	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  void jacobian_volume_apply (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
	  {
		typedef typename X::value_type R;
		const R epsilon=1E-8; // problem: this depends on data type R!
		const int m=lfsv.size();
		const int n=lfsu.size();

		X u(x);
		std::vector<R> down(m,0.0),up(m);

		y.resize(m);
		asImp().alpha_volume(eg,lfsu,u,lfsv,down);
		for (int j=0; j<n; j++) // loop over columns
		  {
			for (int k=0; k<m; k++) up[k]=0.0;
			R delta = epsilon*(1.0+std::abs(u[j]));
			u[j] += delta;
			asImp().alpha_volume(eg,lfsu,u,lfsv,up);
			for (int i=0; i<m; i++)
			  y[i] += ((up[i]-down[i])/delta)*x[j];
			u[j] = x[j];
		  }
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
