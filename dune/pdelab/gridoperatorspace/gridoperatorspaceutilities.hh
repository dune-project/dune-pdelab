// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATORSPACEUTILITIES_HH
#define DUNE_PDELAB_GRIDOPERATORSPACEUTILITIES_HH

#include <vector>

#include <dune/common/tuples.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>

namespace Dune {
  namespace PDELab {

	//! collect types exported by a grid operator space
	template<typename GFSU, typename GFSV, typename B, 
			 typename CU, typename CV>
	struct GridOperatorSpaceTraits
	{
	  typedef GFSU TrialGridFunctionSpace;

	  typedef CU TrialConstraintsType;

	  typedef GFSV TestGridFunctionSpace;

	  typedef CV TestConstraintsType;

	  //! \brief the grid view where grid function is defined upon
	  typedef typename GFSU::Traits::GridViewType GridViewType;

	  //! \brief vector backend
	  typedef B BackendType;

	  //! \brief short cut for size type exported by Backend
	  typedef typename B::size_type SizeType;
	};


    class EmptyTransformation : public ConstraintsTransformation<int,float>
    {
    };


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
       \see SparsityLink
       
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
        typedef StdVectorFlatMatrixBackend Backend;

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
	  static const typename C::value_type& access (const C& c, size_type i, size_type j)
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
      template<typename LFSU, typename LFSV>
      static void pattern_volume (const LA& la, const LFSU& lfsu, const LFSV& lfsv, LocalSparsityPattern& pattern)
      {
      }
      template<typename LFSU, typename LFSV>
      static void pattern_volume_post_skeleton
      ( const LA& la,
        const LFSU& lfsu, const LFSV& lfsv,
        LocalSparsityPattern& pattern)
      {
      }
      template<typename LFSU, typename LFSV>
      static void pattern_skeleton (const LA& la, const LFSU& lfsu_s, const LFSV& lfsv_s, 
                                  const LFSU& lfsu_n, const LFSV& lfsv_n, 
                                   LocalSparsityPattern& pattern_sn,
                                   LocalSparsityPattern& pattern_ns)
      {
      }
      template<typename LFSU, typename LFSV>
      static void pattern_boundary(const LA& la,
                                   const LFSU& lfsu_s, const LFSV& lfsv_s,
                                   LocalSparsityPattern& pattern_ss)
      {
      }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
	  {
	  }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
	  {
	  }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_skeleton (const LA& la, const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           R& r_s, R& r_n)
      {
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_boundary (const LA& la, const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s)
      {
      }

	  template<typename EG, typename LFSV, typename R>
      static void lambda_volume (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
      }
	  template<typename EG, typename LFSV, typename R>
      static void lambda_volume_post_skeleton (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_skeleton(const LA& la, const IG& ig,
                                  const LFSV& lfsv_s, const LFSV& lfsv_n,
                                  R& r_s, R& r_n)
      {
      }
 	  template<typename IG, typename LFSV, typename R>
      static void lambda_boundary (const LA& la, const IG& ig, const LFSV& lfsv, R& r)
      {
      }

	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
	  {
	  }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
	  {
	  }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_skeleton (const LA& la, const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           Y& y_s, Y& y_n)
      {
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_boundary (const LA& la, const IG& ig, 
                                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                           Y& y_s)
      {
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, LocalMatrix<R>& mat)
      {
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, LocalMatrix<R>& mat)
      {
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_skeleton (const LA& la, const IG& ig, 
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                              LocalMatrix<R>& mat_ss, LocalMatrix<R>& mat_sn, 
                              LocalMatrix<R>& mat_ns, LocalMatrix<R>& mat_nn)
      {
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_boundary (const LA& la, const IG& ig, 
                                     const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                     LocalMatrix<R>& mat_ss)
      {
      }
    };
    template<typename LA>
    struct LocalAssemblerCallSwitch<LA,true>
    {
      template<typename LFSU, typename LFSV>
      static void pattern_volume (const LA& la, const LFSU& lfsu, const LFSV& lfsv, LocalSparsityPattern& pattern)
      {
        la.pattern_volume(lfsu,lfsv,pattern);
      }
      template<typename LFSU, typename LFSV>
      static void pattern_volume_post_skeleton
      ( const LA& la,
        const LFSU& lfsu, const LFSV& lfsv,
        LocalSparsityPattern& pattern)
      {
        la.pattern_volume_post_skeleton(lfsu,lfsv,pattern);
      }
      template<typename LFSU, typename LFSV>
      static void pattern_skeleton (const LA& la, const LFSU& lfsu_s, const LFSV& lfsv_s, 
                                  const LFSU& lfsu_n, const LFSV& lfsv_n, 
                                   LocalSparsityPattern& pattern_sn,
                                   LocalSparsityPattern& pattern_ns)
      {
        la.pattern_skeleton(lfsu_s,lfsv_s,lfsu_n,lfsv_n,
                            pattern_sn, pattern_ns);
      }
      template<typename LFSU, typename LFSV>
      static void pattern_boundary(const LA& la,
                                   const LFSU& lfsu_s, const LFSV& lfsv_s,
                                   LocalSparsityPattern& pattern_ss)
      {
        la.pattern_boundary(lfsu_s,lfsv_s,pattern_ss);
      }

	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
	  {
		la.alpha_volume(eg,lfsu,x,lfsv,r);
	  }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
	  {
		la.alpha_volume_post_skeleton(eg,lfsu,x,lfsv,r);
	  }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_skeleton (const LA& la, const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           R& r_s, R& r_n)
      {
        la.alpha_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,r_s,r_n);
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_boundary (const LA& la, const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s)
      {
        la.alpha_boundary(ig,lfsu_s,x_s,lfsv_s,r_s);
      }

	  template<typename EG, typename LFSV, typename R>
      static void lambda_volume (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
        la.lambda_volume(eg,lfsv,r);
      }
	  template<typename EG, typename LFSV, typename R>
      static void lambda_volume_post_skeleton (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
        la.lambda_volume_post_skeleton(eg,lfsv,r);
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_skeleton(const LA& la, const IG& ig,
                                  const LFSV& lfsv_s, const LFSV& lfsv_n,
                                  R& r_s, R& r_n)
      {
        la.lambda_skeleton(ig, lfsv_s, lfsv_n, r_s, r_n);
      }
 	  template<typename IG, typename LFSV, typename R>
      static void lambda_boundary (const LA& la, const IG& ig, const LFSV& lfsv, R& r)
      {
        la.lambda_boundary(ig,lfsv,r);
      }

	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
	  {
		la.jacobian_apply_volume(eg,lfsu,x,lfsv,y);
	  }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
	  {
		la.jacobian_apply_volume_post_skeleton(eg,lfsu,x,lfsv,y);
	  }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_skeleton (const LA& la, const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           Y& y_s, Y& y_n)
      {
        la.jacobian_apply_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,y_s,y_n);
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_boundary (const LA& la, const IG& ig, 
                                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                           Y& y_s)
      {
        la.jacobian_apply_boundary(ig,lfsu_s,x_s,lfsv_s,y_s);
      }

      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, LocalMatrix<R>& mat)
      {
        la.jacobian_volume(eg,lfsu,x,lfsv,mat);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, LocalMatrix<R>& mat)
      {
        la.jacobian_volume_post_skeleton(eg,lfsu,x,lfsv,mat);
      }
 	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_skeleton (const LA& la, const IG& ig, 
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                              LocalMatrix<R>& mat_ss, LocalMatrix<R>& mat_sn, 
                              LocalMatrix<R>& mat_ns, LocalMatrix<R>& mat_nn)
      {
        la.jacobian_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,
                             mat_ss, mat_sn, mat_ns, mat_nn);
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_boundary (const LA& la, const IG& ig, 
                                     const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                     LocalMatrix<R>& mat_ss)
      {
        la.jacobian_boundary(ig,lfsu_s,x_s,lfsv_s,mat_ss);
      }
   };

  } // namespace PDELab
} // namespace Dune

#endif
