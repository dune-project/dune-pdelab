// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATORSPACEUTILITIES_HH
#define DUNE_PDELAB_GRIDOPERATORSPACEUTILITIES_HH

#include <vector>

#include <dune/common/static_assert.hh>
#include <dune/common/tuples.hh>
#include <dune/common/typetraits.hh>

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

    //////////////////////////////////////////////////////////////////////
    //
    //  Stuff for dealing with pattern_skeleton
    //

    //! Determine whether a local operator has the old style pattern_skeleton()
    /**
     * This template determines whether \c LA has an old style (deprecated)
     * pattern_skeleton() method, i.e. one which may be called as
     * \code
la.pattern_skeleton(lfsu, lfsv, lfsu, lfsv, pattern, pattern)
     * \endcode
     * where \c la is of type \c const \c LA, \c lfsu is of type \c const \c
     * LFSU, \c lfsv is of type \c const \c LFSV, and \c pattern is of type
     * LocalSparsityPattern.
     *
     * If such a method exists, this class derives from true_type, otherwise
     * it derives from false_type.
     *
     * \tparam LA   Type of LocalOperator.
     * \tparam LFSU Type of trial LocalFunctionSpace.
     * \tparam LFSV Type of test LocalFunctionSpace.
     *
     * \note The last (unnamed) template parameter is an implementation detail
     *       -- when instanciating this class, you should omit this parameter
     *       from the template parameter list.
     */
    template<typename LA, typename LFSU, typename LFSV, typename = void>
    struct LocalOperatorHasOldPatternSkeleton
      : public false_type
    { };
    template<typename LA, typename LFSU, typename LFSV>
    struct LocalOperatorHasOldPatternSkeleton
    < LA, LFSU, LFSV,
      // The parenthesis after the following sizeof do not enclose an argument
      // list -- they are for scoping only and enclose a comma-expression.
      // This is necessary because pattern_skeleton() will most probably yield
      // an expression of type void, and it is illegal to apply the sizeof
      // operator to void.  However, the comma operator can have void
      // expressions as its arguments.  So we take the function call as the
      // left argument of the comma operator, and some arbitrary value whith a
      // size > 0 as the right argument of the comma operator.  The result is
      // the right argument, and we can apply the sizeof operator.
      //
      // That the theory.  In reality, g++ (4.4 at least, svn for 4.6 is
      // fixed) has a bug which lets it confuse specializations of this
      // template even for different class names as long as the name of the
      // method (pattern_skeleton in this case) matches.  The workaround is to
      // make the second argument of the comma-operator differ for different
      // traits classes.  To minimize the likelyhood of some other programmer
      // writing a traits class testing for a method pattern_skeleton, we use
      // a random value here.
      typename enable_if<sizeof( static_cast<LA*const>(0)->pattern_skeleton
                                 ( *static_cast<LFSU*const>(0),
                                   *static_cast<LFSV*const>(0),
                                   *static_cast<LFSU*const>(0),
                                   *static_cast<LFSV*const>(0),
                                   *static_cast<LocalSparsityPattern*>(0),
                                   *static_cast<LocalSparsityPattern*>(0)),
                                 0xd124f21d )>::type >
      : public true_type
    { };

    //! Determine whether a local operator has the new style pattern_skeleton()
    /**
     * This template determines whether \c LA has an new style
     * pattern_skeleton() method, i.e. one which may be called as
     * \code
la.pattern_skeleton(lfsu, lfsv, lfsu, lfsv, pattern, pattern, pattern, pattern)
     * \endcode
     * where \c la is of type \c const \c LA, \c lfsu is of type \c const \c
     * LFSU, \c lfsv is of type \c const \c LFSV, and \c pattern is of type
     * LocalSparsityPattern.
     *
     * If such a method exists, this class derives from true_type, otherwise
     * it derives from false_type.
     *
     * \tparam LA   Type of LocalOperator.
     * \tparam LFSU Type of trial LocalFunctionSpace.
     * \tparam LFSV Type of test LocalFunctionSpace.
     *
     * \note The last (unnamed) template parameter is an implementation detail
     *       -- when instanciating this class, you should omit this parameter
     *       from the template parameter list.
     */
    template<typename LA, typename LFSU, typename LFSV, typename = void>
    struct LocalOperatorHasNewPatternSkeleton
      : public false_type
    { };
    template<typename LA, typename LFSU, typename LFSV>
    struct LocalOperatorHasNewPatternSkeleton
    < LA, LFSU, LFSV,
      // Read the comment in LocalOperatorHasNewPatternSkeleton.  It explains
      // the construct here and how to work around a related g++ compiler bug.
      typename enable_if<sizeof( static_cast<LA*const>(0)->pattern_skeleton
                                 ( *static_cast<LFSU*const>(0),
                                   *static_cast<LFSV*const>(0),
                                   *static_cast<LFSU*const>(0),
                                   *static_cast<LFSV*const>(0),
                                   *static_cast<LocalSparsityPattern*>(0),
                                   *static_cast<LocalSparsityPattern*>(0),
                                   *static_cast<LocalSparsityPattern*>(0),
                                   *static_cast<LocalSparsityPattern*>(0)),
                                 0xff4bd551 )>::type >
      : public true_type
    { };

    //! \brief Call the old style or the new style pattern_skeleton(),
    //!        depending on which is actually present
    /**
     * \tparam LA      Type of LocalOperator.
     * \tparam LFSU    Type of trial LocalFunctionSpace.
     * \tparam LFSV    Type of test LocalFunctionSpace.
     * \tparam has_old Whether the Local operator has the old style
     *                 pattern_skeleton() method.
     * \tparam has_new Whether the Local operator has the new style
     *                 pattern_skeleton() method.
     *
     * \note \c has_old and \c has_new are implementation details.  You should
     *       omit them when instanciating this class -- the default values
     *       already to the right thing.
     */
    template<typename LA, typename LFSU, typename LFSV,
             bool has_old = LocalOperatorHasOldPatternSkeleton<LA, LFSU, LFSV>::value,
             bool has_new = LocalOperatorHasNewPatternSkeleton<LA, LFSU, LFSV>::value>
    struct PatternSkeletonDeprecationCallSwitch {
      //! Call the local operators pattern_skeleton()
      /**
       * What happens exactly depends on whether the old style or the new
       * style methods are present on the local operator:
       * \li both old style and new style missing: the call is a noop,
       * \li old style present but new style missing: call the old style
       *     pattern_skeleton() method, \c pattern_ss and \c pattern_nn sty
       *     untouched,
       * \li old style missing but new style present: call the new style
       *     pattern_skeleton() method, and finally
       * \li both old style and new style present: throw a static assertion
       *     since this is most probably and error.
       */
      static void pattern_skeleton(const LA& la,
                                   const LFSU& lfsu_s, const LFSV& lfsv_s,
                                   const LFSU& lfsu_n, const LFSV& lfsv_n,
                                   LocalSparsityPattern& pattern_ss,
                                   LocalSparsityPattern& pattern_sn,
                                   LocalSparsityPattern& pattern_ns,
                                   LocalSparsityPattern& pattern_nn)
      { }
    };

    template<typename LA, typename LFSU, typename LFSV>
    struct PatternSkeletonDeprecationCallSwitch<LA, LFSU, LFSV, true, false> {
      static void pattern_skeleton(const LA& la,
                                   const LFSU& lfsu_s, const LFSV& lfsv_s,
                                   const LFSU& lfsu_n, const LFSV& lfsv_n,
                                   LocalSparsityPattern& pattern_ss,
                                   LocalSparsityPattern& pattern_sn,
                                   LocalSparsityPattern& pattern_ns,
                                   LocalSparsityPattern& pattern_nn)
        DUNE_DEPRECATED
      {
        la.pattern_skeleton(lfsu_s, lfsv_s, lfsu_n, lfsv_n,
                            pattern_sn, pattern_ns);
      }
    };

    template<typename LA, typename LFSU, typename LFSV>
    struct PatternSkeletonDeprecationCallSwitch<LA, LFSU, LFSV, false, true> {
      static void pattern_skeleton(const LA& la,
                                   const LFSU& lfsu_s, const LFSV& lfsv_s,
                                   const LFSU& lfsu_n, const LFSV& lfsv_n,
                                   LocalSparsityPattern& pattern_ss,
                                   LocalSparsityPattern& pattern_sn,
                                   LocalSparsityPattern& pattern_ns,
                                   LocalSparsityPattern& pattern_nn)
      {
        la.pattern_skeleton(lfsu_s, lfsv_s, lfsu_n, lfsv_n,
                            pattern_ss, pattern_sn, pattern_ns, pattern_nn);
      }
    };

    template<typename LA, typename LFSU, typename LFSV>
    class PatternSkeletonDeprecationCallSwitch<LA, LFSU, LFSV, true, true> {
      dune_static_assert(AlwaysFalse<LA>::value, "Both the oldstyle "
                         "(deprecated) and the newstyle patter_skeleton() "
                         "methods exist for the local operator");
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
                                   LocalSparsityPattern& pattern_ss,
                                   LocalSparsityPattern& pattern_sn,
                                   LocalSparsityPattern& pattern_ns,
                                   LocalSparsityPattern& pattern_nn)
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
                                   LocalSparsityPattern& pattern_ss,
                                   LocalSparsityPattern& pattern_sn,
                                   LocalSparsityPattern& pattern_ns,
                                   LocalSparsityPattern& pattern_nn)
      {
        PatternSkeletonDeprecationCallSwitch<LA, LFSU, LFSV>::
          pattern_skeleton(la, lfsu_s,lfsv_s,lfsu_n,lfsv_n,
                           pattern_ss, pattern_sn, pattern_ns, pattern_nn);
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
