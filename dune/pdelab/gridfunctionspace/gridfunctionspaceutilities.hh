// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACEUTILITIES_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACEUTILITIES_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include <dune/common/static_assert.hh>

#include"../common/countingptr.hh"
#include"../common/multitypetree.hh"
#include"../common/cpstoragepolicy.hh"
#include"../common/function.hh"

#include"gridfunctionspace.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //===============================================================
    // output: convert grid function space to discrete grid function
    //===============================================================


	/** \brief convert a single component function space into a grid function
     *
     * The functions can be vector-valued.
     *
     * This is just an intermediate solution to provide VTK output.
     *
     * \tparam T Type of GridFunctionSpace
     * \tparam X Type of coefficients vector
     */
	template<typename T, typename X>
	class DiscreteGridFunction
	  : public GridFunctionInterface<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
            T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
            typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType
            >,
          DiscreteGridFunction<T,X>
          >
	{
	  typedef T GFS;

	  typedef GridFunctionInterface<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
          T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
          typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType
          >,
        DiscreteGridFunction<T,X>
        > BaseT;
	
	public:
	  typedef typename BaseT::Traits Traits;
	  
      /** \brief Construct a DiscreteGridFunction
       *
       * \param gfs The GridFunctionsSpace
       * \param x_  The coefficients vector
       */
	  DiscreteGridFunction (const GFS& gfs, const X& x_)
		: pgfs(&gfs), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), yb(gfs.maxLocalSize())
	  {
	  }

      // Evaluate
	  inline void evaluate (const typename Traits::ElementType& e, 
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {  
		lfs.bind(e);
		lfs.vread(xg,xl);
		lfs.localFiniteElement().localBasis().evaluateFunction(x,yb);
		y = 0;
		for (unsigned int i=0; i<yb.size(); i++)
		  y.axpy(xl[i],yb[i]);
 	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView ()
	  {
		return pgfs->gridview();
	  }

	private:
	  CP<GFS const> pgfs;
	  const X& xg;
	  mutable typename GFS::LocalFunctionSpace lfs;
	  mutable std::vector<typename Traits::RangeFieldType> xl;
	  mutable std::vector<typename Traits::RangeType> yb;
	};

	//! \brief convert a single component function space with experimental
	//! global finite elements into a grid function
    /**
     * The functions can be vector-valued.
     *
     * This is just an intermediate solution to provide VTK output.
     *
     * \tparam T Type of GridFunctionSpace.  The LocalBasis must provide the
     *           evaluateFunctionGlobal() method.
     * \tparam X Type of coefficients vector
     */
	template<typename T, typename X>
	class DiscreteGridFunctionGlobal
	  : public GridFunctionInterface<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
            T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
            typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType
            >,
          DiscreteGridFunctionGlobal<T,X>
          >
	{
	  typedef T GFS;

	  typedef GridFunctionInterface<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
          T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
          typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType
          >,
        DiscreteGridFunctionGlobal<T,X>
        > BaseT;
	
	public:
	  typedef typename BaseT::Traits Traits;
	  
      /** \brief Construct a DiscreteGridFunctionGlobal
       *
       * \param gfs The GridFunctionsSpace
       * \param x_  The coefficients vector
       */
	  DiscreteGridFunctionGlobal (const GFS& gfs, const X& x_)
		: pgfs(&gfs), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), yb(gfs.maxLocalSize())
	  {
	  }

      // Evaluate
	  inline void evaluate (const typename Traits::ElementType& e, 
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {  
		lfs.bind(e);
		lfs.vread(xg,xl);
		lfs.localFiniteElement().localBasis().evaluateFunctionGlobal(x,yb,e.geometry());
		y = 0;
		for (unsigned int i=0; i<yb.size(); i++)
		  y.axpy(xl[i],yb[i]);
 	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView ()
	  {
		return pgfs->gridview();
	  }

	private:
	  CP<GFS const> pgfs;
	  const X& xg;
	  mutable typename GFS::LocalFunctionSpace lfs;
	  mutable std::vector<typename Traits::RangeFieldType> xl;
	  mutable std::vector<typename Traits::RangeType> yb;
	};

	//! \brief convert a single component function space with experimental
	//! global finite elements into a grid function representing the curl
    /**
     * The function values should be 3-component vectors.  The Curl will be a
     * 3-component function.
     *
     * This is just an intermediate solution to provide VTK output.
     *
     * \tparam T Type of GridFunctionSpace.  The LocalBasis must provide the
     *           evaluateJacobianGlobal() method.
     * \tparam X Type of coefficients vector
     */
	template<typename T, typename X>
	class DiscreteGridFunctionGlobalCurl
	  : public GridFunctionInterface<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
            T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
            typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType
            >,
          DiscreteGridFunctionGlobalCurl<T,X>
          >
	{
      dune_static_assert(
        T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange == 3,
        "Range dimension of localbasis must be 3 for DiscreteGridFunctionGlobalCurl");
      dune_static_assert(
        T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimDomain == 3,
        "Domain dimension of localbasis must be 3 for DiscreteGridFunctionGlobalCurl");
	  typedef T GFS;
	  typedef GridFunctionInterface<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
          T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
          typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType
          >,
        DiscreteGridFunctionGlobalCurl<T,X>
        > BaseT;

	public:
	  typedef typename BaseT::Traits Traits;
	  
      /** \brief Construct a DiscreteGridFunctionGlobalCurl
       *
       * \param gfs The GridFunctionsSpace
       * \param x_  The coefficients vector
       */
	  DiscreteGridFunctionGlobalCurl (const GFS& gfs, const X& x_)
		: pgfs(&gfs), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), J(gfs.maxLocalSize())
	  {
	  }

      // Evaluate
	  inline void evaluate (const typename Traits::ElementType& e, 
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
		lfs.bind(e);
		lfs.vread(xg,xl);
		lfs.localFiniteElement().localBasis().evaluateJacobianGlobal(x,J,e.geometry());
        y = 0;
		for (unsigned int i=0; i<J.size(); i++) {
		  y[0] += xl[i]*(J[i][2][1] - J[i][1][2]);
		  y[1] += xl[i]*(J[i][0][2] - J[i][2][0]);
		  y[2] += xl[i]*(J[i][1][0] - J[i][0][1]);
        }
 	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView ()
	  {
		return pgfs->gridview();
	  }

	private:
	  CP<GFS const> pgfs;
	  const X& xg;
	  mutable typename GFS::LocalFunctionSpace lfs;
	  mutable std::vector<typename Traits::RangeFieldType> xl;
      mutable std::vector<typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::JacobianType> J;
	};

	//! global finite elements into a grid function representing the curl
    /**
     * The function values should be 2-component vectors.  The Curl will be a
     * 1-component function.  (If the function itself has values in the
     * x-y-plane, the curl will point in z-direction).  It is assumed the the
     * container used to aggregate the components is a FieldVector.
     *
     * This is just an intermediate solution to provide VTK output.
     *
     * \tparam T Type of GridFunctionSpace.  The LocalBasis must provide the
     *           evaluateJacobianGlobal() method.
     * \tparam X Type of coefficients vector
     */
	template<typename T, typename X>
	class DiscreteGridFunctionGlobalCurl2D
	  : public GridFunctionInterface<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
            1,
            FieldVector<
              typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
              1
              >
            >,
          DiscreteGridFunctionGlobalCurl2D<T,X>
          >
	{
      dune_static_assert(
        T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange == 2,
        "Range dimension of localbasis must be 2 for DiscreteGridFunctionGlobalCurl2D");
      dune_static_assert(
        T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimDomain == 2,
        "Domain dimension of localbasis must be 2 for DiscreteGridFunctionGlobalCurl2D");
	  typedef T GFS;
	  typedef GridFunctionInterface<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
          1,
          FieldVector<
            typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
            1
            >
          >,
        DiscreteGridFunctionGlobalCurl2D<T,X>
        > BaseT;

	public:
	  typedef typename BaseT::Traits Traits;
	  
      /** \brief Construct a DiscreteGridFunctionGlobalCurl2D
       *
       * \param gfs The GridFunctionsSpace
       * \param x_  The coefficients vector
       */
	  DiscreteGridFunctionGlobalCurl2D (const GFS& gfs, const X& x_)
		: pgfs(&gfs), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), J(gfs.maxLocalSize())
	  {
	  }

      // Evaluate
	  inline void evaluate (const typename Traits::ElementType& e, 
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
		lfs.bind(e);
		lfs.vread(xg,xl);
		lfs.localFiniteElement().localBasis().evaluateJacobianGlobal(x,J,e.geometry());
        y = 0;
		for (unsigned int i=0; i<J.size(); i++)
		  y += xl[i]*(J[i][1][0] - J[i][0][1]);
 	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView ()
	  {
		return pgfs->gridview();
	  }

	private:
	  CP<GFS const> pgfs;
	  const X& xg;
	  mutable typename GFS::LocalFunctionSpace lfs;
	  mutable std::vector<typename Traits::RangeFieldType> xl;
      mutable std::vector<typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::JacobianType> J;
	};

    /** \brief DiscreteGridFunction with Piola transformation 
     *
     * \copydetails DiscreteGridFunction
     */
	template<typename T, typename X>
	class DiscreteGridFunctionPiola
	  : public GridFunctionInterface<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
            T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
            typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType
            >,
          DiscreteGridFunctionPiola<T,X>
          >
	{
	  typedef T GFS;

	  typedef GridFunctionInterface<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
          T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
          typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType
          >,
        DiscreteGridFunctionPiola<T,X>
        > BaseT;
	
	public:
	  typedef typename BaseT::Traits Traits;
	  
      /** \brief Construct a DiscreteGridFunctionPiola
       *
       * \copydetails DiscreteGridFunction::DiscreteGridFunction(const GFS&,const X&)
       */
	  DiscreteGridFunctionPiola (const GFS& gfs, const X& x_)
		: pgfs(&gfs), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), yb(gfs.maxLocalSize())
	  {
	  }

	  inline void evaluate (const typename Traits::ElementType& e, 
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  { 
        // evaluate shape function on the reference element as before
		lfs.bind(e);
		lfs.vread(xg,xl);
		lfs.localFiniteElement().localBasis().evaluateFunction(x,yb);
        typename Traits::RangeType yhat;
        yhat = 0;
		for (unsigned int i=0; i<yb.size(); i++)
		  yhat.axpy(xl[i],yb[i]);

        // apply Piola transformation
        Dune::FieldMatrix<typename Traits::DomainFieldType,
          GFS::Traits::GridViewType::dimension,GFS::Traits::GridViewType::dimension>
          J = e.geometry().jacobianInverseTransposed(x);
        J.invert();
        y = 0;
        J.umtv(yhat,y);
        y /= J.determinant();
	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView ()
	  {
		return pgfs->gridview();
	  }

	private:
	  CP<GFS const> pgfs;
	  const X& xg;
	  mutable typename GFS::LocalFunctionSpace lfs;
	  mutable std::vector<typename Traits::RangeFieldType> xl;
	  mutable std::vector<typename Traits::RangeType> yb;
	};


    /** \brief DiscreteGridFunction with transformation for edge elements
     *
     * \copydetails DiscreteGridFunction
     */
	template<typename T, typename X>
	class DiscreteGridFunctionEdge
	  : public GridFunctionInterface<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
            T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
            typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType
            >,
          DiscreteGridFunctionEdge<T,X>
          >
	{
	  typedef T GFS;

	  typedef GridFunctionInterface<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
          T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange,
          typename T::Traits::LocalFiniteElementType::Traits::LocalBasisType::Traits::RangeType
          >,
        DiscreteGridFunctionEdge<T,X>
        > BaseT;
	
	public:
	  typedef typename BaseT::Traits Traits;
	  
      /** \brief Construct a DiscreteGridFunctionEdge
       *
       * \copydetails DiscreteGridFunction::DiscreteGridFunction(const GFS&,const X&)
       */
	  DiscreteGridFunctionEdge (const GFS& gfs, const X& x_)
		: pgfs(&gfs), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), yb(gfs.maxLocalSize())
	  {
	  }

	  inline void evaluate (const typename Traits::ElementType& e, 
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  { 
        // evaluate shape function on the reference element as before
		lfs.bind(e);
		lfs.vread(xg,xl);
		lfs.localFiniteElement().localBasis().evaluateFunction(x,yb);
        typename Traits::RangeType yhat;
        yhat = 0;
		for (unsigned int i=0; i<yb.size(); i++)
		  yhat.axpy(xl[i],yb[i]);

        // apply edge transformation
        Dune::FieldMatrix<typename Traits::DomainFieldType,
          GFS::Traits::GridViewType::dimension,GFS::Traits::GridViewType::dimension>
          JTInv = e.geometry().jacobianInverseTransposed(x);
        y = 0;
        JTInv.umv(yhat,y);
        // det J = 1/det JInv
        y /= JTInv.determinant();
	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView ()
	  {
		return pgfs->gridview();
	  }

	private:
	  CP<GFS const> pgfs;
	  const X& xg;
	  mutable typename GFS::LocalFunctionSpace lfs;
	  mutable std::vector<typename Traits::RangeFieldType> xl;
	  mutable std::vector<typename Traits::RangeType> yb;
	};


	// convert a power function space of scalar function spaces into a vector-valued grid function
    // this is just an intermediate solution to provide VTK output
	template<typename T, typename X>
	class VectorDiscreteGridFunction
	  : public GridFunctionInterface<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            typename T::template Child<0>::Type::Traits::LocalFiniteElementType
                     ::Traits::LocalBasisType::Traits::RangeFieldType,
            T::CHILDREN,
            Dune::FieldVector<
              typename T::template Child<0>::Type::Traits::LocalFiniteElementType
                       ::Traits::LocalBasisType::Traits::RangeFieldType,
              T::CHILDREN
              >
            >,
          VectorDiscreteGridFunction<T,X>
          >
	{
	  typedef T GFS;

      typedef GridFunctionInterface<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::template Child<0>::Type::Traits::LocalFiniteElementType
                   ::Traits::LocalBasisType::Traits::RangeFieldType,
          T::CHILDREN,
          Dune::FieldVector<
            typename T::template Child<0>::Type::Traits::LocalFiniteElementType
                     ::Traits::LocalBasisType::Traits::RangeFieldType,
            T::CHILDREN
            >
          >,
        VectorDiscreteGridFunction<T,X>
        > BaseT;

	public:
	  typedef typename BaseT::Traits Traits;
      typedef typename T::template Child<0>::Type ChildType;
      typedef typename ChildType::Traits::LocalFiniteElementType
                       ::Traits::LocalBasisType::Traits::RangeFieldType RF;
      typedef typename ChildType::Traits::LocalFiniteElementType
                       ::Traits::LocalBasisType::Traits::RangeType RT;
	  
	  VectorDiscreteGridFunction (const GFS& gfs, const X& x_)
		: pgfs(&gfs), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), yb(gfs.maxLocalSize())
	  {
	  }

	  inline void evaluate (const typename Traits::ElementType& e, 
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {  
		lfs.bind(e);
		lfs.vread(xg,xl);
        for (int k=0; k<T::CHILDREN; k++)
          {
            lfs.getChild(k).localFiniteElement().localBasis().evaluateFunction(x,yb);
            y[k] = 0.0;
            for (unsigned int i=0; i<yb.size(); i++)
              y[k] += xl[lfs.getChild(k).localIndex(i)]*yb[i];
          }
	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView ()
	  {
		return pgfs->gridview();
	  }

	private:
	  CP<GFS const> pgfs;
	  const X& xg;
	  mutable typename GFS::LocalFunctionSpace lfs;
	  mutable std::vector<RF> xl;
	  mutable std::vector<RT> yb;
	};

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
