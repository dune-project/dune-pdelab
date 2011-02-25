// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACEUTILITIES_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACEUTILITIES_HH

#include <cstdlib>
#include<vector>

#include<dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include"../common/function.hh"
#include <dune/pdelab/common/jacobiantocurl.hh>
#include <dune/pdelab/gridfunctionspace/constraintstransformation.hh> // backward compatibility
#include"gridfunctionspace.hh"
#include <dune/pdelab/gridfunctionspace/ordering.hh> // backward compatibility


namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //===============================================================
    // output: convert grid function space to discrete grid function
    //===============================================================


    /** \brief convert a grid function space and a coefficient vector into a
     *         grid function
     *
     * If a GridFunctionSpace with local-valued finite elements is used, this
     * class can only be used for scalar functions, since for vector-valued
     * local finite elements the values must be transformed, and the
     * transformation depends on the type of element.  For H(div) elements
     * (Raviart-Thomas) look at DiscreteGridFunctionPiola.
     *
     * If a GridFunctionSpace with finite elements using the new global-valued
     * interface is used, this class can be used as-is even for vector-valued
     * functions.
     *
     * If you have a GridFunctionSpace tree of 1-component grid-function
     * spaces, and want to collectively treat them as a vector-valued
     * grid-function, look at VectorDiscreteGridFunction.
     *
     * \tparam T Type of GridFunctionSpace
     * \tparam X Type of coefficients vector
     */
	template<typename T, typename X>
	class DiscreteGridFunction
	  : public TypeTree::LeafNode
      , GridFunctionInterface<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            typename BasisInterfaceSwitch<
              typename FiniteElementInterfaceSwitch<
                typename T::Traits::FiniteElementType
                >::Basis
              >::RangeField,
            BasisInterfaceSwitch<
              typename FiniteElementInterfaceSwitch<
                typename T::Traits::FiniteElementType
                >::Basis
              >::dimRange,
            typename BasisInterfaceSwitch<
              typename FiniteElementInterfaceSwitch<
                typename T::Traits::FiniteElementType
                >::Basis
              >::Range
            >,
          DiscreteGridFunction<T,X>
        >
	{
	  typedef T GFS;

      typedef typename Dune::BasisInterfaceSwitch<
        typename FiniteElementInterfaceSwitch<
          typename T::Traits::FiniteElementType
          >::Basis
        > BasisSwitch;
	  typedef GridFunctionInterface<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename BasisSwitch::RangeField,
          BasisSwitch::dimRange,
          typename BasisSwitch::Range
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
		: pgfs(stackobject_to_shared_ptr(gfs)), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), yb(gfs.maxLocalSize())
	  {
	  }

      // Evaluate
	  inline void evaluate (const typename Traits::ElementType& e,
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
        typedef FiniteElementInterfaceSwitch<
          typename Dune::PDELab::LocalFunctionSpace<GFS>::Traits::FiniteElementType
          > FESwitch;
		lfs.bind(e);
		lfs.vread(xg,xl);
        FESwitch::basis(lfs.finiteElement()).evaluateFunction(x,yb);
		y = 0;
		for (unsigned int i=0; i<yb.size(); i++)
		  y.axpy(xl[i],yb[i]);
 	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView () const
	  {
		return pgfs->gridview();
	  }

	private:
	  shared_ptr<GFS const> pgfs;
	  const X& xg;
	  mutable typename Dune::PDELab::LocalFunctionSpace<GFS> lfs;
	  mutable std::vector<typename Traits::RangeFieldType> xl;
	  mutable std::vector<typename Traits::RangeType> yb;
	};

    /** \brief convert a grid function space and a coefficient vector into a
     *         grid function of the curl
     *
     * This class works only with a GridFunctionSpace with finite elements
     * using the new global-valued interface.
     *
     * \tparam T Type of GridFunctionSpace
     * \tparam X Type of coefficients vector
     */
    template<typename T, typename X>
    class DiscreteGridFunctionCurl :
      public GridFunctionInterface<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename JacobianToCurl<typename T::Traits::FiniteElementType::
                                  Traits::Basis::Traits::Jacobian>::CurlField,
          JacobianToCurl<typename T::Traits::FiniteElementType::Traits::Basis::
                         Traits::Jacobian>::dimCurl,
          typename JacobianToCurl<typename T::Traits::FiniteElementType::
                                  Traits::Basis::Traits::Jacobian>::Curl
          >,
        DiscreteGridFunctionCurl<T,X>
        >
    {
      typedef T GFS;
      typedef typename T::Traits::FiniteElementType::Traits::Basis::Traits::
        Jacobian Jacobian;
      typedef JacobianToCurl<Jacobian> J2C;

    public:
      typedef GridFunctionTraits<
        typename T::Traits::GridViewType,
          typename J2C::CurlField, J2C::dimCurl, typename J2C::Curl
        > Traits;

    private:
      typedef GridFunctionInterface<Traits, DiscreteGridFunctionCurl<T,X> >
        BaseT;

      const GFS &gfs;
      const X &xg;

    public:
      /** \brief Construct a DiscreteGridFunctionCurl
       *
       * \param gfs_ The GridFunctionsSpace
       * \param x_   The coefficients vector
       */
      DiscreteGridFunctionCurl(const GFS& gfs_, const X& xg_) :
        gfs(gfs_), xg(xg_)
      { }

      // Evaluate
      void evaluate (const typename Traits::ElementType& e,
                     const typename Traits::DomainType& x,
                     typename Traits::RangeType& y) const
      {
        static const J2C& j2C = J2C();

        typename GFS::LocalFunctionSpace lfs(gfs);
        lfs.bind(e);
        std::vector<typename Traits::RangeFieldType> xl(lfs.size());
        lfs.vread(xg,xl);
        std::vector<Jacobian> jacobian(lfs.size());
        lfs.finiteElement().basis().evaluateJacobian(x,jacobian);

        y = 0;
        typename Traits::RangeType yb;
        for (std::size_t i=0; i < lfs.size(); i++) {
          j2C(jacobian[i], yb);
          y.axpy(xl[i], yb);
        }
      }

      //! get a reference to the GridView
      const typename Traits::GridViewType& getGridView() const
      { return gfs.gridview(); }
    };

    //! Helper class to calculate the Traits of DiscreteGridFunctionCurl
    /**
     * \tparam GV              Type of the GridView.
     * \tparam RangeFieldType  RangeFieldType of the basis, and resulting
     *                         RangeFieldType of this function.
     * \tparam dimRangeOfBasis Number of components of the function to take
     *                         the curl of, a.k.a. th dimRange of the basis.
     *
     * \note This the non-specialized version of the
     *       DiscreteGridFunctionCurlTraits template.  It must be specialized
     *       for different values of dimRangeOfBasis.  If this non-specialized
     *       version is instantiated, dune_static_assert() will be triggered.
     */
    template<typename GV, typename RangeFieldType, int dimRangeOfBasis>
    struct DiscreteGridFunctionCurlTraits {
      dune_static_assert(AlwaysFalse<GV>::value,
                         "DiscreteGridFunctionCurl (and friends) work in 2D "
                         "and 3D only");
    };
    //! Helper class to calculate the Traits of DiscreteGridFunctionCurl (1D)
    /**
     *  This is the specialization for dimRangeOfBasis == 1.  It takes the
     *  curl of a scalar valued function in a 2D space, i.e. a function with
     *  dimRange == 1 and dimDomain == 2.  The curl itself will have dimRange
     *  == 2.
     */
    template<typename GV, typename RangeFieldType>
    struct DiscreteGridFunctionCurlTraits<GV, RangeFieldType, 1>
      : public GridFunctionTraits<GV,
                                  RangeFieldType, 2,
                                  FieldVector<RangeFieldType, 2> >
    {
      dune_static_assert(GV::dimensionworld == 2,
                         "World dimension of grid must be 2 for the curl of a "
                         "scalar (1D) quantity");
    };
    //! Helper class to calculate the Traits of DiscreteGridFunctionCurl (2D)
    /**
     *  This is the specialization for dimRangeOfBasis == 2.  It takes the
     *  curl of a function with dimRange == 2 and dimDomain == 2.  The curl
     *  itself will have dimRange == 1.
     */
    template<typename GV, typename RangeFieldType>
    struct DiscreteGridFunctionCurlTraits<GV, RangeFieldType, 2>
      : public GridFunctionTraits<GV,
                                  RangeFieldType, 1,
                                  FieldVector<RangeFieldType, 1> >
    {
      dune_static_assert(GV::dimensionworld == 2,
                         "World dimension of grid must be 2 for the curl of a"
                         "2D quantity");
    };
    //! Helper class to calculate the Traits of DiscreteGridFunctionCurl (3D)
    /**
     *  This is the specialization for dimRangeOfBasis == 3.  It takes the
     *  curl of a function with dimRange == 3 and dimDomain == 3.  The curl
     *  itself will have dimRange == 3.
     */
    template<typename GV, typename RangeFieldType>
    struct DiscreteGridFunctionCurlTraits<GV, RangeFieldType, 3>
      : public GridFunctionTraits<GV,
                                  RangeFieldType, 3,
                                  FieldVector<RangeFieldType, 3> >
    {
      dune_static_assert(GV::dimensionworld == 3,
                         "World dimension of grid must be 3 for the curl of a"
                         "3D quantity");
    };

	//! \brief convert a single component function space with experimental
	//! global finite elements into a grid function representing the curl
    /**
     * For dimDomain=dimRange=3 the curl will be a 3-component function.  For
     * dimDomain=2 (x- and y-coordinates present) and dimRange=2 (x- and
     * y-components present) the curl will be a 1-component function
     * (z-component present).  For dimDomain=2 (x- and y-coordinates present)
     * and dimRange=1 (z-component present) the curl will be a 2-component
     * function (x- and y-components present).
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
          DiscreteGridFunctionCurlTraits<
            typename T::Traits::GridViewType,
            typename T::Traits::FiniteElementType::Traits::
               LocalBasisType::Traits::RangeFieldType,
            T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::
               dimRange>,
          DiscreteGridFunctionGlobalCurl<T,X> >
	{
    public:
      typedef DiscreteGridFunctionCurlTraits<
        typename T::Traits::GridViewType,
        typename T::Traits::FiniteElementType::Traits::
          LocalBasisType::Traits::RangeFieldType,
        T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::
          dimRange> Traits;

    private:
	  typedef T GFS;
	  typedef GridFunctionInterface<
        Traits,
        DiscreteGridFunctionGlobalCurl<T,X> > BaseT;
      typedef typename T::Traits::FiniteElementType::Traits::
        LocalBasisType::Traits LBTraits;

	public:
      /** \brief Construct a DiscreteGridFunctionGlobalCurl
       *
       * \param gfs The GridFunctionsSpace
       * \param x_  The coefficients vector
       */
	  DiscreteGridFunctionGlobalCurl (const GFS& gfs, const X& x_)
		: pgfs(stackobject_to_shared_ptr(gfs)), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), J(gfs.maxLocalSize())
	  {
	  }

      // Evaluate
	  inline void evaluate (const typename Traits::ElementType& e,
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
		lfs.bind(e);
		lfs.vread(xg,xl);
        lfs.finiteElement().localBasis().
          evaluateJacobianGlobal(x,J,e.geometry());
        y = 0;
        for (unsigned int i=0; i<J.size(); i++)
          // avoid a "case label value exceeds maximum value for type"
          // warning: since dimRange is an anonymous enum, its type may
          // contain only the values 0 and 1, resulting in a warning.
          switch(unsigned(Traits::dimRange)) {
          case 1:
            y[0] += xl[i] *  J[i][0][1];
            y[1] += xl[i] * -J[i][0][0];
            break;
          case 2:
            y[0] += xl[i]*(J[i][1][0] - J[i][0][1]);
            break;
          case 3:
            y[0] += xl[i]*(J[i][2][1] - J[i][1][2]);
            y[1] += xl[i]*(J[i][0][2] - J[i][2][0]);
            y[2] += xl[i]*(J[i][1][0] - J[i][0][1]);
            break;
          default:
            //how did that get passed all the static asserts?
            std::abort();
          }
 	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView () const
	  {
		return pgfs->gridview();
	  }

	private:
	  shared_ptr<GFS const> pgfs;
	  const X& xg;
	  mutable typename GFS::LocalFunctionSpace lfs;
	  mutable std::vector<typename Traits::RangeFieldType> xl;
      mutable std::vector<typename T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType> J;
	};

    //! \brief convert a single component function space with a grid function
    //! representing the gradient
    /**
     * The function values should be single-component vectors.  The Gradien
     * will be a dimDomain-component function.
     *
     * \tparam T Type of GridFunctionSpace.  The LocalBasis must provide the
     *           evaluateJacobian() method.
     * \tparam X Type of coefficients vector
     */
    template<typename T, typename X>
    class DiscreteGridFunctionGradient
      : public GridFunctionInterface<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            typename T::Traits::FiniteElementType::Traits::LocalBasisType
              ::Traits::RangeFieldType,
            T::Traits::FiniteElementType::Traits::LocalBasisType::Traits
              ::dimDomain,
            FieldVector<
              typename T::Traits::FiniteElementType::Traits
                ::LocalBasisType::Traits::RangeFieldType,
              T::Traits::FiniteElementType::Traits::LocalBasisType::Traits
              ::dimDomain> >,
          DiscreteGridFunctionGradient<T,X> >
    {
      typedef T GFS;
      typedef typename GFS::Traits::FiniteElementType::Traits::
        LocalBasisType::Traits LBTraits;

    public:
      typedef GridFunctionTraits<
        typename GFS::Traits::GridViewType,
        typename LBTraits::RangeFieldType,
        LBTraits::dimDomain,
        FieldVector<
          typename LBTraits::RangeFieldType,
          LBTraits::dimDomain> > Traits;

    private:
      typedef GridFunctionInterface<
        Traits,
        DiscreteGridFunctionGradient<T,X> > BaseT;

    public:
      /** \brief Construct a DiscreteGridFunctionGradient
       *
       * \param gfs The GridFunctionsSpace
       * \param x_  The coefficients vector
       */
      DiscreteGridFunctionGradient (const GFS& gfs, const X& x_)
        : pgfs(stackobject_to_shared_ptr(gfs)), xg(x_)
      { }

      // Evaluate
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        // get and bind local functions space
        typename GFS::LocalFunctionSpace lfs(*pgfs);
        lfs.bind(e);

        // get local coefficients
        std::vector<typename Traits::RangeFieldType> xl(lfs.size());
        lfs.vread(xg,xl);

        // get Jacobian of geometry
        const typename Traits::ElementType::Geometry::Jacobian&
          JgeoIT = e.geometry().jacobianInverseTransposed(x);

        // get local Jacobians/gradients of the shape functions
        std::vector<typename LBTraits::JacobianType> J(lfs.size());
        lfs.finiteElement().localBasis().evaluateJacobian(x,J);

        typename Traits::RangeType gradphi;
        y = 0;
        for(unsigned int i = 0; i < lfs.size(); ++i) {
          // compute global gradient of shape function i
          gradphi = 0;
          JgeoIT.umv(J[i][0], gradphi);

          // sum up global gradients, weighting them with the appropriate coeff
          y.axpy(xl[i], gradphi);
        }
      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return pgfs->gridview();
      }

    private:
      shared_ptr<GFS const> pgfs;
      const X& xg;
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
            typename T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
            T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimRange,
            typename T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
            >,
          DiscreteGridFunctionPiola<T,X>
          >
	{
	  typedef T GFS;

	  typedef GridFunctionInterface<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
          T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimRange,
          typename T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
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
		: pgfs(stackobject_to_shared_ptr(gfs)), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), yb(gfs.maxLocalSize())
	  {
	  }

	  inline void evaluate (const typename Traits::ElementType& e,
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
        // evaluate shape function on the reference element as before
		lfs.bind(e);
		lfs.vread(xg,xl);
        lfs.finiteElement().localBasis().evaluateFunction(x,yb);
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
	  inline const typename Traits::GridViewType& getGridView () const
	  {
		return pgfs->gridview();
	  }

	private:
	  shared_ptr<GFS const> pgfs;
	  const X& xg;
	  mutable typename GFS::LocalFunctionSpace lfs;
	  mutable std::vector<typename Traits::RangeFieldType> xl;
	  mutable std::vector<typename Traits::RangeType> yb;
	};

    /** \brief DiscreteGridFunction for vector-valued functions
     *
     * convert a power function space of scalar function spaces into a
     * vector-valued grid function this is just an intermediate
     * solution to provide VTK output
     *
     * \tparam T Type of PowerGridFunctionSpace
     * \tparam X Type of coefficients vector
     */
	template<typename T, typename X>
	class VectorDiscreteGridFunction
	  : public GridFunctionInterface<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            typename T::template Child<0>::Type::Traits::FiniteElementType
                     ::Traits::LocalBasisType::Traits::RangeFieldType,
            T::CHILDREN,
            Dune::FieldVector<
              typename T::template Child<0>::Type::Traits::FiniteElementType
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
          typename T::template Child<0>::Type::Traits::FiniteElementType
                   ::Traits::LocalBasisType::Traits::RangeFieldType,
          T::CHILDREN,
          Dune::FieldVector<
            typename T::template Child<0>::Type::Traits::FiniteElementType
                     ::Traits::LocalBasisType::Traits::RangeFieldType,
            T::CHILDREN
            >
          >,
        VectorDiscreteGridFunction<T,X>
        > BaseT;

	public:
	  typedef typename BaseT::Traits Traits;
      typedef typename T::template Child<0>::Type ChildType;
      typedef typename ChildType::Traits::FiniteElementType
                       ::Traits::LocalBasisType::Traits::RangeFieldType RF;
      typedef typename ChildType::Traits::FiniteElementType
                       ::Traits::LocalBasisType::Traits::RangeType RT;

	  VectorDiscreteGridFunction (const GFS& gfs, const X& x_)
		: pgfs(stackobject_to_shared_ptr(gfs)), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), yb(gfs.maxLocalSize())
	  {
	  }

	  inline void evaluate (const typename Traits::ElementType& e,
							const typename Traits::DomainType& x,
							typename Traits::RangeType& y) const
	  {
		lfs.bind(e);
		lfs.vread(xg,xl);
        for (unsigned int k=0; k<T::CHILDREN; k++)
          {
            lfs.child(k).finiteElement().localBasis().
              evaluateFunction(x,yb);
            y[k] = 0.0;
            for (unsigned int i=0; i<yb.size(); i++)
              y[k] += xl[lfs.child(k).localIndex(i)]*yb[i];
          }
	  }

      //! get a reference to the GridView
	  inline const typename Traits::GridViewType& getGridView () const
	  {
		return pgfs->gridview();
	  }

	private:
	  shared_ptr<GFS const> pgfs;
	  const X& xg;
	  mutable LocalFunctionSpace<GFS> lfs;
	  mutable std::vector<RF> xl;
	  mutable std::vector<RT> yb;
	};

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
