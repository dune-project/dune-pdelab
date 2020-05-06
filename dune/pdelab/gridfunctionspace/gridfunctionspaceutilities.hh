// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACEUTILITIES_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACEUTILITIES_HH

#include <cstdlib>
#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/shared_ptr.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include"../common/function.hh"
#include <dune/pdelab/common/jacobiantocurl.hh>
#include"gridfunctionspace.hh"
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

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
      : public GridFunctionBase<
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
      typedef GridFunctionBase<
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
        : pgfs(stackobject_to_shared_ptr(gfs))
        , lfs(gfs)
        , lfs_cache(lfs)
        , x_view(x_)
        , xl(gfs.maxLocalSize())
        , yb(gfs.maxLocalSize())
      {
      }

      /** \brief Construct a DiscreteGridFunction
       *
       * \param gfs shared pointer to the GridFunctionsSpace
       * \param x_  shared pointer to the coefficients vector
       */
      DiscreteGridFunction (std::shared_ptr<const GFS> gfs, std::shared_ptr<const X> x_)
        : pgfs(gfs)
        , lfs(*gfs)
        , lfs_cache(lfs)
        , x_view(*x_)
        , xl(gfs->maxLocalSize())
        , yb(gfs->maxLocalSize())
        , px(x_) // FIXME: The LocalView should handle a shared_ptr correctly!
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
        lfs_cache.update();
        x_view.bind(lfs_cache);
        x_view.read(xl);
        x_view.unbind();
        FESwitch::basis(lfs.finiteElement()).evaluateFunction(x,yb);
        y = 0;
        for (unsigned int i=0; i<yb.size(); i++)
        {
          y.axpy(xl[i],yb[i]);
        }
      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return pgfs->gridView();
      }

    private:

      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      std::shared_ptr<GFS const> pgfs;
      mutable LFS lfs;
      mutable LFSCache lfs_cache;
      mutable XView x_view;
      mutable std::vector<typename Traits::RangeFieldType> xl;
      mutable std::vector<typename Traits::RangeType> yb;
      std::shared_ptr<const X> px; // FIXME: dummy pointer to make sure we take ownership of X
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
      public GridFunctionBase<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename JacobianToCurl<typename T::Traits::FiniteElementType::
                                  Traits::LocalBasisType::Traits::JacobianType>::CurlField,
          JacobianToCurl<typename T::Traits::FiniteElementType::Traits::LocalBasisType::
                         Traits::JacobianType>::dimCurl,
          typename JacobianToCurl<typename T::Traits::FiniteElementType::
                                  Traits::LocalBasisType::Traits::JacobianType>::Curl
          >,
        DiscreteGridFunctionCurl<T,X>
        >
    {
      typedef T GFS;
      typedef typename T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::
        JacobianType Jacobian;
      typedef JacobianToCurl<Jacobian> J2C;

    public:
      typedef GridFunctionTraits<
        typename T::Traits::GridViewType,
          typename J2C::CurlField, J2C::dimCurl, typename J2C::Curl
        > Traits;

      /** \brief Construct a DiscreteGridFunctionCurl
       *
       * \param gfs  The GridFunctionsSpace
       * \param x_   The coefficients vector
       */
      DiscreteGridFunctionCurl(const GFS& gfs, const X& x_)
      : pgfs(stackobject_to_shared_ptr(gfs))
      , lfs(gfs)
      , lfs_cache(lfs)
      , x_view(x_)
      , xl(gfs.maxLocalSize())
      , jacobian(gfs.maxLocalSize())
      , yb(gfs.maxLocalSize())
      , px(stackobject_to_shared_ptr(x_))
      {}

      /** \brief Construct a DiscreteGridFunctionCurl
       *
       * \param gfs shared pointer to the GridFunctionsSpace
       * \param x_  shared pointer to the coefficients vector
       */
      DiscreteGridFunctionCurl(std::shared_ptr<const GFS> gfs, std::shared_ptr<const X> x_)
      : pgfs(gfs)
      , lfs(*gfs)
      , lfs_cache(lfs)
      , x_view(*x_)
      , xl(gfs->maxLocalSize())
      , jacobian(gfs->maxLocalSize())
      , yb(gfs->maxLocalSize())
      , px(x_)
      {}

      // Evaluate
      void evaluate (const typename Traits::ElementType& e,
                     const typename Traits::DomainType& x,
                     typename Traits::RangeType& y) const
      {
        static const J2C& j2C = J2C();

        lfs.bind();
        lfs_cache.update();
        x_view.bind(lfs_cache);
        x_view.read(xl);
        x_view.unbind();

        lfs.finiteElement().basis().evaluateJacobian(x,jacobian);

        y = 0;
        for (std::size_t i=0; i < lfs.size(); i++) {
          j2C(jacobian[i], yb);
          y.axpy(xl[i], yb);
        }
      }

      //! get a reference to the GridView
      const typename Traits::GridViewType& getGridView() const
      { return pgfs->gridView(); }

    private:
      typedef GridFunctionBase<Traits, DiscreteGridFunctionCurl<T,X> >
        BaseT;
      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      std::shared_ptr<GFS const> pgfs;
      mutable LFS lfs;
      mutable LFSCache lfs_cache;
      mutable XView x_view;
      mutable std::vector<typename Traits::RangeFieldType> xl;
      mutable std::vector<Jacobian> jacobian;
      mutable std::vector<typename Traits::RangeType> yb;
      std::shared_ptr<const X> px; // FIXME: dummy pointer to make sure we take ownership of X
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
     *       version is instantiated, static_assert() will be triggered.
     */
    template<typename GV, typename RangeFieldType, int dimRangeOfBasis>
    struct DiscreteGridFunctionCurlTraits {
      static_assert(AlwaysFalse<GV>::value,
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
      static_assert(GV::dimensionworld == 2,
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
      static_assert(GV::dimensionworld == 2,
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
      static_assert(GV::dimensionworld == 3,
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
      : public GridFunctionBase<
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
      typedef GridFunctionBase<
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
      : pgfs(stackobject_to_shared_ptr(gfs))
      , lfs(gfs)
      , lfs_cache(lfs)
      , x_view(x_)
      , xl(gfs.maxLocalSize())
      , J(gfs.maxLocalSize())
      , px(stackobject_to_shared_ptr(x_))
      {}

      /** \brief Construct a DiscreteGridFunctionGlobalCurl
       *
       * \param gfs shared pointer to the GridFunctionsSpace
       * \param x_  shared pointer to the coefficients vector
       */
      DiscreteGridFunctionGlobalCurl(std::shared_ptr<const GFS> gfs, std::shared_ptr<const X> x_)
      : pgfs(gfs)
      , lfs(*gfs)
      , lfs_cache(lfs)
      , x_view(*x_)
      , xl(gfs->maxLocalSize())
      , J(gfs->maxLocalSize())
      , px(x_)
      {}

      // Evaluate
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        lfs.bind(e);
        lfs_cache.update();
        x_view.bind(lfs_cache);
        x_view.read(xl);
        x_view.unbind();

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
            //how did that pass all the static asserts?
            std::abort();
          }
      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return pgfs->gridView();
      }

    private:
      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      std::shared_ptr<GFS const> pgfs;
      mutable LFS lfs;
      mutable LFSCache lfs_cache;
      mutable XView x_view;
      mutable std::vector<typename Traits::RangeFieldType> xl;
      mutable std::vector<typename T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType> J;
      std::shared_ptr<const X> px; // FIXME: dummy pointer to make sure we take ownership of X
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
      : public GridFunctionBase<
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
      typedef GridFunctionBase<
        Traits,
        DiscreteGridFunctionGradient<T,X> > BaseT;

    public:
      /** \brief Construct a DiscreteGridFunctionGradient
       *
       * \param gfs The GridFunctionsSpace
       * \param x_  The coefficients vector
       */
      DiscreteGridFunctionGradient (const GFS& gfs, const X& x_)
        : pgfs(stackobject_to_shared_ptr(gfs))
        , lfs(gfs)
        , lfs_cache(lfs)
        , x_view(x_)
        , xl(lfs.size())
      { }

      /** \brief Construct a DiscreteGridFunctionGradient
       *
       * \param gfs shared pointer to the GridFunctionsSpace
       * \param x_  shared pointer to the coefficients vector
       */
      DiscreteGridFunctionGradient (std::shared_ptr<const GFS> gfs, std::shared_ptr<const X> x_)
        : pgfs(gfs)
        , lfs(*gfs)
        , lfs_cache(lfs)
        , x_view(*x_)
        , xl(lfs.size())
      { }

      // Evaluate
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        // get and bind local functions space
        lfs.bind(e);
        lfs_cache.update();
        x_view.bind(lfs_cache);

        // get local coefficients
        xl.resize(lfs.size());
        x_view.read(xl);
        x_view.unbind();

        // get Jacobian of geometry
        const typename Traits::ElementType::Geometry::JacobianInverseTransposed
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
        return pgfs->gridView();
      }

    private:
      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      std::shared_ptr<GFS const> pgfs;
      mutable LFS lfs;
      mutable LFSCache lfs_cache;
      mutable XView x_view;
      mutable std::vector<typename Traits::RangeFieldType> xl;
    };

    /** \brief DiscreteGridFunction with Piola transformation
     *
     * \copydetails DiscreteGridFunction
     */
    template<typename T, typename X>
    class DiscreteGridFunctionPiola
      : public GridFunctionBase<
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

      typedef GridFunctionBase<
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
        : pgfs(stackobject_to_shared_ptr(gfs))
        , lfs(gfs)
        , lfs_cache(lfs)
        , x_view(x_)
        , xl(pgfs->maxLocalSize())
        , yb(pgfs->maxLocalSize())
      {
      }

      /** \brief Construct a DiscreteGridFunctionPiola
       *
       * \param gfs shared pointer to the GridFunctionsSpace
       * \param x_  shared pointer to the coefficients vector
       */
      DiscreteGridFunctionPiola (std::shared_ptr<const GFS> gfs, std::shared_ptr<const X> x_)
        : pgfs(gfs)
        , lfs(*gfs)
        , lfs_cache(lfs)
        , x_view(*x_)
        , xl(pgfs->maxLocalSize())
        , yb(pgfs->maxLocalSize())
      {
      }

      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        // evaluate shape function on the reference element as before
        lfs.bind(e);
        lfs_cache.update();
        x_view.bind(lfs_cache);
        x_view.read(xl);
        x_view.unbind();

        lfs.finiteElement().localBasis().evaluateFunction(x,yb);
        typename Traits::RangeType yhat;
        yhat = 0;
        for (unsigned int i=0; i<yb.size(); i++)
          yhat.axpy(xl[i],yb[i]);

        // apply Piola transformation
        typename Traits::ElementType::Geometry::JacobianInverseTransposed
          J = e.geometry().jacobianInverseTransposed(x);
        J.invert();
        y = 0;
        J.umtv(yhat,y);
        y /= J.determinant();
      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return pgfs->gridView();
      }

    private:

      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      std::shared_ptr<GFS const> pgfs;
      mutable LFS lfs;
      mutable LFSCache lfs_cache;
      mutable XView x_view;
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
     * \tparam dimR Force a different number of components for the resulting
     *              GridFunction than the PowerGridFunctionSpace.
     */
#ifdef __clang__
    // clang is too stupid to correctly apply the constexpr qualifier of staticDegree in this context
    template<typename T, typename X, std::size_t dimR = TypeTree::StaticDegree<T>::value>
#else
    template<typename T, typename X, std::size_t dimR = TypeTree::StaticDegree<T>::value>
#endif
    class VectorDiscreteGridFunction
      : public GridFunctionBase<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            typename T::template Child<0>::Type::Traits::FiniteElementType
                     ::Traits::LocalBasisType::Traits::RangeFieldType,
            dimR,
            Dune::FieldVector<
              typename T::template Child<0>::Type::Traits::FiniteElementType
                       ::Traits::LocalBasisType::Traits::RangeFieldType,
              dimR
              >
            >,
          VectorDiscreteGridFunction<T,X>
          >
    {
      typedef T GFS;

      typedef GridFunctionBase<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::template Child<0>::Type::Traits::FiniteElementType
                   ::Traits::LocalBasisType::Traits::RangeFieldType,
          dimR,
          Dune::FieldVector<
            typename T::template Child<0>::Type::Traits::FiniteElementType
                     ::Traits::LocalBasisType::Traits::RangeFieldType,
            dimR
            >
          >,
        VectorDiscreteGridFunction<T,X,dimR>
        > BaseT;

    public:
      typedef typename BaseT::Traits Traits;
      typedef typename T::template Child<0>::Type ChildType;
      typedef typename ChildType::Traits::FiniteElementType
                       ::Traits::LocalBasisType::Traits::RangeFieldType RF;
      typedef typename ChildType::Traits::FiniteElementType
                       ::Traits::LocalBasisType::Traits::RangeType RT;

      //! construct
      /**
       * \param gfs   GridFunctionSpace.
       * \param x_    Coefficient vector.
       * \param start Number of first child of gfs to use.
       */
      VectorDiscreteGridFunction(const GFS& gfs, const X& x_,
                                 std::size_t start = 0)
      : pgfs(stackobject_to_shared_ptr(gfs))
      , lfs(gfs)
      , lfs_cache(lfs)
      , x_view(x_)
      , xl(gfs.maxLocalSize())
      , yb(gfs.maxLocalSize())
      {
        for(std::size_t i = 0; i < dimR; ++i)
          remap[i] = i + start;
      }

      /** \brief Construct a VectorDiscreteGridFunction
       *
       * \param gfs shared pointer to the GridFunctionsSpace
       * \param x_  shared pointer to the coefficients vector
       * \param start Number of first child of gfs to use.
       */
      VectorDiscreteGridFunction (std::shared_ptr<const GFS> gfs, std::shared_ptr<const X> x_,
                                 std::size_t start = 0)
        : pgfs(gfs)
        , lfs(*gfs)
        , lfs_cache(lfs)
        , x_view(*x_)
        , xl(pgfs->maxLocalSize())
        , yb(pgfs->maxLocalSize())
      {
         for(std::size_t i = 0; i < dimR; ++i)
            remap[i] = i + start;
      }

      //! construct
      /**
       * \param gfs    GridFunctionSpace.
       * \param x_     Coefficient vector.
       * \param remap_ Subscriptable entity (i.e. a container, array, or
       *               pointer) with at least dimR entries.  The relevant
       *               entries are copied.
       *
       * \note If \c i denotes a component of the resulting grid function,
       *       then remap_[i] denotes the corresponding child of the
       *       gridfunctionspace.
       */
      template<class Remap>
      VectorDiscreteGridFunction(const GFS& gfs, const X& x_,
                                 const Remap &remap_)
      : pgfs(stackobject_to_shared_ptr(gfs))
      , lfs(gfs)
      , lfs_cache(lfs)
      , x_view(x_)
      , xl(gfs.maxLocalSize())
      , yb(gfs.maxLocalSize())
      ,	px(stackobject_to_shared_ptr(x_))
      {
        for(std::size_t i = 0; i < dimR; ++i)
          remap[i] = remap_[i];
      }

      /** \brief Construct a VectorDiscreteGridFunction
       *
       * \param gfs shared pointer to the GridFunctionsSpace
       * \param x_  shared pointer to the coefficients vector
       * \param remap_ Subscriptable entity (i.e. a container, array, or
       *               pointer) with at least dimR entries.  The relevant
       *               entries are copied.
       *
       * \note If \c i denotes a component of the resulting grid function,
       *       then remap_[i] denotes the corresponding child of the
       *       gridfunctionspace.
       */
      template<class Remap>
      VectorDiscreteGridFunction (std::shared_ptr<const GFS> gfs, std::shared_ptr<const X> x_,
                                 const Remap &remap_)
        : pgfs(gfs)
        , lfs(*gfs)
        , lfs_cache(lfs)
        , x_view(*x_)
        , xl(pgfs->maxLocalSize())
        , yb(pgfs->maxLocalSize())
        , px(x_)
      {
        for(std::size_t i = 0; i < dimR; ++i)
          remap[i] = remap_[i];
      }

      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        lfs.bind(e);
        lfs_cache.update();
        x_view.bind(lfs_cache);
        x_view.read(xl);
        x_view.unbind();
        for (unsigned int k=0; k < dimR; k++)
          {
            lfs.child(remap[k]).finiteElement().localBasis().
              evaluateFunction(x,yb);
            y[k] = 0.0;
            for (unsigned int i=0; i<yb.size(); i++)
              y[k] += xl[lfs.child(remap[k]).localIndex(i)]*yb[i];
          }
      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return pgfs->gridView();
      }

    private:
      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      std::shared_ptr<GFS const> pgfs;
      std::size_t remap[dimR];
      mutable LFS lfs;
      mutable LFSCache lfs_cache;
      mutable XView x_view;
      mutable std::vector<RF> xl;
      mutable std::vector<RT> yb;
      std::shared_ptr<const X> px; // FIXME: dummy pointer to make sure we take ownership of X
    };

    /** \brief Equivalent of DiscreteGridFunctionGradient for vector-valued functions
     *
     * \tparam T Type of PowerGridFunctionSpace
     * \tparam X Type of coefficients vector
     */
    template<typename T, typename X>
    class VectorDiscreteGridFunctionGradient
      : public GridFunctionBase<
                 GridFunctionTraits<
                   typename T::Traits::GridViewType,
                   typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
                   //T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimDomain,
                   TypeTree::StaticDegree<T>::value,
                   Dune::FieldMatrix<
                     typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
                     TypeTree::StaticDegree<T>::value,
                     T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimDomain
                   >
                 >,
                 VectorDiscreteGridFunctionGradient<T,X>
               >
    {
      typedef T GFS;

      typedef GridFunctionBase<
                GridFunctionTraits<
                  typename T::Traits::GridViewType,
                  typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
                  //T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimDomain,
                  TypeTree::StaticDegree<T>::value,
                  Dune::FieldMatrix<
                    typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
                    TypeTree::StaticDegree<T>::value,
                    T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimDomain>
                  >,
                  VectorDiscreteGridFunctionGradient<T,X>
                > BaseT;

    public:
      typedef typename BaseT::Traits Traits;
      typedef typename T::template Child<0>::Type ChildType;
      typedef typename ChildType::Traits::FiniteElementType::Traits::LocalBasisType::Traits LBTraits;

      typedef typename LBTraits::RangeFieldType RF;
      typedef typename LBTraits::JacobianType JT;

      /** \brief Construct a VectorDiscreteGridFunctionGradient
       *
       * \param gfs    GridFunctionSpace.
       * \param x_     Coefficient vector.
       */
      VectorDiscreteGridFunctionGradient (const GFS& gfs, const X& x_)
        : pgfs(stackobject_to_shared_ptr(gfs))
        , lfs(gfs)
        , lfs_cache(lfs)
        , x_view(x_)
        , xl(gfs.maxLocalSize())
        , J(gfs.maxLocalSize())
      {
      }

      /** \brief Construct a VectorDiscreteGridFunctionGradient
       *
       * \param gfs shared pointer to the GridFunctionsSpace
       * \param x_  shared pointer to the coefficients vector
       */
      VectorDiscreteGridFunctionGradient (std::shared_ptr<const GFS> gfs, std::shared_ptr<const X> x_)
        : pgfs(gfs)
        , lfs(*gfs)
        , lfs_cache(lfs)
        , x_view(*x_)
        , xl(pgfs->maxLocalSize())
        , J(pgfs->maxLocalSize())
      {
      }

      inline void evaluate(const typename Traits::ElementType& e,
          const typename Traits::DomainType& x,
          typename Traits::RangeType& y) const
      {
        // get and bind local functions space
        lfs.bind(e);
        lfs_cache.update();
        x_view.bind(lfs_cache);
        x_view.read(xl);
        x_view.unbind();

        // get Jacobian of geometry
        const typename Traits::ElementType::Geometry::JacobianInverseTransposed
          JgeoIT = e.geometry().jacobianInverseTransposed(x);

        y = 0.0;

        // Loop over PowerLFS and calculate gradient for each child separately
        for(unsigned int k = 0; k != TypeTree::degree(lfs); ++k)
        {
          // get local Jacobians/gradients of the shape functions
          std::vector<typename LBTraits::JacobianType> J(lfs.child(k).size());
          lfs.child(k).finiteElement().localBasis().evaluateJacobian(x,J);

          Dune::FieldVector<RF,LBTraits::dimDomain> gradphi;
          for (typename LFS::Traits::SizeType i=0; i<lfs.child(k).size(); i++)
          {
            gradphi = 0;
            JgeoIT.umv(J[i][0], gradphi);

            y[k].axpy(xl[lfs.child(k).localIndex(i)], gradphi);
          }
        }
      }


      //! \brief get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return pgfs->gridView();
      }

    private:
      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      std::shared_ptr<GFS const> pgfs;
      mutable LFS lfs;
      mutable LFSCache lfs_cache;
      mutable XView x_view;
      mutable std::vector<RF> xl;
      mutable std::vector<JT> J;
      std::shared_ptr<const X> px; // FIXME: dummy pointer to make sure we take ownership of X
    };

    /** \brief Helper class to compute a single derivative of scalar basis functions.

        \tparam Mat  The matrix type of the geometry transformation.
        \tparam RF   The underlying field type of Mat,
                     needed especially for the template specialization.
        \tparam size The size of Mat,
                     needed especially for the template specialization

        * The k-th derivative of scalar basis functions is calculated from the
        * matrix-vector product of the geometry transformation and the gradient
        * on the reference element by picking out the k-th component.

    */
    template<typename Mat, typename RF, std::size_t size>
    struct SingleDerivativeComputationHelper {
      /**
       * \param[in] mat The jacobian of the geometry transformation.
       * \param[in] t   The gradient of the shape function on the reference element.
      */
      template<typename T>
      static inline RF compute_derivative(const Mat& mat, const T& t, const unsigned int k)
      {
        // setting this to zero is just a test if the template specialization work
        Dune::FieldVector<RF,size> grad_phi(0.0);
        mat.umv(t,grad_phi);
        return grad_phi[k];
        // return 0.0;
      }
    };

    /** \brief Template specialization for Dune::FieldMatrix.
     *
     * This is a template specialization if the matrix type of the geometry transformation
     * is equal to a Dune::FieldMatrix. The k-th component of the matrix-vector product
     * can then be calculated as a scalar product of the k-th row of the geometry
     * transformation with the gradient on the reference element.
     *
     */
    template<typename RF, std::size_t size>
    struct SingleDerivativeComputationHelper<Dune::FieldMatrix<RF,size,size>,RF,size> {
      /**
       * \param[in] mat The jacobian of the geometry transformation,
       *                has to be a Dune::FieldMatrix.
       * \param[in] t   The gradient of the shape function on the reference element.
       */
      template<typename T>
      static inline RF compute_derivative(const Dune::FieldMatrix<RF,size,size>& mat, const T& t, const unsigned int k)
      {
        return mat[k]*t;
      }
    };

    /** \brief Template specialization for Dune::DiagonalMatrix.
     *
     * This is a template specialization if the matrix type of the geometry transformation
     * is equal to a Dune::DiagonalMatrix. The k-th component of the matrix-vector product
     * can then be calculated as the product of the k-th diagonal element of the
     * geometry transformation with the k-th derivative of the gradient on the reference
     * element.
     * This specialization especially occurs for YaspGrid.
     *
     */
    template<typename RF, std::size_t size>
    struct SingleDerivativeComputationHelper<Dune::DiagonalMatrix<RF,size>,RF,size> {
      /**
       * \param[in] mat The jacobian of the geometry transformation,
       *                has to be a Dune::DiagonalMatrix.
       * \param[in] t   The gradient of the shape function on the reference element.
       */
      template<typename T>
      static inline RF compute_derivative(const Dune::DiagonalMatrix<RF,size>& mat, const T& t, const unsigned int k)
      {
        return mat[k][k]*t[k];
      }
    };

    /** \brief Compute divergence of vector-valued functions.

        \tparam T Type of VectorGridFunctionSpace.
        \tparam X Type of coefficients vector.

        * The grid function space should be a vector grid function space
        * consisting of scalar-valued function spaces. The divergence will be
        * a single component function.

    */
    template<typename T, typename X>
    class VectorDiscreteGridFunctionDiv
      : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<
        typename T::Traits::GridViewType,
        typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
        T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimRange,
        typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType>,
      VectorDiscreteGridFunctionDiv<T,X> >
    {
      typedef T GFS;

      typedef Dune::PDELab::GridFunctionBase<
        Dune::PDELab::GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
          T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimRange,
          typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType>,
        VectorDiscreteGridFunctionDiv<T,X> > BaseT;
    public :
      typedef typename BaseT::Traits Traits;
      typedef typename T::template Child<0>::Type ChildType;
      typedef typename ChildType::Traits::FiniteElementType::Traits::LocalBasisType::Traits LBTraits;

      typedef typename LBTraits::RangeFieldType RF;
      typedef typename LBTraits::JacobianType JT;

      /** \brief Construct a VectorDiscreteGridFunctionDiv
       *
       * \param gfs    GridFunctionSpace.
       * \param x_     Coefficient vector.
       */
      VectorDiscreteGridFunctionDiv(const GFS& gfs, const X& x_)
        : pgfs(stackobject_to_shared_ptr(gfs))
        , lfs(gfs)
        , lfs_cache(lfs)
        , x_view(x_)
        , xl(gfs.maxLocalSize())
        , J(gfs.maxLocalSize())
      {
        static_assert(LBTraits::dimDomain == TypeTree::StaticDegree<T>::value,
                           "dimDomain and number of children has to be the same");
      }

      /** \brief Construct a VectorDiscreteGridFunctionDiv
       *
       * \param gfs shared pointer to the GridFunctionsSpace
       * \param x_  shared pointer to the coefficients vector
       */
      VectorDiscreteGridFunctionDiv (std::shared_ptr<const GFS> gfs, std::shared_ptr<const X> x_)
        : pgfs(gfs)
        , lfs(*gfs)
        , lfs_cache(lfs)
        , x_view(*x_)
        , xl(pgfs->maxLocalSize())
        , J(pgfs->maxLocalSize())
      {
        static_assert(LBTraits::dimDomain == TypeTree::StaticDegree<T>::value,
                           "dimDomain and number of children has to be the same");
      }

      inline void evaluate(const typename Traits::ElementType& e,
                           const typename Traits::DomainType& x,
                           typename Traits::RangeType& y) const
      {
        // get and bind local functions space
        lfs.bind(e);
        lfs_cache.update();
        x_view.bind(lfs_cache);
        x_view.read(xl);
        x_view.unbind();

        // get Jacobian of geometry
        const typename Traits::ElementType::Geometry::JacobianInverseTransposed
          JgeoIT = e.geometry().jacobianInverseTransposed(x);

        const typename Traits::ElementType::Geometry::JacobianInverseTransposed::size_type N =
          Traits::ElementType::Geometry::JacobianInverseTransposed::rows;

        y = 0.0;

        // loop over VectorGFS and calculate k-th derivative of k-th child
        for(unsigned int k=0; k != TypeTree::degree(lfs); ++k) {

          // get local Jacobians/gradients of the shape functions
          std::vector<typename LBTraits::JacobianType> J(lfs.child(k).size());
          lfs.child(k).finiteElement().localBasis().evaluateJacobian(x,J);

          RF d_k_phi;
          for(typename LFS::Traits::SizeType i=0; i<lfs.child(k).size(); i++) {
            // compute k-th derivative of k-th child
            d_k_phi =
              SingleDerivativeComputationHelper<
                typename Traits::ElementType::Geometry::JacobianInverseTransposed,
              typename Traits::ElementType::Geometry::JacobianInverseTransposed::field_type,
              N>::template compute_derivative<typename LBTraits::JacobianType::row_type>
              (JgeoIT,J[i][0],k);

            y += xl[lfs.child(k).localIndex(i)] * d_k_phi;
          }
        }
      }

      //! \brief get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView() const
      {
        return pgfs->gridView();
      }

    private :
      typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
      typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      std::shared_ptr<GFS const> pgfs;
      mutable LFS lfs;
      mutable LFSCache lfs_cache;
      mutable XView x_view;
      mutable std::vector<RF> xl;
      mutable std::vector<JT> J;
      std::shared_ptr<const X> px;
    }; // end class VectorDiscreteGridFunctionDiv

    /** \brief Compute curl of vector-valued functions.

        \tparam T Type of VectorGridFunctionSpace.
        \tparam X Type of coefficients vector.
        \tparam dimR The number of components the curl should be taken of.

        \note This is the non-specialized version of VectorDiscreteGridFunctionCurl.
              There is a specialized version for the values of dimR equal to 2 or 3.
              If this non-specialized version is instantiated, a static_assert()
              will be triggered.

     */
    template<typename T, typename X, std::size_t dimR = TypeTree::StaticDegree<T>::value>
    class VectorDiscreteGridFunctionCurl
    {
      typedef T GFS;
    public :

      VectorDiscreteGridFunctionCurl(const GFS& gfs, const X& x)
      {
        static_assert(AlwaysFalse<typename GFS::Traits::GridViewType>::value,
                      "Curl computation can only be done in two or three dimensions");
      }

      VectorDiscreteGridFunctionCurl(std::shared_ptr<const GFS> gfs, std::shared_ptr<const X> x_)
      {
        static_assert(AlwaysFalse<typename GFS::Traits::GridViewType>::value,
                      "Curl computation can only be done in two or three dimensions");
      }

    };

    /** \brief Compute curl of vector-valued functions (3D).

        \tparam T Type of VectorGridFunctionSpace.
        \tparam X Type of coefficients vector.

        * This is the specialized version for dimR == 3. It takes the curl of a
        * 3D-valued function and the result will be a 3-component vector
        * consisting of scalar-valued functions.

     */
    template<typename T, typename X>
    class VectorDiscreteGridFunctionCurl<T,X,3>
      : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<
        typename T::Traits::GridViewType,
        typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
        //T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimRange,
        TypeTree::StaticDegree<T>::value,
        Dune::FieldVector<
          typename T::template Child<0>::Type::Traits::FiniteElementType
          ::Traits::LocalBasisType::Traits::RangeFieldType,
          //T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimRange
          TypeTree::StaticDegree<T>::value
          >
        >,
      VectorDiscreteGridFunctionCurl<T,X>
      >
    {
      typedef T GFS;

      typedef Dune::PDELab::GridFunctionBase<
        Dune::PDELab::GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
          //T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimRange,
          TypeTree::StaticDegree<T>::value,
          Dune::FieldVector<
            typename T::template Child<0>::Type::Traits::FiniteElementType
            ::Traits::LocalBasisType::Traits::RangeFieldType,
            //T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimRange
            TypeTree::StaticDegree<T>::value
            >
          >,
        VectorDiscreteGridFunctionCurl<T,X> > BaseT;

    public :
      typedef typename BaseT::Traits Traits;
      typedef typename T::template Child<0>::Type ChildType;
      typedef typename ChildType::Traits::FiniteElementType::Traits::LocalBasisType::Traits LBTraits;

      typedef typename LBTraits::RangeFieldType RF;
      typedef typename LBTraits::JacobianType JT;

      /** \brief Construct a VectorDiscreteGridFunctionCurl
       *
       * \param gfs    GridFunctionSpace.
       * \param x_     Coefficient vector.
       */
      VectorDiscreteGridFunctionCurl(const GFS& gfs, const X& x_)
        : pgfs(stackobject_to_shared_ptr(gfs))
        , lfs(gfs)
        , lfs_cache(lfs)
        , x_view(x_)
        , xl(gfs.maxLocalSize())
        , J(gfs.maxLocalSize())
      {
        static_assert(LBTraits::dimDomain == TypeTree::StaticDegree<T>::value,
                           "dimDomain and number of children has to be the same");
      }

      /** \brief Construct a VectorDiscreteGridFunctionCurl
       *
       * \param gfs shared pointer to the GridFunctionsSpace
       * \param x_  shared pointer to the coefficients vector
       */
      VectorDiscreteGridFunctionCurl (std::shared_ptr<const GFS> gfs, std::shared_ptr<const X> x_)
        : pgfs(gfs)
        , lfs(*gfs)
        , lfs_cache(lfs)
        , x_view(*x_)
        , xl(pgfs->maxLocalSize())
        , J(pgfs->maxLocalSize())
      {
        static_assert(LBTraits::dimDomain == TypeTree::StaticDegree<T>::value,
                           "dimDomain and number of children has to be the same");
      }

      inline void evaluate(const typename Traits::ElementType& e,
                           const typename Traits::DomainType& x,
                           typename Traits::RangeType& y) const
      {
        // get and bind local functions space
        lfs.bind(e);
        lfs_cache.update();
        x_view.bind(lfs_cache);
        x_view.read(xl);
        x_view.unbind();

        // get Jacobian of geometry
        const typename Traits::ElementType::Geometry::JacobianInverseTransposed
          JgeoIT = e.geometry().jacobianInverseTransposed(x);

        const typename Traits::ElementType::Geometry::JacobianInverseTransposed::size_type N =
          Traits::ElementType::Geometry::JacobianInverseTransposed::rows;

        y = 0.0;

        // some handy variables for the curl in 3D
        int i1, i2;

        // loop over childs of VectorGFS
        for(unsigned int k=0; k != TypeTree::degree(lfs); ++k) {

          // get local Jacobians/gradients of the shape functions
          std::vector<typename LBTraits::JacobianType> J(lfs.child(k).size());
          lfs.child(k).finiteElement().localBasis().evaluateJacobian(x,J);

          // pick out the right derivative and component for the curl computation
          i1 = (k+1)%3;
          i2 = (k+2)%3;

          RF d_k_phi;
          for(typename LFS::Traits::SizeType i=0; i<lfs.child(k).size(); i++) {
            // compute i2-th derivative of k-th child
            d_k_phi =
              SingleDerivativeComputationHelper<
                typename Traits::ElementType::Geometry::JacobianInverseTransposed,
              typename Traits::ElementType::Geometry::JacobianInverseTransposed::field_type,
              N>::template compute_derivative<typename LBTraits::JacobianType::row_type>
              (JgeoIT,J[i][0],i2);

            y[i1] += xl[lfs.child(k).localIndex(i)] * d_k_phi;

            // compute i1-th derivative of k-th child
            d_k_phi =
              SingleDerivativeComputationHelper<
                typename Traits::ElementType::Geometry::JacobianInverseTransposed,
              typename Traits::ElementType::Geometry::JacobianInverseTransposed::field_type,
              N>::template compute_derivative<typename LBTraits::JacobianType::row_type>
              (JgeoIT,J[i][0],i1);

            y[i2] -= xl[lfs.child(k).localIndex(i)] * d_k_phi;
          }
        }
      }

      //! \brief get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView() const
      {
        return pgfs->gridView();
      }

    private :
      typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
      typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      std::shared_ptr<GFS const> pgfs;
      mutable LFS lfs;
      mutable LFSCache lfs_cache;
      mutable XView x_view;
      mutable std::vector<RF> xl;
      mutable std::vector<JT> J;
      std::shared_ptr<const X> px;
    }; // end class VectorDiscreteGridFunctionCurl (3D)

    /** \brief Compute curl of vector-valued functions (2D).

        \tparam T Type of VectorGridFunctionSpace.
        \tparam X Type of coefficients vector.

        * This is the specialized version for dimR == 2. It takes the curl of a
        * 2D-valued function and the result will be a single component
        * scalar-valued function.

     */
    template<typename T, typename X>
    class VectorDiscreteGridFunctionCurl<T,X,2>
      : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<
        typename T::Traits::GridViewType,
        typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
        T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimRange,
        typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType>,
      VectorDiscreteGridFunctionDiv<T,X> >
    {
      typedef T GFS;

      typedef Dune::PDELab::GridFunctionBase<
        Dune::PDELab::GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
          T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimRange,
          typename T::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType>,
        VectorDiscreteGridFunctionDiv<T,X> > BaseT;
    public :
      typedef typename BaseT::Traits Traits;
      typedef typename T::template Child<0>::Type ChildType;
      typedef typename ChildType::Traits::FiniteElementType::Traits::LocalBasisType::Traits LBTraits;

      typedef typename LBTraits::RangeFieldType RF;
      typedef typename LBTraits::JacobianType JT;

      /** \brief Construct a VectorDiscreteGridFunctionCurl
       *
       * \param gfs    GridFunctionSpace.
       * \param x_     Coefficient vector.
       */
      VectorDiscreteGridFunctionCurl(const GFS& gfs, const X& x_)
        : pgfs(stackobject_to_shared_ptr(gfs))
        , lfs(gfs)
        , lfs_cache(lfs)
        , x_view(x_)
        , xl(gfs.maxLocalSize())
        , J(gfs.maxLocalSize())
      {
        static_assert(LBTraits::dimDomain == TypeTree::StaticDegree<T>::value,
                           "dimDomain and number of children has to be the same");
      }

      /** \brief Construct a VectorDiscreteGridFunctionCurl
       *
       * \param gfs shared pointer to the GridFunctionsSpace
       * \param x_  shared pointer to the coefficients vector
       */
      VectorDiscreteGridFunctionCurl (std::shared_ptr<const GFS> gfs, std::shared_ptr<const X> x_)
        : pgfs(gfs)
        , lfs(*gfs)
        , lfs_cache(lfs)
        , x_view(*x_)
        , xl(pgfs->maxLocalSize())
        , J(pgfs->maxLocalSize())
      {
        static_assert(LBTraits::dimDomain == TypeTree::StaticDegree<T>::value,
                           "dimDomain and number of children has to be the same");
      }

      inline void evaluate(const typename Traits::ElementType& e,
                           const typename Traits::DomainType& x,
                           typename Traits::RangeType& y) const
      {
        // get and bind local functions space
        lfs.bind(e);
        lfs_cache.update();
        x_view.bind(lfs_cache);
        x_view.read(xl);
        x_view.unbind();

        // get Jacobian of geometry
        const typename Traits::ElementType::Geometry::JacobianInverseTransposed
          JgeoIT = e.geometry().jacobianInverseTransposed(x);

        const typename Traits::ElementType::Geometry::JacobianInverseTransposed::size_type N =
          Traits::ElementType::Geometry::JacobianInverseTransposed::rows;

        y = 0.0;

        // some handy variables for the curl computation in 2D
        RF sign = -1.0;
        int i2;

        // loop over childs of VectorGFS
        for(unsigned int k=0; k != TypeTree::degree(lfs); ++k) {

          // get local Jacobians/gradients of the shape functions
          std::vector<typename LBTraits::JacobianType> J(lfs.child(k).size());
          lfs.child(k).finiteElement().localBasis().evaluateJacobian(x,J);

          RF d_k_phi;
          // pick out the right derivative
          i2 = 1-k;
          for(typename LFS::Traits::SizeType i=0; i<lfs.child(k).size(); i++) {
            // compute i2-th derivative of k-th child
            d_k_phi =
              SingleDerivativeComputationHelper<
                typename Traits::ElementType::Geometry::JacobianInverseTransposed,
              typename Traits::ElementType::Geometry::JacobianInverseTransposed::field_type,
              N>::template compute_derivative<typename LBTraits::JacobianType::row_type>
              (JgeoIT,J[i][0],i2);

            y += sign * xl[lfs.child(k).localIndex(i)] * d_k_phi;
          }
          sign *= -1.0;
        }
      }

      //! \brief get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView() const
      {
        return pgfs->gridView();
      }

    private :
      typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
      typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      std::shared_ptr<GFS const> pgfs;
      mutable LFS lfs;
      mutable LFSCache lfs_cache;
      mutable XView x_view;
      mutable std::vector<RF> xl;
      mutable std::vector<JT> J;
      std::shared_ptr<const X> px;
    }; // end class VectorDiscreteGridFunctionCurl (2D)

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACEUTILITIES_HH
