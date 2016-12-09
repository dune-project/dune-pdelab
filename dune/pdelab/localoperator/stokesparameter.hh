#ifndef DUNE_PDELAB_LOCALOPERATOR_STOKESPARAMETER_HH
#define DUNE_PDELAB_LOCALOPERATOR_STOKESPARAMETER_HH

#include <dune/common/parametertree.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>

namespace Dune {
  namespace PDELab {

    /**
     * These are the boundary condition types as to be returned by
     * the employed boundary type function.
     *
     * Possible types:
     *
     * <ul>
     *
     * <li>\a DoNothing : No boundary conditions.
     *
     * <li>\a VelocityDirichlet : Dirichlet conditions for velocity.
     *
     * <li>\a StressNeumann : Natural Neumann conditions for the
     * impulse flux. These are equivalent to a fixed pressure
     * condition \b if \f$ \forall i : n \cdot \nabla v_i = 0 \f$.
     *
     * <li>\a SlipVelocity : Smooth transition between slip and no-slip
     * condition - only works for DG!
     *
     * </ul>
     */
    struct StokesBoundaryCondition {
      enum Type {
        DoNothing = 0,
        VelocityDirichlet = 1,
        StressNeumann = 2,
        SlipVelocity = 3
      };
    };

    /**
     * Traits class for the parameter class of a Navier-Stokes local operator.
     */
    template<typename GV, typename RF>
    struct NavierStokesParameterTraits
    {
      //! \brief the grid view
      typedef GV GridView;

      //! \brief Enum for domain dimension
      enum {
        //! \brief dimension of the domain
        dimDomain = GV::dimension
      };

      //! \brief Export type for domain field
      typedef typename GV::Grid::ctype DomainField;

      //! \brief domain type
      typedef Dune::FieldVector<DomainField,dimDomain> Domain;

      //! \brief domain type
      typedef Dune::FieldVector<DomainField,dimDomain-1> IntersectionDomain;

      //! \brief Export type for range field
      typedef RF RangeField;

      //! \brief deformation range type
      typedef Dune::FieldVector<RF,GV::dimensionworld> VelocityRange;

      //! \brief pressure range type
      typedef Dune::FieldVector<RF,1> PressureRange;

      //! \brief boundary type value
      typedef StokesBoundaryCondition BoundaryCondition;

      //! \brief grid element type
      typedef typename GV::Traits::template Codim<0>::Entity Element;
      //! \brief grid intersection type
      typedef typename GV::Intersection Intersection;
    };

    namespace {

      /**
       * Compile-time switch allowing the evaluation of a
       * vector valued grid function.
       */
      template<typename GF, typename Entity, typename Domain>
      typename GF::Traits::RangeType
      evaluateVelocityGridFunction(const GF& gf,
                                   const Entity& e,
                                   const Domain& x)
      {
        static_assert(int(GF::Traits::dimRange) == int(Domain::dimension),"dimension of function range does not match grid dimension");
        typename GF::Traits::RangeType y;
        gf.evaluate(e,x,y);
        return y;
      }

      /**
       * Compile-time switch allowing the evaluation of a
       * power grid function.
       */
      template<typename GF, typename Entity, typename Domain>
      FieldVector<typename GF::template Child<0>::Type::Traits::RangeFieldType,TypeTree::StaticDegree<GF>::value>
        evaluateVelocityGridFunction(const GF& gf,
                                     const Entity& e,
                                     const Domain& x)
      {
        static_assert(Domain::dimension == TypeTree::StaticDegree<GF>::value,"dimension of function range does not match grid dimension");
        FieldVector<typename GF::template Child<0>::Type::Traits::RangeFieldType,TypeTree::StaticDegree<GF>::value> y;
        typename GF::template Child<0>::Type::Traits::RangeType cy;
        for (int d = 0; d < Domain::dimension; ++d)
          {
            gf.child(d).evaluate(e,x,cy);
            y[d] = cy;
          }
        return y;
      }

    }

    /**
     * Default implementation for the parameter class to be
     * used with the Taylor-Hood Navier-Stokes local operator.
     *
     * This is designed to work with the TaylorHoodNavierStokes,
     * TaylorHoodNavierStokesJacobian and NavierStokesMass local
     * operator classes.
     *
     * \tparam GV    GridView.
     * \tparam RF    The range field type of the Navier-Stokes solution.
     * \tparam F     External force term function (vector-valued).
     * \tparam B     Boundary type function returning an element
     *               of StokesBoundaryCondition.
     * \tparam V     Dirichlet velocity function.
     * \tparam J     Neumann stress boundary function (vector- or scalar-valued).
     *               Scalar values will be interpreted as the magnitude of a vector
     *               oriented in outer normal direction.
     */
    template <typename GV, typename RF, typename F, typename B, typename V, typename J, bool navier = false, bool tensor = false>
    class NavierStokesDefaultParameters
    {
    public:

      static const bool assemble_navier = navier;
      static const bool assemble_full_tensor = tensor;

      //! Type traits
      typedef NavierStokesParameterTraits<GV,RF> Traits;

      //! Constructor
      NavierStokesDefaultParameters(const Dune::ParameterTree& config,
                                    F& f,
                                    B& b,
                                    V& v,
                                    J& j)
        : _rho(config.get<RF>("rho"))
        , _mu(config.get<RF>("mu"))
        , _f(f)
        , _b(b)
        , _v(v)
        , _j(j)
      {}

      NavierStokesDefaultParameters(const RF& mu,
                                    const RF& rho,
                                    F& f,
                                    B& b,
                                    V& v,
                                    J& j)
        : _rho(rho)
        , _mu(mu)
        , _f(f)
        , _b(b)
        , _v(v)
        , _j(j)
      {}


      //! source term
      template<typename EG>
      typename Traits::VelocityRange
      f(const EG& e, const typename Traits::Domain& x) const
      {
        typename F::Traits::RangeType fvalue;
        return evaluateVelocityGridFunction(_f,e,x);
      }

      //! boundary condition type from local intersection coordinate
      template<typename IG>
      typename Traits::BoundaryCondition::Type
      bctype(const IG& is,
             const typename Traits::IntersectionDomain& x) const
      {
        typename B::Traits::RangeType y;
        _b.evaluate(is,x,y);
        return y;
      }

      //! Dynamic viscosity value from local cell coordinate
      template<typename EG>
      typename Traits::RangeField
      mu(const EG& e, const typename Traits::Domain& x) const
      {
        return _mu;
      }

      //! Dynamic viscosity value from local intersection coordinate
      template<typename IG>
      typename Traits::RangeField
      mu(const IG& ig, const typename Traits::IntersectionDomain& x) const
      {
        return _mu;
      }

      //! Density value from local cell coordinate
      template<typename EG>
      typename Traits::RangeField
      rho(const EG& eg, const typename Traits::Domain& x) const
      {
        return _rho;
      }

      //! Density value from local intersection coordinate
      template<typename IG>
      typename Traits::RangeField
      rho(const IG& ig, const typename Traits::IntersectionDomain& x) const
      {
        return _rho;
      }

      //! Dirichlet boundary condition value from local cell coordinate
      template<typename EG>
      typename Traits::VelocityRange
      g(const EG& e, const typename Traits::Domain& x) const
      {
        typename V::Traits::RangeType y;
        _v.evaluate(e,x,y);
        return y;
      }

      //! pressure source term
      template<typename EG>
      typename Traits::RangeField
      g2(const EG& e, const typename Traits::Domain& x) const
      {
        return 0;
      }

#ifdef DOXYGEN

      //! Neumann boundary condition (stress)
      template<typename IG>
      typename Traits::VelocityRange>
      j(const IG& ig,
        const typename Traits::IntersectionDomain& x,
        const typename Traits::Domain& normal) const;

#else // DOXYGEN

      //! Neumann boundary condition (stress) - version for scalar function
      template<typename IG>
      typename std::enable_if<
        J::Traits::dimRange == 1 &&
        (GV::dimension > 1) &&
        AlwaysTrue<IG>::value, // required to force lazy evaluation
        typename Traits::VelocityRange
          >::type
          j(const IG& ig,
            const typename Traits::IntersectionDomain& x,
            typename Traits::Domain normal) const
      {
        typename J::Traits::RangeType r;
        auto e = ig.inside();
        _j.evaluate(e,ig.geometryInInside().global(x),r);
        normal *= r;
        return normal;
      }

      //! Neumann boundary condition (stress) - version for vector-valued function
      template<typename IG>
      typename std::enable_if<
        J::Traits::dimRange == GV::dimension &&
        AlwaysTrue<IG>::value, // required to force lazy evaluation
        typename Traits::VelocityRange
        >::type
        j(const IG& ig,
          const typename Traits::IntersectionDomain& x,
          const typename Traits::Domain& normal) const
      {
        auto e = ig.inside();
        typename J::Traits::RangeType y;
        _j.evaluate(e,ig.geometryInInside().global(x),y);
        return y;
      }

#endif // DOXYGEN

      void setTime(RF time)
      {
        _f.setTime(time);
        _b.setTime(time);
        _v.setTime(time);
        _j.setTime(time);
      }

    private:
      const RF _rho;
      const RF _mu;
      F& _f;
      B& _b;
      V& _v;
      J& _j;
    };


    /**
     * Stokes velocity boundary constraints function
     */
    template<typename PRM>
    class StokesVelocityDirichletConstraints
      : public Dune::PDELab::DirichletConstraintsParameters
    {
    private:
      const PRM & prm_;

    public:

      /** \brief Constructor */
      StokesVelocityDirichletConstraints (const PRM & _prm)
        : prm_(_prm) { }

      /** Predicate identifying Dirichlet boundaries for velocity. */
      template<typename I>
      bool isDirichlet(const I & intersection,
                       const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        StokesBoundaryCondition::Type bctype = prm_.bctype(intersection,coord);
        return (bctype == StokesBoundaryCondition::VelocityDirichlet);
      }

    };

    /**
     * Stokes pressure boundary constraints function
     */
    template<typename PRM>
    class StokesPressureDirichletConstraints
      : public Dune::PDELab::DirichletConstraintsParameters
    {
    private:
      const PRM & prm_;

    public:

      /** \brief Constructor */
      StokesPressureDirichletConstraints (const PRM & _prm)
        : prm_(_prm) { }

      /** Predicate identifying Dirichlet boundaries for velocity. */
      template<typename I>
      bool isDirichlet(const I & intersection,
                       const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord) const
      { return false; }
    };



#ifndef DOXYGEN

    /**
     * Common base class for NavierStokesParameters -> GridFunction adapters.
     */
    template<typename PRM, int rangeDim>
    class NavierStokesFunctionAdapterBase
      : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<
        typename PRM::Traits::GridView,
        typename PRM::Traits::RangeField,
        rangeDim,
        Dune::FieldVector<typename PRM::Traits::RangeField,rangeDim>
        >,
      NavierStokesFunctionAdapterBase<PRM,rangeDim>
      >
    {
    public:
      //! Traits class
      typedef Dune::PDELab::GridFunctionTraits<
      typename PRM::Traits::GridView,
      typename PRM::Traits::RangeField,
      rangeDim,
      Dune::FieldVector<typename PRM::Traits::RangeField,rangeDim>
      > Traits;

      //! Constructor
      NavierStokesFunctionAdapterBase(PRM& prm)
        : _prm(prm)
      {}

      void setTime(const double time)
      {
        _prm.setTime(time);
      }

      const PRM& parameters() const
      {
        return _prm;
      }

      //! Access to underlying grid view
      const typename Traits::GridViewType& getGridView () const
      {
        return _prm.gridView();
      }

    private:
      PRM& _prm;
    };


#endif // DOXYGEN

    /**
     * Adapter that extracts force density Dirichlet boundary
     * conditions from parameter class
     *
     * \tparam PRM   Model of NavierStokesDefaultParameters
     */
    template<typename PRM>
    class NavierStokesVelocityFunctionAdapter
      : public NavierStokesFunctionAdapterBase<PRM,PRM::Traits::dimDomain>
    {

      //! Base class
      typedef NavierStokesFunctionAdapterBase<PRM,PRM::Traits::dimDomain> Base;

      using Base::parameters;

    public:

      typedef typename Base::Traits Traits;

      //! Constructor
      NavierStokesVelocityFunctionAdapter(PRM& prm)
        : Base(prm)
      {}

      //! Evaluate dirichlet function
      void evaluate (const typename Traits::ElementType& e,
                     const typename Traits::DomainType& x,
                     typename Traits::RangeType& y) const
      {
        y = parameters().g(e,x);
      }
    };



#if 0
    /** \brief Factory for a Dirichlet function which can be used
        for interpolation. */
    template < typename PRM >
    class NavierStokesDirichletFunctionAdapterFactory
    {
    public:
      typedef Dune::PDELab::CompositeGridFunction<
      NavierStokesDirichletFunctionAdapter<PRM>,
      NavierStokesPressureDirichletFunctionAdapter<PRM> >
      BoundaryDirichletFunction;

      NavierStokesDirichletFunctionAdapterFactory(PRM & prm)
        : v(prm), p(prm), df(v,p)
      {}

      BoundaryDirichletFunction & dirichletFunction()
      {
        return df;
      }

    private:
      NavierStokesVelocityDirichletFunctionAdapter<PRM> v;
      NavierStokesPressureDirichletFunctionAdapter<PRM> p;
      BoundaryDirichletFunction df;
    };
#endif


  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_STOKESPARAMETER_HH
