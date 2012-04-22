#ifndef DUNE_PDELAB_LOCALOPERATOR_STOKESPARAMETER_HH
#define DUNE_PDELAB_LOCALOPERATOR_STOKESPARAMETER_HH

#include <dune/pdelab/constraints/constraintsparameters.hh>
#include <dune/common/parametertree.hh>
#include <dune/pdelab/common/function.hh>

namespace Dune {
    namespace PDELab {

        /**
           These are the boundary condition types as to be returned by
           the employed boundary type function. 

           Possible types:

           <ul>

           <li>\a DoNothing : Do not evaluate boundary integrals.

           <li>\a VelocityDirichlet : Dirichlet conditions for velocity.

           <li>\a PressureDirichlet : Natural Neumann conditions for the
           impulse flux. These are equivalent to a fixed pressure
           condition \b if \f$ \forall i : n \cdot \nabla v_i = 0 \f$.

           </ul>
         */
        struct StokesBoundaryCondition {
            enum Type {
                DoNothing = 0,
                VelocityDirichlet = 1,
                PressureDirichlet = 2,
                SlipVelocity = 3
            };
        };

      /** \brief Traits class for the parameter class of a Navier
       * Stokes local operator  */
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

        //! grid types
        typedef typename GV::Traits::template Codim<0>::Entity Element;
        typedef typename GV::Intersection Intersection;
      };


      /** \brief Compile-time switch allowing the evaluation of a
          vector valued grid function. */
      template< typename GF, typename Switch = void >
      struct evaluateVelocityGridFunction{
        template< typename Entity, typename Domain, typename Range >
        static void apply(const GF & gf, const Entity & e, const Domain & x, Range & y)
        {
          gf.evaluate(e,x,y);
        }
      };

      /** \brief Compile-time switch allowing the evaluation of a
          power grid function. */
      template< typename GF>
      struct evaluateVelocityGridFunction
      <GF,typename Dune::enable_if<
            Dune::is_same<typename GF::ImplementationTag,
                          typename Dune::PDELab::PowerGridFunctionTag
                          >::value 
            >::type>
      {
        template< typename Entity, typename Domain, typename Range >
        static void apply(const GF & gf, const Entity & e, const Domain & x, Range & y){
          typename GF::template Child<0>::Type::Traits::RangeType sy;
          for(int c=0; c<Domain::dimension; ++c){
            gf.child(c).evaluate(e,x,sy);
            y[c] = sy;
          }
        }
      };

      /** \brief Default implementation for the parameter class to be
          used with the Taylor-Hood Navier-Stokes local operator. 

          This is designed to work with the TaylorHoodNavierStokes,
          TaylorHoodNavierStokesJacobian and NavierStokesMass local
          operator classes.

          \tparam BF    The boundary type function returning an element
                        of StokesBoundaryCondition.
          \tparam NF    The Neumann stress flux boundary function.
          \tparam DVF   The Dirichlet velocity function.
          \tparam RF    The range field type of the Navier-Stokes solution.
      */
      template <class GV, class BF, class NF, class DVF, class RF>
      class TaylorHoodNavierStokesDefaultParameters
      {
      public:

        //! Type traits
        typedef NavierStokesParameterTraits<GV,RF> Traits;

        //! Constructor
        TaylorHoodNavierStokesDefaultParameters
        (Dune::ParameterTree config, const BF & _bf, const NF & _nf, const DVF & _dvf):
          rho_(config.get<double>("rho")), 
          mu_(config.get<double>("mu")), 
          bf_(_bf), nf_(_nf), dvf_(_dvf)
        {}

        /** \brief Density evaluated on a codim 1 geometry. */
        template<typename IG>
        RF rho(const IG & ig, const typename Traits::IntersectionDomain & x) const
        { return rho_; }

        /** \brief Density evaluated on a codim 0 geometry. */
        template<typename EG>
        RF rho(const EG & eg, const typename Traits::Domain & x) const
        { return rho_; }

        /** \brief Viscosity evaluated on a codim 1 geometry. */
        template<typename IG>
        RF mu(const IG & ig, const typename Traits::IntersectionDomain & x) const
        { return mu_; }

        /** \brief Viscosity evaluated on a codim 0 geometry. */
        template<typename EG>
        RF mu(const EG & eg, const typename Traits::Domain & x) const
        { return mu_; }

        /** \brief General source term representing a source of
            momentum. */
        template<typename EG>
        typename Traits::VelocityRange
        source(const EG & eg, const typename Traits::Domain & x) const
        { return typename Traits::VelocityRange(0); }

        /** \brief Boundary condition type */
        template<typename IG>
        typename Traits::BoundaryCondition::Type
        bcType(const IG & ig, const typename Traits::IntersectionDomain & x) const
        {
          typename Traits::BoundaryCondition::Type y;
          bf_.evaluate(ig,x,y);
          return y;
        }
        
        /** \brief Dirichlet velocity */
        template<typename EG>
        typename Traits::VelocityRange
        velocityDirichlet(const EG & eg, const typename Traits::Domain & x)
        {
          typename Traits::VelocityRange y;
          evaluateVelocityGridFunction<typename DVF::template Child<0>::Type>::
            apply(dvf_.template child<0>(),eg,x,y);
          return y;
        }

        /** \brief Dirichlet pressure */
        template<typename EG>
        typename Traits::PressureRange
        pressureDirichlet(const EG & eg, const typename Traits::Domain & x)
        {
          typename Traits::PressureRange y(0);
          return y;
        }
        
        /** \brief Neumann momentum flux */
        template<typename IG>
        typename Traits::VelocityRange
        stress
        (const IG & ig, 
         const typename Traits::IntersectionDomain & x, 
         typename Traits::Domain normal) const
        {
          typename NF::Traits::RangeType r;
          nf_.evaluate(*(ig.inside()),ig.geometryInInside().global(x),r);
          normal *= r;
          return normal;
        }

      private:
        const RF rho_;
        const RF mu_;
        const BF & bf_;
        const NF & nf_;
        const DVF & dvf_;
      };


      /** \brief Stokes velocity boundary constraints function */
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
        bool isDirichlet(
                         const I & intersection,
                         const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                         ) const
        {
          StokesBoundaryCondition::Type bctype =prm_.bcType(intersection,coord);
          return (bctype == StokesBoundaryCondition::VelocityDirichlet);
        }
      };


      /** \brief Stokes pressure boundary constraints function */
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
        bool isDirichlet(
                         const I & intersection,
                         const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                         ) const
        {
          StokesBoundaryCondition::Type bctype =prm_.bcType(intersection,coord);
          return (bctype == StokesBoundaryCondition::PressureDirichlet);
        }
      };
     


      /*! Adapter that extracts Dirichlet boundary conditions from
        parameter class.

        \tparam T               Model of TaylorHoodNavierStokesDefaultParameters
        \tparam rangeDim        Dimension of range of Dirichlet function
      */
      template<typename T, int rangeDim>
      class NavierStokesDirichletFunctionAdapterBase : 
        public Dune::PDELab::GridFunctionBase<
        Dune::PDELab::GridFunctionTraits<
          typename T::Traits::GridView,
          typename T::Traits::RangeField,
          rangeDim,Dune::FieldVector<typename T::Traits::RangeField,rangeDim> 
          >,
            NavierStokesDirichletFunctionAdapterBase<T,rangeDim> 
            >
      {
      public:
        //! Traits class
        typedef Dune::PDELab::GridFunctionTraits<
        typename T::Traits::GridView,
        typename T::Traits::RangeField,
        rangeDim,Dune::FieldVector<typename T::Traits::RangeField,rangeDim> > Traits;

        //! Constructor 
        NavierStokesDirichletFunctionAdapterBase (T& t_) : t(t_) {}

        void setTime(const double time_){
          t.setTime(time_);
        }

        //! Access to underlying grid view
        inline const typename Traits::GridViewType& getGridView () const { return t.gridView(); }
  
      protected:
        T& t;
      };

      /*! Adapter that extracts force density Dirichlet boundary
        conditions from parameter class

        \tparam T Model of TaylorHoodNavierStokesDefaultParameters
      */
      template<typename T>
      class NavierStokesVelocityDirichletFunctionAdapter : 
        public NavierStokesDirichletFunctionAdapterBase<T,T::Traits::dimDomain>
      {
      public:
        //! Base class
        typedef NavierStokesDirichletFunctionAdapterBase<T,T::Traits::dimDomain> Base;
        //! Constructor 
        NavierStokesVelocityDirichletFunctionAdapter ( T& t_) : Base(t_) {}

        //! Evaluate dirichlet function
        inline void evaluate (const typename Base::Traits::ElementType& e, 
                              const typename Base::Traits::DomainType& x, 
                              typename Base::Traits::RangeType& y) const
        { y = Base::t.velocityDirichlet (e,x); }
      };


      /*! Adapter that extracts pressure Dirichlet boundary conditions
        from parameter class

        \tparam T Model of TaylorHoodNavierStokesDefaultParameters
      */
      template<typename T>
      class NavierStokesPressureDirichletFunctionAdapter : public NavierStokesDirichletFunctionAdapterBase<T,1>
      {
      public:
        //! Base class
        typedef NavierStokesDirichletFunctionAdapterBase<T,1> Base;
        //! Constructor 
        NavierStokesPressureDirichletFunctionAdapter ( T& t_) : Base(t_) {}

        //! Evaluate dirichlet function
        inline void evaluate (const typename Base::Traits::ElementType& e, 
                              const typename Base::Traits::DomainType& x, 
                              typename Base::Traits::RangeType& y) const
        { 
          y = Base::t.pressureDirichlet(e,x);
        }
      };
   
      /** \brief Factory for a Dirichlet function which can be used
          for interpolation. */
      template < typename PRM >
      class NavierStokesDirichletFunctionAdapterFactory
      {
      public:
        typedef Dune::PDELab::CompositeGridFunction
        <NavierStokesVelocityDirichletFunctionAdapter<PRM>,
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



    }
}

#endif
