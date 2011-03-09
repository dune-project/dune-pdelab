#include <gridoperatorutilities.hh>
#include <assemblerutilties.hh>

namespace Dune {
  namespace PDELab {

    //! \defgroup GridOperator Grid Operator
    //! \ingroup PDELab
    //! \{

    /** \brief The global assembler which performs the traversing of
        the integration parts.

        The global assembler does only provide the local function
        spaces and the integration parts. It does not perform the
        actual integration or modify any of the assembling objects
        like the pattern, residual, and jacobian matrix.

     */
    class AssemblerInterface{
    public:
      template<class LocalAssemblerEngine>
      void assemble(LocalAssemblerEngine & local_assembler_engine);
    };


    /** \brief The local assembler engine which handles the
        integration parts as provided by the global assemblers.

     */
    class LocalAssemblerEngine{
    public:
      //! The type of the local assembler
      typedef LocalAssemblerInterface LocalAssembler;
      
      //! Access to the superior local assembler object
      const LocalAssembler & localAssembler();

      /** @name Query methods 

          Query methods indicating which assembling methods need to be
          called by the global assembler.
      
      @{
      */
      bool requireSkeleton() const;
      bool requireSkeletonTwoSided() const;
      bool requireUVVolume() const;
      bool requireVVolume() const;
      bool requireUVSkeleton() const;
      bool requireVSkeleton() const;
      bool requireUVBoundary() const;
      bool requireVBoundary() const;
      bool requireUVEnrichedCoupling() const;
      bool requireVEnrichedCoupling() const;
      bool requireUVVolumePostSkeleton() const;
      bool requireVVolumePostSkeleton() const;
      //! @}

      /**
         @name Assembling methods

         All local function spaces as provided in these methods are
         already bound to the grid cell corresponding to the entity
         part \a eg or intersection part \a ig .

         @{
       */

      /** Assembling method which is called for a given grid cell. It
      is called before the local function spaces are bound to the cell
      and the coefficients for the local trial function space are
      extracted.
      
      \return Indicate whether assembling of this cell may be aborted
      after the call of this method. This may avoid unneccessary costs
      due to binding of the local function spaces etc. 
      */
      template<typename EG>
      bool assembleCell(const EG & eg);

      /** Assembling for a codim 0 entity part for trial and test
      local function spaces. 
      */
      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv);

      /** Assembling for a codim 0 entity part for test local function
      spaces. 
      */
      template<typename EG, typename LFSV>
      void assembleVVolume(const EG & eg, const LFSV & lfsv);

      /** Assembling for an interior codim 1 entity part for trial and
      test local function spaces. 
      */
      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void assembleUVSkeleton(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                              const LFSU_N & lfsu_n, const LFSV_N & lfsv_n);

      /** Assembling for an interior codim 1 entity part for test
      local function spaces. 
      */
      template<typename IG, typename LFSV_S, typename LFSV_N>
      void assembleVSkeleton(const IG & ig, const LFSV_S & lfsv_s, const LFSV_N & lfsv_n);

      /** Assembling for a boundary codim 1 entity part for trial and
      test local function spaces. 
      */
      template<typename IG, typename LFSU_S, typename LFSV_S>
      void assembleUVBoundary(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s);

      /** Assembling for a boundary codim 1 entity part for test local
      function spaces. 
      */
      template<typename IG, typename LFSV_S>
      void assembleVBoundary(const IG & ig, const LFSV_S & lfsv_s);

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N, 
               typename LFSU_C, typename LFSV_C>
      void assembleUVEnrichedCoupling(const IG & ig,
                                      const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                                      const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                                      const LFSU_C & lfsu_c, const LFSV_C & lfsv_c);

      template<typename IG, typename LFSV_S, typename LFSV_N, typename LFSV_C>
      void assembleVEnrichedCoupling(const IG & ig,
                                     const LFSV_S & lfsv_s,
                                     const LFSV_N & lfsv_n,
                                     const LFSV_C & lfsv_c);

      /** Assembling for a codim 0 entity part for trial and test
      local function spaces which is called after the intersection
      parts of the current cell have been handled.
      */
      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSU & lfsu, const LFSV & lfsv);

      /** Assembling for a codim 0 entity part for test local function
      spaces which is called after the intersection parts of the
      current cell have been handled.
      */
      template<typename EG, typename LFSV>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv);

      //! Called directly before assembling
      void preAssembly();

      //! Called last thing after assembling
      void postAssembly();

      /**
         @}
       */

      /**
         @name Notifications

         Notification methods called by the global assembler when
         binding and unbinding the local function spaces.

         @{
       */
      void onBindLFSUV(const EG & eg, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s);
      void onBindLFSV(const EG & eg, const LFSV_S & lfsv_s);
      void onBindLFSUVInside(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s);
      void onBindLFSVInside(const IG & ig, const LFSV_S & lfsv_s);
      void onBindLFSUVOutside(const IG & ig, const LFSU_N & lfsu_n, const LFSV_N & lfsv_n);
      void onBindLFSVOutside(const IG & ig, const LFSV_N & lfsv_n);
      void onBindLFSUVCoupling(const IG & ig, 
                               const LFSU_Coupling & lfsu_coupling, const LFSV_Coupling & lfsv_coupling);
      void onBindLFSVCoupling(const IG & ig, const LFSV_Coupling & lfsv_coupling);

      void onUnbindLFSUV(const EG & eg, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s);
      void onUnbindLFSV(const EG & eg, const LFSV_S & lfsv_s);
      void onUnbindLFSUVInside(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s);
      void onUnbindLFSVInside(const IG & ig, const LFSV_S & lfsv_s);
      void onUnbindLFSUVOutside(const IG & ig, const LFSU_N & lfsu_n, const LFSV_N & lfsv_n);
      void onUnbindLFSVOutside(const IG & ig, const LFSV_N & lfsv_n);
      void onUnbindLFSUVCoupling(const IG & ig, 
                                 const LFSU_Coupling & lfsu_coupling, const LFSV_Coupling & lfsv_coupling);
      void onUnbindLFSVCoupling(const IG & ig, const LFSV_Coupling & lfsv_coupling);

      /** @} */

      /**
         @name Loading trial space coefficients

         Tells the engine to load the local coefficients for the given
         local function space.

         @{
       */
      void loadCoefficientsLFSUInside(const LFSU_S & lfsu_s);
      void loadCoefficientsLFSUOutside(const LFSU_N & lfsu_n);
      void loadCoefficientsLFSUCoupling(const LFSU_Coupling & lfsu_coupling);
      /** @} */


      /**
         @name Assign the assembler target objects

         These methods assign the objects into which the assembler
         should assemble i.e. the solution vector from which the local
         coefficients are to be extracted.

         @{
       */
      void setSolution(const X& x);
      void setPattern(const P& p);
      void setJacobian(const J & j);
      void setResidual(const R& r);
      /** @} */

    };

    /** \brief The local assembler which provides the engines that
        drive the global assembler.

        The local assembler provides engines for the standard
        operations of the grid operator. This includes setting up the
        pattern, computing the residual and the jacobian matrix. 

        It also provides a standard interface which may be used by
        implementations of time stepping methods.

     */
    template<typename B, typename CU, typename CV>
    class LocalAssemblerInterface : public LocalAssemblerBase<B,CU,CV>
    {
    public:

      template<class TT>
      void setTime(TT time);

      template<class RF>
      void setWeight(RF weight);

      template<typename TT>
      void preStep (TT time, TT dt, std::size_t stages);

      void postStep ();

      template<typename TT>
      void preStage (TT time, std::size_t stage);

      void postStage ();

      template<typename TT>
      TT suggestTimestep (TT dt) const;

      LocalPatternAssemblerEngine & localPatternAssemblerEngine(P & p);
      LocalResidualAssemblerEngine & localResidualAssemblerEngine(R & r, const X & x);
      LocalJacobianAssemblerEngine & localJacobianAssemblerEngine(A & a, const X & x);
      LocalResidualJacobianAssemblerEngine & localResidualJacobianAssemblerEngine(R & r, A & a, const X & x);

      

      class LocalPatternAssemblerEngine : public LocalAssemblerEngine {};
      class LocalResidualAssemblerEngine : public LocalAssemblerEngine {};
      class LocalJacobianAssemblerEngine : public LocalAssemblerEngine {};
      class LocalResidualJacobianAssemblerEngine : public LocalAssemblerEngine {};
    };

    template<typename GFSU, typename GFSV,
             typename MB, typename DF, typename RF, typename JF>
    class GridOperatorInterface{
    public:

      typedef GridOperatorTraits
      <GFSU,GFSV,MB,DF,RF,JF,CU,CV,AssemblerInterface,LocalAssemblerInterface> Traits;

      template<typename P>
      void fill_pattern (P& globalpattern) const;

      template<typename X, typename R>
      void residual (const X& x, R& r) const;

      template<typename X, typename A>
      void jacobian (const X& x, A& a) const;

      Assembler & assembler();

      LocalAssemblerInterface & localAssembler();

      const GFSU& trialGridFunctionSpace() const;

      const GFSV& testGridFunctionSpace() const;

      typename GFSU::Traits::SizeType globalSizeU () const;

      typename GFSV::Traits::SizeType globalSizeV () const;

      //! Interpolate xnew from f, taking unconstrained values from xold.
      /**
       * \note the exact type of F will depend on the GridOperator and
       *       may be a more complicated object than a simple GridFunction
       *       for scenarios like MultiDomain or grid-glue.
       */
      template<typename F>
      void interpolate(const typename Traits::Domain& xold,
                       const F& f,
                       typename Traits::Domain& xnew);

      //! Set up the passed-in tuple of GridOperators to cooperate, e.g.
      //! for a time-stepping method. The caller guarantees that the
      //! GridOperators will always be invoked in the order that they
      //! appear in the tuple.
      template<typename GridOperatorTuple>
      static void setupGridOperators(GridOperatorTuple& tuple);

    };
    //! \}
  };
};

/**
\page GridOperatorDocPage Grid Operator
@ingroup GridOperator

\section GridOperatorDocIntroduction Introduction

In the PDELab concept, the continuous PDE problem is reduced to an
algebraic problem: 

Find \f$ \mathbf{u}\in\mathbf{U} \f$ such that \f$
\mathcal{R}(\mathbf{u}) = \mathbf{0} \f$. 

For instationary problems a corresponding algebraic problem is setup
for each time step or stage of a time step.

The grid operator object represents the operator mapping \f$
\mathcal{R} : \mathbf{U} \to \mathbf{V} \f$ . It is evaluated via the
Dune::PDELab::GridOperatorInterface::residual() member method and its
derivatives are aquired with the
Dune::PDELab::GridOperatorInterface::jacobian() method. 

Evaluating the grid operator and its jacobian matrix entails
integrations over the computational domain during which the
corresponding algebraic objects are assembled incrementally. The
assembling is provided by two objects corresponding to the
Dune::PDELab::AssemblerInterface and the
Dune::PDELab::LocalAssemblerInterface . The former, the global
assembler, provides the geometric objects, representing parts of the
computational domain to be integrated. The latter, the local
assembler, calls the local operator in an appropriate way to compute
the local integrals and afterwards writes the results into the algebra
objects. 

The separation of these two tasks into different objects has
significant advantages to a monolithic approach. It allows different
implementations of the assembler interface to be used interchangeably
in common implementations of time stepping schemes. The latter may
link in between the global and the local assembler and thus apply the
necessary modifications to the local integrations (like multiplying
with the time step with or Runge-Kutta coefficients). 

Furthermore, the separation provides a junction for caching
objects. Such objects would provide the interface of the local
assembler, while actually wrapping the true local assembler object.

*/
