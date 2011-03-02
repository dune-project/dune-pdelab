#ifndef DUNE_UDG_RESIDUALENGINE_HH
#define DUNE_UDG_RESIDUALENGINE_HH

namespace Dune{
  namespace UDG{

    /**
       \brief The local assembler engine for UDG sub triangulations which
       assembles the residual vector

       \tparam LA The local udg assembler

    */
    template<typename OSLA>
    class OneStepLocalResidualAssemblerEngine
    {
    public:
      //! The type of the wrapping local assembler
      typedef OSLA OneStepLocalAssembler;

      typedef OSLA::LocalAssemblerDT0 LocalAssemblerDT0;
      typedef OSLA::LocalAssemblerDT1 LocalAssemblerDT1;

      //! The type of the residual vector
      typedef typename OSLA::Residual Residual;
      typedef typename Residual::ElementType ResidualElement;

      //! The type of the solution vector
      typedef typename OSLA::Solution Solution;
      typedef typename Solution::ElementType SolutionElement;

      //! The type for real numbers
      typedef typename OSLA::Real Real;

      //! The sub triangulation type
      typedef typename LA::SubTriangulation SubTriangulation;

      /**
         \brief Constructor 

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      OneStepLocalResidualAssemblerEngine(const LocalAssembler & local_assembler_)
        : la(local_assembler_), lae0(invalid_lae0), lae1(invalid_lae1), 
          invalid_residual(static_cast<Residual*>(0)), 
          invalid_solutions(static_cast<Solution*>(0)),
          const_residual(invalid_residual), solutions(invalid_solutions), stage(0)
      {}

      //! Query methods for the global grid assembler
      //! @{
      bool requireSkeleton() const 
      { return lae0->requireSkeleton() || lae1->requireSkeleton(); }
      bool requireSkeletonTwoSided() const
      { return lae0->requireSkeletonTwoSided() || lae1->requireSkeletonTwoSided(); }
      bool requireUVVolume() const
      { return lae0->doPatternVolume() || lae1->doPatternVolume(); }
      bool requireVVolume() const
      { return lae0->requireVVolume() || lae1->requireVVolume(); }
      bool requireUVSkeleton() const
      { return lae0->requireUVSkeleton() || lae1->requireUVSkeleton(); }
      bool requireVSkeleton() const
      { return lae0->requireVSkeleton() || lae1->requireVSkeleton(); }
      bool requireUVBoundary() const
      { return lae0->requireUVBoundary() || lae1->requireUVBoundary(); }
      bool requireVBoundary() const
      { return lae0->requireVBoundary() || lae1->requireVBoundary(); }
      bool requireUVEnrichedCoupling() const
      { return lae0->requireUVEnrichedCoupling() || lae1->requireUVEnrichedCoupling(); }
      bool requireVEnrichedCoupling() const
      { return lae0->requireVEnrichedCoupling() || lae1->requireVEnrichedCoupling(); }
      bool requireUVVolumePostSkeleton() const
      { return lae0->requireUVVolumePostSkeleton() || lae1->requireUVVolumePostSkeleton();}
      bool requireVVolumePostSkeleton() const
      { return lae0->requireVVolumePostSkeleton() || lae1->requireVVolumePostSkeleton(); }
      //! @}


      //! Public access to the wrapping local assembler
      const LocalAssembler & localAssembler(){ return la; }

      //! Set current const residual vector. Should be called prior to
      //! assembling.
      void setConstResidual(Residual & const_residual_){
        const_residual = &const_residual_;
      }

      //! Set current solution vector. Should be called prior to
      //! assembling.
      void setSolutions(const Solutions & solutions_){
        solutions = &solutions_;
      }

      //! Set current stage number
      void setStage(const int stage_){
        stage = stage_;
      }
    
      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG, typename LFSU>
      void onBindLFSU(const EG & eg, const LFSU & lfsu){
        lae0->onBindLFSU(eg,lfsu);
        lae1->onBindLFSU(eg,lfsu);
      }

      template<typename EG, typename LFSV>
      void onBindLFSV(const EG & eg, const LFSV & lfsv){
        lae0->onBindLFSV(eg,lfsv);
        lae1->onBindLFSV(eg,lfsv);
      }

      template<typename IG, typename LFSU>
      void onBindLFSUInside(const IG & ig, const LFSU & lfsu){
        lae0->onBindLFSUInside(eg,lfsu);
        lae1->onBindLFSUInside(eg,lfsu);
      }

      template<typename IG, typename LFSU>
      void onBindLFSUOutside(const IG & ig, const LFSU & lfsun){
        lae0->onBindLFSUOutside(eg,lfsun);
        lae1->onBindLFSUOutside(eg,lfsun);
      }

      template<typename IG, typename LFSV>
      void onBindLFSVInside(const IG & ig, const LFSV & lfsv){
        lae0->onBindLFSVInside(eg,lfsv);
        lae1->onBindLFSVInside(eg,lfsv);
      }

      template<typename IG, typename LFSV>
      void onBindLFSVOutside(const IG & ig, const LFSV & lfsvn){
        lae0->onBindLFSVOutside(eg,lfsvn);
        lae1->onBindLFSVOutside(eg,lfsvn);
      }

      template<typename LFSU>
      void onBindLFSUCoupling(const LFSU & lfsu){}
      template<typename LFSU>
      void onUnbindLFSUCoupling(const LFSU & lfsu){}
      template<typename LFSV>
      void onBindLFSVCoupling(const LFSV & lfsv){}
      template<typename LFSV>
      void onUnbindLFSVCoupling(const LFSV & lfsv){}

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded 
      //! @{
      template<typename EG, typename LFSU>
      void onUnbindLFSU(const EG & eg, const LFSU & lfsu){
        lae0->onUnbindLFSU(eg,lfsu);
        lae1->onUnbindLFSU(eg,lfsu);
      }

      template<typename EG, typename LFSV>
      void onUnbindLFSV(const EG & eg, const LFSV & lfsv){
        lae0->onUnbindLFSV(eg,lfsv);
        lae1->onUnbindLFSV(eg,lfsv);
      }

      template<typename IG, typename LFSU>
      void onUnbindLFSUInside(const IG & ig, const LFSU & lfsu){
        lae0->onUnbindLFSUInside(ig,lfsu);
        lae1->onUnbindLFSUInside(ig,lfsu);
      }

      template<typename IG, typename LFSU>
      void onUnbindLFSUOutside(const IG & ig, const LFSU & lfsun){
        lae0->onUnbindLFSUOutside(ig,lfsun);
        lae1->onUnbindLFSUOutside(ig,lfsun);
      }

      template<typename IG, typename LFSV>
      void onUnbindLFSVInside(const IG & ig, const LFSV & lfsv){
        lae0->onUnbindLFSVInside(ig,lfsv);
        lae1->onUnbindLFSVInside(ig,lfsv);
      }

      template<typename IG, typename LFSV>
      void onUnbindLFSVOutside(const IG & ig, const LFSV & lfsvn){
        lae0->onUnbindLFSVOutside(ig,lfsvn);
        lae1->onUnbindLFSVOutside(ig,lfsvn);
      }

      //! @}

      //! Methods for loading of the local function's
      //! coefficients. These methods are blocked. The loading of the
      //! coefficients is done in each assemble call.
      //!@{
      template<typename LFSU>
      void loadCoefficientsLFSUInside(const LFSU & lfsu_s){}
      template<typename LFSU>
      void loadCoefficientsLFSUOutside(const LFSU & lfsu_n){}
      template<typename LFSU>
      void loadCoefficientsLFSUCoupling(const LFSU & lfsu_c){}
      //! @}


      //! Notifier functions, called immediately before and after assembling
      //! @{
      void preAssembly()
      {
        // Initialize constant part of residual
        *const_residual = 0.;

        // Set residual vector of subordinate local assemblers
        lae0->setResidual(*const_residual);
        lae1->setResidual(*const_residual);

        // Extract the coefficients of the time step scheme
        a.resize(stage);
        b.resize(stage);
        d.resize(stage);
        for (size_t i=0; i<stage; ++i){ 
          a[i] = la.method.a(stage,i);
          b[i] = la.method.b(stage,i);
          d[i] = la.method.d(i);
          do0[i] = ( std::abs(a[i]) < 1E-6 );
          do1[i] = ( std::abs(b[i]) < 1E-6 );
        }

        // prepare local operators for stage
        lae0->preStage(la.time+la.method.d(stage)*la.dt,stage);
        lae1->preStage(la.time+la.method.d(stage)*la.dt,stage);

      }
      void postAssembly()
      { 
        
      }
      //! @}

      //! Assembling methods
      //! @{
      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        for (unsigned s=0; s<stage; ++s){
          // Reset the time in the local assembler
          la->setTime(la.time+d[s]*la.dt);

          // Set assembling weights
          lae0->localAssembler().setWeight(b[s]*la.dt);
          lae1->localAssembler().setWeight(a[s]);

          // Compute residual for current solution
          lae0->setSolution(solutions[s]);
          lae0->loadCoefficientsLFSUInside(lfsu);
          lae0->assembleUVVolume(eg,lfsu,lfsv);

          lae1->setSolution(solutions[s]);
          lae1->loadCoefficientsLFSUInside(lfsu);
          lae1->assembleUVVolume(eg,lfsu,lfsv);
        }
      }
      
      template<typename EG, typename LFSV>
      void assembleVVolume(const EG & eg, const LFSV & lfsv)
      {
        for (unsigned s=0; s<stage; ++s){
          // Reset the time in the local assembler
          la->setTime(la.time+d[s]*la.dt);

          // Set assembling weights
          lae0->localAssembler().setWeight(b[s]*la.dt);
          lae1->localAssembler().setWeight(a[s]);

          // Compute residual for current solution
          lae0->assembleVVolume(eg,lfsv);
          lae1->assembleVVolume(eg,lfsv);
        }
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void assembleUVSkeleton(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                              const LFSU_N & lfsu_n, const LFSV_N & lfsv_n)
      {
        for (unsigned s=0; s<stage; ++s){
          // Reset the time in the local assembler
          la->setTime(la.time+d[s]*la.dt);

          // Set assembling weights
          lae0->localAssembler().setWeight(b[s]*la.dt);
          lae1->localAssembler().setWeight(a[s]);

          // Compute residual for current solution
          lae0->setSolution(solutions[s]);
          lae0->loadCoefficientsLFSUInside(lfsu_s);
          lae0->loadCoefficientsLFSUOutside(lfsu_n);
          lae0->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);

          lae1->setSolution(solutions[s]);
          lae1->loadCoefficientsLFSUInside(lfsu_s);
          lae1->loadCoefficientsLFSUOutside(lfsu_n);
          lae1->assembleUVSkelton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
        }
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void assembleVSkeleton(const IG & ig, const LFSV_S & lfsv_s, const LFSV_N & lfsv_n)
      {
        for (unsigned s=0; s<stage; ++s){
          // Reset the time in the local assembler
          la->setTime(la.time+d[s]*la.dt);

          // Set assembling weights
          lae0->localAssembler().setWeight(b[s]*la.dt);
          lae1->localAssembler().setWeight(a[s]);

          // Compute residual for current solution
          lae0->assembleVSkeleton(ig,lfsv_s,lfsv_n);
          lae1->assembleVSkelton(ig,lfsv_s,lfsv_n);
        }
      }

      template<typename IG, typename LFSU_S, typename LFSV_S>
      void assembleUVBoundary(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s)
      {
        for (unsigned s=0; s<stage; ++s){
          // Reset the time in the local assembler
          la->setTime(la.time+d[s]*la.dt);

          // Set assembling weights
          lae0->localAssembler().setWeight(b[s]*la.dt);
          lae1->localAssembler().setWeight(a[s]);

          // Compute residual for current solution
          lae0->setSolution(solutions[s]);
          lae0->loadCoefficientsLFSUInside(lfsu_s);
          lae0->assembleUVBoundary(ig,lfsu_s,lfsv_s);

          lae1->setSolution(solutions[s]);
          lae1->loadCoefficientsLFSUInside(lfsu_s);
          lae1->assembleUVSkelton(ig,lfsu_s,lfsv_s);
        }
      }

      template<typename IG, typename LFSV_S>
      void assembleVBoundary(const IG & ig, const LFSV_S & lfsv_s)
      {
        for (unsigned s=0; s<stage; ++s){
          // Reset the time in the local assembler
          la->setTime(la.time+d[s]*la.dt);

          // Set assembling weights
          lae0->localAssembler().setWeight(b[s]*la.dt);
          lae1->localAssembler().setWeight(a[s]);

          // Compute residual for current solution
          lae0->assembleVBoundary(ig,lfsv_s);
          lae1->assembleVBoundary(ig,lfsv_s);
        }
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N, 
               typename LFSU_C, typename LFSV_C>
      static void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                                             const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                                             const LFSU_C & lfsu_c, const LFSV_C & lfsv_c)
      {
        for (unsigned s=0; s<stage; ++s){
          // Reset the time in the local assembler
          la->setTime(la.time+d[s]*la.dt);

          // Set assembling weights
          lae0->localAssembler().setWeight(b[s]*la.dt);
          lae1->localAssembler().setWeight(a[s]);

          // Compute residual for current solution
          lae0->setSolution(solutions[s]);
          lae0->loadCoefficientsLFSUInside(lfsu_s);
          lae0->loadCoefficientsLFSUOutside(lfsu_n);
          lae0->loadCoefficientsLFSUCoupling(lfsu_c);
          lae0->assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_c,lfsv_c);

          lae1->setSolution(solutions[s]);
          lae1->loadCoefficientsLFSUInside(lfsu_s);
          lae1->loadCoefficientsLFSUOutside(lfsu_n);
          lae1->loadCoefficientsLFSUCoupling(lfsu_c);
          lae1->assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_c,lfsv_c);
        }
      }

      template<typename IG, typename LFSV_S, typename LFSV_N, typename LFSV_C>
      static void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSV_S & lfsv_s,
                                            const LFSV_N & lfsv_n,
                                            const LFSV_C & lfsv_c) 
      {
        for (unsigned s=0; s<stage; ++s){
          // Reset the time in the local assembler
          la->setTime(la.time+d[s]*la.dt);

          // Set assembling weights
          lae0->localAssembler().setWeight(b[s]*la.dt);
          lae1->localAssembler().setWeight(a[s]);

          // Compute residual for current solution
          lae0->assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_c);
          lae1->assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_c);
        }
      }

      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        for (unsigned s=0; s<stage; ++s){
          // Reset the time in the local assembler
          la->setTime(la.time+d[s]*la.dt);

          // Set assembling weights
          lae0->localAssembler().setWeight(b[s]*la.dt);
          lae1->localAssembler().setWeight(a[s]);

          // Compute residual for current solution
          lae0->setSolution(solutions[s]);
          lae0->loadCoefficientsLFSUInside(lfsu);
          lae0->assembleUVVolumePostSkeleton(eg,lfsu,lfsv);

          lae1->setSolution(solutions[s]);
          lae1->loadCoefficientsLFSUInside(lfsu);
          lae1->assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
        }
      }

      template<typename EG, typename LFSV>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv)
      {
        for (unsigned s=0; s<stage; ++s){
          // Reset the time in the local assembler
          la->setTime(la.time+d[s]*la.dt);

          // Set assembling weights
          lae0->localAssembler().setWeight(b[s]*la.dt);
          lae1->localAssembler().setWeight(a[s]);

          // Compute residual for current solution
          lae0->assembleVVolumePostSkeleton(eg,lfsv);
          lae1->assembleVVolumePostSkeleton(eg,lfsv);
        }
      }
      //! @}

    private:
      //! Reference to the wrapping local assembler object which
      //! constructed this engine
      const LocalAssembler & la;

      typedef typename LocalAssemblerDT0::LocalPatternAssemblerEngine PatternEngineDT0;
      typedef typename LocalAssemblerDT1::LocalPatternAssemblerEngine PatternEngineDT1;

      PatternEngineDT0 * const invalid_lae0;
      PatternEngineDT0 * lae0;
      PatternEngineDT1 * const invalid_lae1;;
      PatternEngineDT1 * lae1;

      //! Default value indicating an invalid residual pointer
      Residual * const invalid_residual;

      //! Default value indicating an invalid solution pointer
      Solutions * const invalid_solutions;

      //! Pointer to the current constant part residual vector in
      //! which to assemble
      Residual * const_residual;

      //! Pointer to the current residual vector in which to assemble
      const Solutions * solutions;

      //! The current stage number
      int stage;

      //! Coefficients of time stepping scheme
      std::vector<Real> a;
      std::vector<Real> b;
      std::vector<Real> d;
      std::vector<bool> do0;
      std::vector<bool> do1;


      //! The local vectors and matrices as required for assembling
      //! @{
      typedef Dune::PDELab::TrialSpaceTag LocalTrialSpaceTag;
      typedef Dune::PDELab::TestSpaceTag LocalTestSpaceTag;
      Dune::PDELab::LocalVector<SolutionElement, LocalTrialSpaceTag> xl;
      Dune::PDELab::LocalVector<ResidualElement, LocalTestSpaceTag> rl0;
      Dune::PDELab::LocalVector<ResidualElement, LocalTestSpaceTag> rl1;
      Dune::PDELab::LocalVector<SolutionElement, LocalTrialSpaceTag> xn;
      Dune::PDELab::LocalVector<ResidualElement, LocalTestSpaceTag> rn;
      //! @}

    
    }; // End of class OneStepLocalResidualAssemblerEngine

  };
};
#endif
