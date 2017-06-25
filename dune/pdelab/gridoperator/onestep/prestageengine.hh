#ifndef DUNE_PDELAB_GRIDOPERATOR_ONESTEP_PRESTAGEENGINE_HH
#define DUNE_PDELAB_GRIDOPERATOR_ONESTEP_PRESTAGEENGINE_HH

#include <dune/pdelab/gridoperator/onestep/enginebase.hh>
#include <cmath>
#include <vector>

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for one step methods which
       assembles the constant part of the residual vector

       \tparam LA The local one step assembler

    */
    template<typename OSLA>
    class OneStepLocalPreStageAssemblerEngine
      : public OneStepLocalAssemblerEngineBase<OSLA,
                                               typename OSLA::LocalAssemblerDT0::LocalResidualAssemblerEngine,
                                               typename OSLA::LocalAssemblerDT1::LocalResidualAssemblerEngine
                                               >
    {

      typedef OneStepLocalAssemblerEngineBase<OSLA,
                                              typename OSLA::LocalAssemblerDT0::LocalResidualAssemblerEngine,
                                              typename OSLA::LocalAssemblerDT1::LocalResidualAssemblerEngine
                                              > BaseT;

      using BaseT::la;
      using BaseT::lae0;
      using BaseT::lae1;
      using BaseT::implicit;
      using BaseT::setLocalAssemblerEngineDT0;
      using BaseT::setLocalAssemblerEngineDT1;

    public:
      //! The type of the wrapping local assembler
      typedef OSLA LocalAssembler;

      //! Types of the subordinate assemblers and engines
      //! @{
      typedef typename OSLA::LocalAssemblerDT0 LocalAssemblerDT0;
      typedef typename OSLA::LocalAssemblerDT1 LocalAssemblerDT1;

      typedef typename LocalAssemblerDT0::LocalResidualAssemblerEngine ResidualEngineDT0;
      typedef typename LocalAssemblerDT1::LocalResidualAssemblerEngine ResidualEngineDT1;
      //! @}

      //! The type of the residual vector
      typedef typename OSLA::Traits::Residual Residual;
      typedef typename Residual::ElementType ResidualElement;

      //! The type of the solution vector
      typedef typename OSLA::Traits::Solution Solution;
      typedef typename Solution::ElementType SolutionElement;

      //! The type for real numbers
      typedef typename OSLA::Real Real;

      //! The type of the solution container
      typedef std::vector<Solution*> Solutions;

      /**
         \brief Constructor

         \param [in] la_ The local assembler object which
         creates this engine
      */
      OneStepLocalPreStageAssemblerEngine(LocalAssembler & la_)
        : BaseT(la_)
        , invalid_residual(nullptr)
        , invalid_solutions(nullptr)
        , const_residual_0(invalid_residual)
        , const_residual_1(invalid_residual)
        , solutions(invalid_solutions)
      {}

      //! Query methods for the global grid assembler
      //! @{
      bool requireSkeleton() const
      { return lae0->requireSkeleton() || lae1->requireSkeleton(); }
      //! @}

      //! Set current solution vector. Must be called before
      //! setConstResidual()! Should be called prior to assembling.
      void setSolutions(const Solutions & solutions_)
      {
        solutions = &solutions_;
      }

      //! Set current const residual vector. Should be called prior to
      //! assembling.
      void setConstResiduals(Residual & const_residual_0_, Residual & const_residual_1_)
      {
        const_residual_0 = &const_residual_0_;
        const_residual_1 = &const_residual_1_;

        // Initialize the engines of the two wrapped local assemblers
        assert(solutions != invalid_solutions);
        setLocalAssemblerEngineDT0(la.la0.localResidualAssemblerEngine(*const_residual_0,*((*solutions)[0])));
        setLocalAssemblerEngineDT1(la.la1.localResidualAssemblerEngine(*const_residual_1,*((*solutions)[0])));
      }

      //! Set current const residual vector. Should be called prior to
      //! assembling.
      void setConstResidual(Residual & const_residual_)
      {
        const_residual_0 = &const_residual_;
        const_residual_1 = &const_residual_;

        // Initialize the engines of the two wrapped local assemblers
        assert(solutions != invalid_solutions);
        setLocalAssemblerEngineDT0(la.la0.localResidualAssemblerEngine(*const_residual_0,*((*solutions)[0])));
        setLocalAssemblerEngineDT1(la.la1.localResidualAssemblerEngine(*const_residual_1,*((*solutions)[0])));
      }

      //! Methods for binding the local function space.
      //! These methods are empty. The binding of the local function space
      //! is done after setting the solution in the assembleUVVolume(),
      //! assembleUVSkeleton(), assembleUVBoundary(), assembleUVProcessor(),
      //! assembleUVEnrichedCoupling() and assembleUVVolumePostSkeleton() calls.
      //! @{
      template<typename EG, typename LFSU, typename LFSV>
      void onBindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
      {}
      template<typename IG, typename LFSU_S, typename LFSV_S>
      void onBindLFSUVInside(const IG & ig,
                             const LFSU_S & lfsu_s, const LFSV_S & lfsv_s)
      {}
      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void onBindLFSUVOutside(const IG & ig,
                              const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                              const LFSU_N & lfsu_n, const LFSV_N & lfsv_n)
      {}
      //! @}

      //! Methods for loading of the local function's coefficients.
      //! These methods are empty. The loading of the coefficients
      //! is done after setting the solution in the assembleUVVolume(),
      //! assembleUVSkeleton(), assembleUVBoundary(), assembleUVProcessor(),
      //! assembleUVEnrichedCoupling() and assembleUVVolumePostSkeleton() calls.
      //! @{
      template<typename LFSU>
      void loadCoefficientsLFSUInside(const LFSU & lfsu_s)
      {}
      template<typename LFSU>
      void loadCoefficientsLFSUOutside(const LFSU & lfsu_n)
      {}
      template<typename LFSU>
      void loadCoefficientsLFSUCoupling(const LFSU & lfsu_c)
      {}
      //! @}

      //! Method setting time for la1 local assembler.
      //! This function must be called for explicit methods
      //! before jacobian_engine->assemble.. was called
      void setTimeInLastStage()
      {
        la.la1.setTime(la.time+la.osp_method->d(la.stage)*la.dt);
      }

      //! Notifier functions, called immediately before and after assembling
      //! @{
      void preAssembly()
      {
        lae0->preAssembly();
        lae1->preAssembly();

        *const_residual_0 = 0.0;
        *const_residual_1 = 0.0;

        // Extract the coefficients of the time step scheme
        a.resize(la.stage);
        b.resize(la.stage);
        d.resize(la.stage);
        do0.resize(la.stage);
        do1.resize(la.stage);
        for (int i=0; i<la.stage; ++i){
          a[i] = la.osp_method->a(la.stage,i);
          b[i] = la.osp_method->b(la.stage,i);
          d[i] = la.osp_method->d(i);
          do0[i] = ( std::abs(b[i]) > 1E-6 );
          do1[i] = ( std::abs(a[i]) > 1E-6 );
        }

        // prepare local operators for stage
        la.la0.preStage(la.time+la.osp_method->d(la.stage)*la.dt,la.stage);
        la.la1.preStage(la.time+la.osp_method->d(la.stage)*la.dt,la.stage);
      }

      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
      {
        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
      }
      //! @}

      //! @ Assembling methods
      //! @{

      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        for (int s=0; s<la.stage; ++s)
        {
          // Reset the time in the local assembler
          la.la0.setTime(la.time+d[s]*la.dt);
          la.la1.setTime(la.time+d[s]*la.dt);

          if(do0[s])
          {
            lae0->setSolution(*((*solutions)[s]));
            lae0->onBindLFSUV(eg,lfsu,lfsv);
            lae0->loadCoefficientsLFSUInside(lfsu);
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->assembleUVVolume(eg,lfsu,lfsv);
          }

          if(do1[s])
          {
            lae1->setSolution(*((*solutions)[s]));
            lae1->onBindLFSUV(eg,lfsu,lfsv);
            lae1->loadCoefficientsLFSUInside(lfsu);
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->assembleUVVolume(eg,lfsu,lfsv);
          }
        }
      }

      template<typename EG, typename LFSV>
      void assembleVVolume(const EG & eg, const LFSV & lfsv)
      {
        for (int s=0; s<la.stage; ++s)
        {
          // Reset the time in the local assembler
          la.la0.setTime(la.time+d[s]*la.dt);
          la.la1.setTime(la.time+d[s]*la.dt);

          if(do0[s])
          {
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->assembleVVolume(eg,lfsv);
          }

          if(do1[s])
          {
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->assembleVVolume(eg,lfsv);
          }

        }
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void assembleUVSkeleton(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                              const LFSU_N & lfsu_n, const LFSV_N & lfsv_n)
      {
        for (int s=0; s<la.stage; ++s)
        {
          // Reset the time in the local assembler
          la.la0.setTime(la.time+d[s]*la.dt);
          la.la1.setTime(la.time+d[s]*la.dt);

          if(do0[s])
          {
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->setSolution(*((*solutions)[s]));
            lae0->onBindLFSUVInside(ig,lfsu_s,lfsv_s);
            lae0->onBindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
            lae0->loadCoefficientsLFSUInside(lfsu_s);
            lae0->loadCoefficientsLFSUOutside(lfsu_n);
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
          }

          if(do1[s])
          {
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->setSolution(*((*solutions)[s]));
            lae1->onBindLFSUVInside(ig,lfsu_s,lfsv_s);
            lae1->onBindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
            lae1->loadCoefficientsLFSUInside(lfsu_s);
            lae1->loadCoefficientsLFSUOutside(lfsu_n);
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
          }
        }
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void assembleVSkeleton(const IG & ig, const LFSV_S & lfsv_s, const LFSV_N & lfsv_n)
      {
        for (int s=0; s<la.stage; ++s)
        {
          // Reset the time in the local assembler
          la.la0.setTime(la.time+d[s]*la.dt);
          la.la1.setTime(la.time+d[s]*la.dt);

          if(do0[s])
          {
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->assembleVSkeleton(ig,lfsv_s,lfsv_n);
          }

          if(do1[s])
          {
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->assembleVSkeleton(ig,lfsv_s,lfsv_n);
          }
        }
      }

      template<typename IG, typename LFSU_S, typename LFSV_S>
      void assembleUVBoundary(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s)
      {
        for (int s=0; s<la.stage; ++s)
        {
          // Reset the time in the local assembler
          la.la0.setTime(la.time+d[s]*la.dt);
          la.la1.setTime(la.time+d[s]*la.dt);

          if(do0[s])
          {
            lae0->setSolution(*((*solutions)[s]));
            lae0->onBindLFSUVInside(ig,lfsu_s,lfsv_s);
            lae0->loadCoefficientsLFSUInside(lfsu_s);
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->assembleUVBoundary(ig,lfsu_s,lfsv_s);
          }

          if(do1[s])
          {
            lae1->setSolution(*((*solutions)[s]));
            lae1->onBindLFSUVInside(ig,lfsu_s,lfsv_s);
            lae1->loadCoefficientsLFSUInside(lfsu_s);
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->assembleUVBoundary(ig,lfsu_s,lfsv_s);
          }
        }
      }

      template<typename IG, typename LFSV_S>
      void assembleVBoundary(const IG & ig, const LFSV_S & lfsv_s)
      {
        for (int s=0; s<la.stage; ++s)
        {
          // Reset the time in the local assembler
          la.la0.setTime(la.time+d[s]*la.dt);
          la.la1.setTime(la.time+d[s]*la.dt);

          if(do0[s])
          {
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->assembleVBoundary(ig,lfsv_s);
          }

          if(do1[s])
          {
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->assembleVBoundary(ig,lfsv_s);
          }
        }
      }

      template<typename IG, typename LFSU_S, typename LFSV_S>
      void assembleUVProcessor(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s)
      {
        for (int s=0; s<la.stage; ++s)
        {
          // Reset the time in the local assembler
          la.la0.setTime(la.time+d[s]*la.dt);
          la.la1.setTime(la.time+d[s]*la.dt);

          if(do0[s])
          {
            lae0->setSolution(*((*solutions)[s]));
            lae0->onBindLFSUVInside(ig,lfsu_s,lfsv_s);
            lae0->loadCoefficientsLFSUInside(lfsu_s);
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->assembleUVProcessor(ig,lfsu_s,lfsv_s);
          }

          if(do1[s])
          {
            lae1->setSolution(*((*solutions)[s]));
            lae1->onBindLFSUVInside(ig,lfsu_s,lfsv_s);
            lae1->loadCoefficientsLFSUInside(lfsu_s);
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->assembleUVProcessor(ig,lfsu_s,lfsv_s);
          }
        }
      }

      template<typename IG, typename LFSV_S>
      void assembleVProcessor(const IG & ig, const LFSV_S & lfsv_s)
      {
        for (int s=0; s<la.stage; ++s)
        {
          // Reset the time in the local assembler
          la.la0.setTime(la.time+d[s]*la.dt);
          la.la1.setTime(la.time+d[s]*la.dt);

          if(do0[s])
          {
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->assembleVProcessor(ig,lfsv_s);
          }

          if(do1[s])
          {
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->assembleVProcessor(ig,lfsv_s);
          }
        }
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N,
               typename LFSU_C, typename LFSV_C>
      void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                                             const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                                             const LFSU_C & lfsu_c, const LFSV_C & lfsv_c)
      {
        for (int s=0; s<la.stage; ++s)
        {
          // Reset the time in the local assembler
          la.la0.setTime(la.time+d[s]*la.dt);
          la.la1.setTime(la.time+d[s]*la.dt);

          if(do0[s])
          {
            lae0->setSolution(*((*solutions)[s]));
            lae0->loadCoefficientsLFSUInside(lfsu_s);
            lae0->loadCoefficientsLFSUOutside(lfsu_n);
            lae0->loadCoefficientsLFSUCoupling(lfsu_c);
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_c,lfsv_c);
          }

          if(do1[s])
          {
            lae1->setSolution(*((*solutions)[s]));
            lae1->loadCoefficientsLFSUInside(lfsu_s);
            lae1->loadCoefficientsLFSUOutside(lfsu_n);
            lae1->loadCoefficientsLFSUCoupling(lfsu_c);
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_c,lfsv_c);
          }
        }
      }

      template<typename IG, typename LFSV_S, typename LFSV_N, typename LFSV_C>
      void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSV_S & lfsv_s,
                                            const LFSV_N & lfsv_n,
                                            const LFSV_C & lfsv_c)
      {
        for (int s=0; s<la.stage; ++s)
        {
          // Reset the time in the local assembler
          la.la0.setTime(la.time+d[s]*la.dt);
          la.la1.setTime(la.time+d[s]*la.dt);

          if(do0[s])
          {
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_c);
          }

          if(do1[s])
          {
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_c);
          }

        }
      }

      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        for (int s=0; s<la.stage; ++s)
        {
          // Reset the time in the local assembler
          la.la0.setTime(la.time+d[s]*la.dt);
          la.la1.setTime(la.time+d[s]*la.dt);

          if(do0[s])
          {
            lae0->setSolution(*((*solutions)[s]));
            lae0->onBindLFSUV(eg,lfsu,lfsv);
            lae0->loadCoefficientsLFSUInside(lfsu);
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
          }

          if(do1[s])
          {
            lae1->setSolution(*((*solutions)[s]));
            lae1->onBindLFSUV(eg,lfsu,lfsv);
            lae1->loadCoefficientsLFSUInside(lfsu);
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
          }

        }
      }

      template<typename EG, typename LFSV>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv)
      {
        for (int s=0; s<la.stage; ++s)
        {
          // Reset the time in the local assembler
          la.la0.setTime(la.time+d[s]*la.dt);
          la.la1.setTime(la.time+d[s]*la.dt);

          if(do0[s])
          {
            la.la0.setWeight(b[s]*la.dt_factor0);
            lae0->assembleVVolumePostSkeleton(eg,lfsv);
          }

          if(do1[s])
          {
            la.la1.setWeight(a[s]*la.dt_factor1);
            lae1->assembleVVolumePostSkeleton(eg,lfsv);
          }
        }
      }
      //! @}

    private:

      //! Default value indicating an invalid residual pointer
      Residual * const invalid_residual;

      //! Default value indicating an invalid solution pointer
      Solutions * const invalid_solutions;

      //! Pointer to the current constant part residual vector in
      //! which to assemble the residual corresponding to the operator
      //! representing the time derivative of order zero and one.
      //! @{
      Residual * const_residual_0;
      Residual * const_residual_1;
      //! @}

      //! Pointer to the current residual vector in which to assemble
      const Solutions * solutions;

      //! Coefficients of time stepping scheme
      std::vector<Real> a;
      std::vector<Real> b;
      std::vector<Real> d;
      std::vector<bool> do0;
      std::vector<bool> do1;

    }; // End of class OneStepLocalPreStageAssemblerEngine

  }
}
#endif // DUNE_PDELAB_GRIDOPERATOR_ONESTEP_PRESTAGEENGINE_HH
