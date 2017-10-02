#ifndef DUNE_PDELAB_GRIDOPERATOR_ONESTEP_ENGINEBASE_HH
#define DUNE_PDELAB_GRIDOPERATOR_ONESTEP_ENGINEBASE_HH

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for UDG sub triangulations which
       assembles the residual vector

       \tparam LA The local udg assembler

    */
    template<typename OSLA, typename LAE0, typename LAE1>
    class OneStepLocalAssemblerEngineBase
    {
    public:
      //! The type of the wrapping local assembler
      typedef OSLA OneStepLocalAssembler;

      typedef typename LAE0::Traits Traits;

      template<typename TrialConstraintsContainer, typename TestConstraintsContainer>
      bool needsConstraintsCaching(const TrialConstraintsContainer& cu, const TestConstraintsContainer& cv) const
      {
        return (lae0->needsConstraintsCaching(cu,cv) or lae1->needsConstraintsCaching(cu,cv));
      }


      //! Types of the subordinate assemblers and engines
      //! @{
      typedef typename OSLA::LocalAssemblerDT0 LocalAssemblerDT0;
      typedef typename OSLA::LocalAssemblerDT1 LocalAssemblerDT1;

      typedef LAE0 LocalAssemblerEngineDT0;
      typedef LAE1 LocalAssemblerEngineDT1;
      //! @}

      //! The type for real numbers
      typedef typename OSLA::Real Real;

      typedef OSLA LocalAssembler;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      OneStepLocalAssemblerEngineBase(const LocalAssembler & local_assembler_)
        : invalid_lae0(nullptr),
          invalid_lae1(nullptr),
          la(local_assembler_),
          lae0(invalid_lae0), lae1(invalid_lae1),
          implicit(true)
      {}

      //! Query methods for the global grid assembler
      //! @{
      bool requireSkeleton() const
      { return implicit && (lae0->requireSkeleton() || lae1->requireSkeleton()); }
      bool requireSkeletonTwoSided() const
      { return lae0->requireSkeletonTwoSided() || lae1->requireSkeletonTwoSided(); }
      bool requireUVVolume() const
      { return lae0->requireUVVolume() || lae1->requireUVVolume(); }
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
      bool requireUVProcessor() const
      { return lae0->requireUVProcessor() || lae1->requireUVProcessor(); }
      bool requireVProcessor() const
      { return lae0->requireVProcessor() || lae1->requireVProcessor(); }
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

      LocalAssemblerEngineDT0& localAssemblerEngineDT0()
      {
        return *lae0;
      }

      const LocalAssemblerEngineDT0& localAssemblerEngineDT0() const
      {
        return *lae0;
      }

      LocalAssemblerEngineDT1& localAssemblerEngineDT1()
      {
        return *lae1;
      }

      const LocalAssemblerEngineDT1& localAssemblerEngineDT1() const
      {
        return *lae1;
      }

      auto partition() const
      {
        return localAssemblerEngineDT0().partition();
      }

      void setLocalAssemblerEngineDT0(LocalAssemblerEngineDT0& lae0_)
      {
        lae0 = &lae0_;
      }

      void setLocalAssemblerEngineDT1(LocalAssemblerEngineDT1& lae1_)
      {
        lae1 = &lae1_;
      }

      const typename OneStepLocalAssembler::Traits::TrialGridFunctionSpaceConstraints& trialConstraints() const
      {
        return localAssemblerEngineDT0().trialConstraints();
      }

      const typename OneStepLocalAssembler::Traits::TestGridFunctionSpaceConstraints& testConstraints() const
      {
        return localAssemblerEngineDT0().testConstraints();
      }

      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG, typename LFSU, typename LFSV>
      void onBindLFSUV(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        lae0->onBindLFSUV(eg,lfsu,lfsv);
        lae1->onBindLFSUV(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void onBindLFSV(const EG & eg, const LFSV & lfsv)
      {
        lae0->onBindLFSV(eg,lfsv);
        lae1->onBindLFSV(eg,lfsv);
      }

      template<typename IG, typename LFSU, typename LFSV>
      void onBindLFSUVInside(const IG & ig, const LFSU & lfsu, const LFSV & lfsv)
      {
        lae0->onBindLFSUVInside(ig,lfsu,lfsv);
        lae1->onBindLFSUVInside(ig,lfsu,lfsv);
      }

      template<typename IG,
               typename LFSU_S, typename LFSV_S,
               typename LFSU_N, typename LFSV_N>
      void onBindLFSUVOutside(const IG & ig,
                              const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                              const LFSU_N & lfsu_n, const LFSV_N & lfsv_n)
      {
        lae0->onBindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
        lae1->onBindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
      }

      template<typename IG, typename LFSV>
      void onBindLFSVInside(const IG & ig, const LFSV & lfsv)
      {
        lae0->onBindLFSVInside(ig,lfsv);
        lae1->onBindLFSVInside(ig,lfsv);
      }

      template<typename IG,
               typename LFSV_S,
               typename LFSV_N>
      void onBindLFSVOutside(const IG & ig,
                             const LFSV_S & lfsv_s,
                             const LFSV_N & lfsv_n)
      {
        lae0->onBindLFSVOutside(ig,lfsv_s,lfsv_n);
        lae1->onBindLFSVOutside(ig,lfsv_s,lfsv_n);
      }

      template<typename IG,
               typename LFSU_S, typename LFSV_S,
               typename LFSU_N, typename LFSV_N,
               typename LFSU_C, typename LFSV_C>
      void onBindLFSUVCoupling(const IG & ig,
                               const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                               const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                               const LFSU_C & lfsu_c, const LFSV_C & lfsv_c)
      {
        lae0->onBindLFSUVCoupling(ig,
                                  lfsu_s,lfsv_s,
                                  lfsu_n,lfsv_n,
                                  lfsu_c,lfsv_c);
        lae1->onBindLFSUVCoupling(ig,
                                  lfsu_s,lfsv_s,
                                  lfsu_n,lfsv_n,
                                  lfsu_c,lfsv_c);
      }

      template<typename IG,
               typename LFSV_S,
               typename LFSV_N,
               typename LFSV_C>
      void onBindLFSVCoupling(const IG & ig,
                              const LFSV_S & lfsv_s,
                              const LFSV_N & lfsv_n,
                              const LFSV_C & lfsv_c)
      {
        lae0->onBindLFSVCoupling(ig,lfsv_s,lfsv_n,lfsv_c);
        lae1->onBindLFSVCoupling(ig,lfsv_s,lfsv_n,lfsv_c);
      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG, typename LFSU, typename LFSV>
      void onUnbindLFSUV(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        lae0->onUnbindLFSUV(eg,lfsu, lfsv);
        lae1->onUnbindLFSUV(eg,lfsu, lfsv);
      }

      template<typename EG, typename LFSV>
      void onUnbindLFSV(const EG & eg, const LFSV & lfsv)
      {
        lae0->onUnbindLFSV(eg,lfsv);
        lae1->onUnbindLFSV(eg,lfsv);
      }

      template<typename IG, typename LFSU, typename LFSV>
      void onUnbindLFSUVInside(const IG & ig, const LFSU & lfsu, const LFSV & lfsv)
      {
        lae0->onUnbindLFSUVInside(ig,lfsu, lfsv);
        lae1->onUnbindLFSUVInside(ig,lfsu, lfsv);
      }

      template<typename IG,
               typename LFSU_S, typename LFSV_S,
               typename LFSU_N, typename LFSV_N>
      void onUnbindLFSUVOutside(const IG & ig,
                                const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                                const LFSU_N & lfsu_n, const LFSV_N & lfsv_n)
      {
        lae0->onUnbindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
        lae1->onUnbindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
      }

      template<typename IG, typename LFSV>
      void onUnbindLFSVInside(const IG & ig, const LFSV & lfsv)
      {
        lae0->onUnbindLFSVInside(ig,lfsv);
        lae1->onUnbindLFSVInside(ig,lfsv);
      }

      template<typename IG,
               typename LFSV_S,
               typename LFSV_N>
      void onUnbindLFSVOutside(const IG & ig,
                               const LFSV_S & lfsv_s,
                               const LFSV_N & lfsv_n)
      {
        lae0->onUnbindLFSVOutside(ig,lfsv_s,lfsv_n);
        lae1->onUnbindLFSVOutside(ig,lfsv_s,lfsv_n);
      }

      template<typename IG,
               typename LFSU_S, typename LFSV_S,
               typename LFSU_N, typename LFSV_N,
               typename LFSU_C, typename LFSV_C>
      void onUnbindLFSUVCoupling(const IG & ig,
                                 const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                                 const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                                 const LFSU_C & lfsu_c, const LFSV_C & lfsv_c)
      {
        lae0->onUnbindLFSUVCoupling(ig,
                                    lfsu_s,lfsv_s,
                                    lfsu_n,lfsv_n,
                                    lfsu_c,lfsv_c);
        lae1->onUnbindLFSUVCoupling(ig,
                                    lfsu_s,lfsv_s,
                                    lfsu_n,lfsv_n,
                                    lfsu_c,lfsv_c);
      }

      template<typename IG,
               typename LFSV_S,
               typename LFSV_N,
               typename LFSV_C>
      void onUnbindLFSVCoupling(const IG & ig,
                                const LFSV_S & lfsv_s,
                                const LFSV_N & lfsv_n,
                                const LFSV_C & lfsv_c)
      {
        lae0->onUnbindLFSVCoupling(ig,lfsv_s,lfsv_n,lfsv_c);
        lae1->onUnbindLFSVCoupling(ig,lfsv_s,lfsv_n,lfsv_c);
      }

      //! @}

      //! Methods for loading of the local function's
      //! coefficients.
      //!@{
      template<typename LFSU>
      void loadCoefficientsLFSUInside(const LFSU & lfsu_s)
      {
        lae0->loadCoefficientsLFSUInside(lfsu_s);
        lae1->loadCoefficientsLFSUInside(lfsu_s);
      }
      template<typename LFSU>
      void loadCoefficientsLFSUOutside(const LFSU & lfsu_n)
      {
        lae0->loadCoefficientsLFSUOutside(lfsu_n);
        lae1->loadCoefficientsLFSUOutside(lfsu_n);
      }
      template<typename LFSU>
      void loadCoefficientsLFSUCoupling(const LFSU & lfsu_c)
      {
        lae0->loadCoefficientsLFSUCoupling(lfsu_c);
        lae1->loadCoefficientsLFSUCoupling(lfsu_c);
      }
      //! @}

      //! @name Assembling methods
      //! @{

      template<typename EG>
      bool assembleCell(const EG & eg)
      {
        bool rv = true;
        rv &= lae0->assembleCell(eg);
        rv &= lae1->assembleCell(eg);
        return rv;
      }

      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        if(implicit)
          lae0->assembleUVVolume(eg,lfsu,lfsv);
        lae1->assembleUVVolume(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void assembleVVolume(const EG & eg, const LFSV & lfsv)
      {
        if(implicit)
          lae0->assembleVVolume(eg,lfsv);
        lae1->assembleVVolume(eg,lfsv);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void assembleUVSkeleton(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                              const LFSU_N & lfsu_n, const LFSV_N & lfsv_n)
      {
        if(implicit)
          lae0->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
        lae1->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void assembleVSkeleton(const IG & ig, const LFSV_S & lfsv_s, const LFSV_N & lfsv_n)
      {
        if(implicit)
          lae0->assembleVSkeleton(ig,lfsv_s,lfsv_n);
        lae1->assembleVSkeleton(ig,lfsv_s,lfsv_n);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S>
      void assembleUVBoundary(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s)
      {
        if(implicit)
          lae0->assembleUVBoundary(ig,lfsu_s,lfsv_s);
        lae1->assembleUVBoundary(ig,lfsu_s,lfsv_s);
      }

      template<typename IG, typename LFSV_S>
      void assembleVBoundary(const IG & ig, const LFSV_S & lfsv_s)
      {
        if(implicit)
          lae0->assembleVBoundary(ig,lfsv_s);
        lae1->assembleVBoundary(ig,lfsv_s);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S>
      void assembleUVProcessor(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s)
      {
        if(implicit)
          lae0->assembleUVProcessor(ig,lfsu_s,lfsv_s);
        lae1->assembleUVProcessor(ig,lfsu_s,lfsv_s);
      }

      template<typename IG, typename LFSV_S>
      void assembleVProcessor(const IG & ig, const LFSV_S & lfsv_s)
      {
        if(implicit)
          lae0->assembleVProcessor(ig,lfsv_s);
        lae1->assembleVProcessor(ig,lfsv_s);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N,
               typename LFSU_C, typename LFSV_C>
      void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                                             const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                                             const LFSU_C & lfsu_c, const LFSV_C & lfsv_c)
      {
        if(implicit)
          lae0->assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_c,lfsv_c);
        lae1->assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_c,lfsv_c);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N, typename LFSV_C>
      void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSV_S & lfsv_s,
                                            const LFSV_N & lfsv_n,
                                            const LFSV_C & lfsv_c)
      {
        if(implicit)
          lae0->assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_c);
      }

      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        if(implicit)
          lae0->assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv)
      {
        if(implicit)
          lae0->assembleVVolumePostSkeleton(eg,lfsv);
      }
      //! @}

    private:

      LocalAssemblerEngineDT0 * const invalid_lae0;
      LocalAssemblerEngineDT1 * const invalid_lae1;

    protected:

      //! Reference to the wrapping local assembler object which
      //! constructed this engine
      const LocalAssembler & la;

      LocalAssemblerEngineDT0 * lae0;
      LocalAssemblerEngineDT1 * lae1;

      bool implicit;

    }; // End of class OneStepLocalAssemblerEngineBase

  }
}

#endif // DUNE_PDELAB_GRIDOPERATOR_ONESTEP_ENGINEBASE_HH
