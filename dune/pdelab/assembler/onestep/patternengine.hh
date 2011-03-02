#ifndef DUNE_UDG_PATTERNENGINE_HH
#define DUNE_UDG_PATTERNENGINE_HH
namespace Dune{
  namespace UDG{

    /**
       \brief The local assembler engine for OneStep sub triangulations which
       creates the matrix pattern 

       \tparam LA The local udg assembler

    */
    template<typename OSLA>
    class OneStepLocalPatternAssemblerEngine
    {
    public:
      //! The type of the wrapping local assembler
      typedef OSLA OneStepLocalAssembler;

      typedef OSLA::LocalAssemblerDT0 LocalAssemblerDT0;
      typedef OSLA::LocalAssemblerDT1 LocalAssemblerDT1;

      //! The type of the matrix pattern container
      typedef typename LocalAssemblerDT0::Pattern Pattern;
      typedef Dune::PDELab::LocalSparsityPattern LocalPattern;

      /**
         \brief Constructor 

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      OneStepLocalPatternAssemblerEngine(const LocalAssembler & local_assembler_)
        : local_assembler(local_assembler_), 
          invalid_lae0(static_cast<PatternEngineDT0*>(0)), lae0(invalid_lae0),
          invalid_lae1(static_cast<PatternEngineDT1*>(0)), lae1(invalid_lae1), 
          invalid_pattern(static_cast<Pattern*>(0)), pattern(invalid_pattern)
      {}

      //! Public access to the wrapping local assembler
      const LocalAssembler & localAssembler(){ return local_assembler; }

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setPattern(Pattern & pattern_){

        // Set pointer to global pattern
        pattern = &pattern_;

        // Initialize the engines of the two wrapped local assemblers
        lae0 = & local_assembler.la0.localPatternAssemblerEngine(pattern_);
        lae1 = & local_assembler.la1.localPatternAssemblerEngine(pattern_);
      }

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
        lae0->onBindLFSUInside(ig,lfsu);
        lae1->onBindLFSUInside(ig,lfsu);
      }

      template<typename IG, typename LFSU>
      void onBindLFSUOutside(const IG & ig, const LFSU & lfsun){
        lae0->onBindLFSUOutside(ig,lfsun);
        lae1->onBindLFSUOutside(ig,lfsun);
      }

      template<typename IG, typename LFSV>
      void onBindLFSVInside(const IG & ig, const LFSV & lfsv){
        lae0->onBindLFSVInside(ig,lfsv);
        lae1->onBindLFSVInside(ig,lfsv);
      }

      template<typename IG, typename LFSV>
      void onBindLFSVOutside(const IG & ig, const LFSV & lfsvn){
        lae0->onBindLFSVOutside(ig,lfsvn);
        lae1->onBindLFSVOutside(ig,lfsvn);
      }
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
        lae0->onUnbindLFSUInside(eg,lfsu);
        lae1->onUnbindLFSUInside(eg,lfsu);
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

      //! Notifier functions, called immediately before and after assembling
      //! @{
      void preAssembly(){}
      void postAssembly(){}
      //! @}

      //! Assembling methods
      //! @{
      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        lae0->assembleUVVolume(eg,lfsu,lfsv);
        lae1->assembleUVVolume(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void assembleVVolume(const EG & eg, const LFSV & lfsv)
      {
        lae0->assembleVVolume(eg,lfsv);
        lae1->assembleVVolume(eg,lfsv);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void assembleUVSkeleton(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                              const LFSU_N & lfsu_n, const LFSV_N & lfsv_n)
      {
        lae0->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
        lae1->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void assembleVSkeleton(const IG & ig, const LFSV_S & lfsv_s, const LFSV_N & lfsv_n)
      {
        lae0->assembleVSkeleton(ig,lfsv_s,lfsv_n);
        lae1->assembleVSkeleton(ig,lfsv_s,lfsv_n);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S>
      void assembleUVBoundary(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s)
      {
        lae0->assembleUVBoundary(ig,lfsu_s,lfsv_s);
        lae1->assembleUVBoundary(ig,lfsu_s,lfsv_s);
      }

      template<typename IG, typename LFSV_S>
      void assembleVBoundary(const IG & ig, const LFSV_S & lfsv_s)
      {
        lae0->assembleVBoundary(ig,lfsv_s);
        lae1->assembleVBoundary(ig,lfsv_s);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S,, typename LFSU_N,, typename LFSV_N,
               typename LFSU_C, typename LFSV_C>
      static void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                                             const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                                             const LFSU_C & lfsu_coupling, const LFSV_C & lfsv_coupling)
      {
        lae0->assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_coupling,lfsv_coupling);
        lae1->assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_coupling,lfsv_coupling);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N, typename LFSV_C>
      static void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSV_S & lfsv_s,
                                            const LFSV_N & lfsv_n,
                                            const LFSV_C & lfsv_coupling) 
      {
        lae0->assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_coupling);
        lae1->assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_coupling);
      }

      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        lae0->assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
        lae1->assembleUVVolumePostSkeleton(eg,lfsu,lsfv);
      }

      template<typename EG, typename LFSV>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv)
      {
        lae0->assembleVVolumePostSkeleton(eg,lfsv);
        lae1->assembleVVolumePostSkeleton(eg,lsfv);
      }

      //! @}


    private:

      //! Reference to the wrapping local assembler object
      const LocalAssembler & local_assembler;

      typedef typename LocalAssemblerDT0::LocalPatternAssemblerEngine PatternEngineDT0;
      typedef typename LocalAssemblerDT1::LocalPatternAssemblerEngine PatternEngineDT1;

      PatternEngineDT0 * const invalid_lae0;
      PatternEngineDT0 * lae0;
      PatternEngineDT1 * const invalid_lae1;;
      PatternEngineDT1 * lae1;
      
      //! Default value indicating an invalid solution pointer
      Pattern * const invalid_pattern;

      //! Pointer to the current matrix pattern container
      Pattern * pattern;

    }; // End of class OneStepLocalPatternAssemblerEngine

  };
};
#endif
