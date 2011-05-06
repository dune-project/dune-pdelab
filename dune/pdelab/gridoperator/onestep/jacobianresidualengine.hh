#ifndef DUNE_PDELAB_ONESTEP_JACOBIANRESIDUALENGINE_HH
#define DUNE_PDELAB_ONESTEP_JACOBIANRESIDUALENGINE_HH

#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for DUNE grids which
       assembles the residual vector

       \tparam LA The local assembler

    */
    template<typename OSLA>
    class OneStepLocalJacobianResidualAssemblerEngine
      : public LocalAssemblerEngineBase
    {
    public:
      //! The type of the wrapping local assembler
      typedef OSLA OneStepLocalAssembler;

      //! The type for real numbers
      typedef typename OSLA::Real Real;

      typedef OSLA LocalAssembler;
      typedef typename LocalAssembler::LocalPreStageAssemblerEngine PreStageEngine;
      typedef typename LocalAssembler::LocalResidualAssemblerEngine ResidualEngine;
      typedef typename LocalAssembler::LocalJacobianAssemblerEngine JacobianEngine;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      OneStepLocalJacobianResidualAssemblerEngine
      (const LocalAssembler & local_assembler_,
       PreStageEngine & prestage_engine_,
       ResidualEngine & residual_engine_,
       JacobianEngine & jacobian_engine_)
        : la(local_assembler_),
          prestage_engine(prestage_engine_),
          residual_engine(residual_engine_),
          jacobian_engine(jacobian_engine_)
      {}

      //! Query methods for the global grid assembler
      //! @{
      bool requireSkeleton() const
      { return prestage_engine.requireSkeleton() || 
          residual_engine.requireSkeleton() || 
          jacobian_engine.requireSkeleton(); }
      bool requireSkeletonTwoSided() const
      { return prestage_engine.requireSkeletonTwoSided() || 
          residual_engine.requireSkeletonTwoSided() || 
          jacobian_engine.requireSkeletonTwoSided(); }
      bool requireUVVolume() const
      { return prestage_engine.requireUVVolume() || 
          residual_engine.requireUVVolume() || 
          jacobian_engine.requireUVVolume(); }
      bool requireVVolume() const
      { return prestage_engine.requireVVolume() || 
          residual_engine.requireVVolume() || 
          jacobian_engine.requireVVolume(); }
      bool requireUVSkeleton() const
      { return prestage_engine.requireUVSkeleton() || 
          residual_engine.requireUVSkeleton() || 
          jacobian_engine.requireUVSkeleton(); }
      bool requireVSkeleton() const
      { return prestage_engine.requireVSkeleton() || 
          residual_engine.requireVSkeleton() || 
          jacobian_engine.requireVSkeleton(); }
      bool requireUVBoundary() const
      { return prestage_engine.requireUVBoundary() || 
          residual_engine.requireUVBoundary() || 
          jacobian_engine.requireUVBoundary(); }
      bool requireVBoundary() const
      { return prestage_engine.requireVBoundary() || 
          residual_engine.requireVBoundary() || 
          jacobian_engine.requireVBoundary(); }
      bool requireUVVolumePostSkeleton() const
      { return prestage_engine.requireUVVolumePostSkeleton() || 
          residual_engine.requireUVVolumePostSkeleton() || 
          jacobian_engine.requireUVVolumePostSkeleton(); }
      bool requireVVolumePostSkeleton() const
      { return prestage_engine.requireVVolumePostSkeleton() || 
          residual_engine.requireVVolumePostSkeleton() || 
          jacobian_engine.requireVVolumePostSkeleton(); }

      //! @}

      //! Public access to the wrapping local assembler
      const LocalAssembler & localAssembler() const { return la; }

      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG, typename LFSU, typename LFSV>
      void onBindLFSUV(const EG & eg, const LFSU & lfsu, const LFSV & lfsv){
        prestage_engine.onBindLFSUV(eg,lfsu,lfsv);
         residual_engine.onBindLFSUV(eg,lfsu,lfsv);
        jacobian_engine.onBindLFSUV(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void onBindLFSV(const EG & eg, const LFSV & lfsv){
        prestage_engine.onBindLFSV(eg,lfsv);
         residual_engine.onBindLFSV(eg,lfsv);
        jacobian_engine.onBindLFSV(eg,lfsv);
      }

      template<typename IG, typename LFSU, typename LFSV>
      void onBindLFSUVInside(const IG & ig, const LFSU & lfsu, const LFSV & lfsv){
        prestage_engine.onBindLFSUVInside(ig,lfsu,lfsv);
         residual_engine.onBindLFSUVInside(ig,lfsu,lfsv);
        jacobian_engine.onBindLFSUVInside(ig,lfsu,lfsv);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void onBindLFSUVOutside(const IG & ig,
                              const LFSU_S & lfsus, const LFSV_S & lfsvs,
                              const LFSU_N & lfsun, const LFSV_N & lfsvn)
      {
        prestage_engine.onBindLFSUVOutside(ig,lfsus,lfsvs,lfsun,lfsvn);
         residual_engine.onBindLFSUVOutside(ig,lfsus,lfsvs,lfsun,lfsvn);
        jacobian_engine.onBindLFSUVOutside(ig,lfsus,lfsvs,lfsun,lfsvn);
      }

      template<typename IG, typename LFSV>
      void onBindLFSVInside(const IG & ig, const LFSV & lfsv){
        prestage_engine.onBindLFSVInside(ig,lfsv);
         residual_engine.onBindLFSVInside(ig,lfsv);
        jacobian_engine.onBindLFSVInside(ig,lfsv);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void onBindLFSVOutside(const IG & ig,
                             const LFSV_S & lfsvs,
                             const LFSV_N & lfsvn)
      {
        prestage_engine.onBindLFSVOutside(ig,lfsvs,lfsvn);
         residual_engine.onBindLFSVOutside(ig,lfsvs,lfsvn);
        jacobian_engine.onBindLFSVOutside(ig,lfsvs,lfsvn);
      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG, typename LFSV>
      void onUnbindLFSV(const EG & eg, const LFSV & lfsv){
        prestage_engine.onUnbindLFSV(eg,lfsv);
         residual_engine.onUnbindLFSV(eg,lfsv);
        jacobian_engine.onUnbindLFSV(eg,lfsv);
      }

      template<typename IG, typename LFSV>
      void onUnbindLFSVInside(const IG & ig, const LFSV & lfsv){
        prestage_engine.onUnbindLFSVInside(ig,lfsv);
         residual_engine.onUnbindLFSVInside(ig,lfsv);
        jacobian_engine.onUnbindLFSVInside(ig,lfsv);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void onUnbindLFSVOutside(const IG & ig,
                               const LFSV_S & lfsvs,
                               const LFSV_N & lfsvn)
      {
        prestage_engine.onUnbindLFSVOutside(ig,lfsvs,lfsvn);
         residual_engine.onUnbindLFSVOutside(ig,lfsvs,lfsvn);
        jacobian_engine.onUnbindLFSVOutside(ig,lfsvs,lfsvn);
      }

      template<typename EG, typename LFSU, typename LFSV>
      void onUnbindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
      {
        prestage_engine.onUnbindLFSUV(eg,lfsu,lfsv);
         residual_engine.onUnbindLFSUV(eg,lfsu,lfsv);
        jacobian_engine.onUnbindLFSUV(eg,lfsu,lfsv);
      }

      template<typename IG, typename LFSU, typename LFSV>
      void onUnbindLFSUVInside(const IG& ig, const LFSU& lfsu, const LFSV& lfsv)
      {
        prestage_engine.onUnbindLFSUVInside(ig,lfsu,lfsv);
         residual_engine.onUnbindLFSUVInside(ig,lfsu,lfsv);
        jacobian_engine.onUnbindLFSUVInside(ig,lfsu,lfsv);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void onUnbindLFSUVOutside(const IG& ig,
                                const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                                const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
      {
        prestage_engine.onUnbindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
         residual_engine.onUnbindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
        jacobian_engine.onUnbindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
      }


      //! @}

      //! Methods for loading of the local function's coefficients
      //! @{
      template<typename LFSU>
      void loadCoefficientsLFSUInside(const LFSU & lfsu_s){
        prestage_engine.loadCoefficientsLFSUInside(lfsu_s);
         residual_engine.loadCoefficientsLFSUInside(lfsu_s);
        jacobian_engine.loadCoefficientsLFSUInside(lfsu_s);
      }
      template<typename LFSU>
      void loadCoefficientsLFSUOutside(const LFSU & lfsu_n){
        prestage_engine.loadCoefficientsLFSUOutside(lfsu_n);
         residual_engine.loadCoefficientsLFSUOutside(lfsu_n);
        jacobian_engine.loadCoefficientsLFSUOutside(lfsu_n);
      }
      template<typename LFSU>
      void loadCoefficientsLFSUCoupling(const LFSU & lfsu_c)
      {
        prestage_engine.loadCoefficientsLFSUCoupling(lfsu_c);
         residual_engine.loadCoefficientsLFSUCoupling(lfsu_c);
        jacobian_engine.loadCoefficientsLFSUCoupling(lfsu_c);
      }
      //! @}

      //! Notifier functions, called immediately before and after assembling
      //! @{

      void preAssembly(){
        prestage_engine.preAssembly();
         residual_engine.preAssembly();
        jacobian_engine.preAssembly();
      }

      void postAssembly(){
        prestage_engine.postAssembly();
         residual_engine.postAssembly();
        jacobian_engine.postAssembly();
      }

      //! @}

      //! Assembling methods
      //! @{

      template<typename EG>
      bool assembleCell(const EG & eg)
      {
        const bool abort_a = prestage_engine.assembleCell(eg);
        const bool abort_b = residual_engine.assembleCell(eg);
        const bool abort_c = jacobian_engine.assembleCell(eg);
        return abort_a && abort_b && abort_c;
      }

      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        prestage_engine.assembleUVVolume(eg,lfsu,lfsv);
         residual_engine.assembleUVVolume(eg,lfsu,lfsv);
        jacobian_engine.assembleUVVolume(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void assembleVVolume(const EG & eg, const LFSV & lfsv)
      {
        prestage_engine.assembleVVolume(eg,lfsv);
         residual_engine.assembleVVolume(eg,lfsv);
        jacobian_engine.assembleVVolume(eg,lfsv);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void assembleUVSkeleton(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                              const LFSU_N & lfsu_n, const LFSV_N & lfsv_n)
      {
        prestage_engine.assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
         residual_engine.assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
        jacobian_engine.assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void assembleVSkeleton(const IG & ig, const LFSV_S & lfsv_s, const LFSV_N & lfsv_n)
      {
        prestage_engine.assembleVSkeleton(ig,lfsv_s,lfsv_n);
         residual_engine.assembleVSkeleton(ig,lfsv_s,lfsv_n);
        jacobian_engine.assembleVSkeleton(ig,lfsv_s,lfsv_n);
      }

      template<typename IG, typename LFSU, typename LFSV>
      void assembleUVBoundary(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s)
      {
        prestage_engine.assembleUVBoundary(ig,lfsu_s,lfsv_s);
         residual_engine.assembleUVBoundary(ig,lfsu_s,lfsv_s);
        jacobian_engine.assembleUVBoundary(ig,lfsu_s,lfsv_s);
      }

      template<typename IG, typename LFSV>
      void assembleVBoundary(const IG & ig, const LFSV & lfsv_s)
      {
        prestage_engine.assembleVBoundary(ig,lfsv_s);
         residual_engine.assembleVBoundary(ig,lfsv_s);
        jacobian_engine.assembleVBoundary(ig,lfsv_s);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N, 
               typename LFSU_C, typename LFSV_C>
      void assembleUVEnrichedCoupling(const IG & ig,
                                      const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                                      const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                                      const LFSU_C & lfsu_coupling, const LFSV_C & lfsv_coupling)
      {
        prestage_engine.assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_coupling,lfsv_coupling);
         
          residual_engine.assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_coupling,lfsv_coupling);
        jacobian_engine.assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_coupling,lfsv_coupling);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N, typename LFSV_C>
      void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSV_S & lfsv_s,
                                            const LFSV_N & lfsv_n,
                                            const LFSV_C & lfsv_coupling)
      {
        prestage_engine.assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_coupling);
         residual_engine.assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_coupling);
        jacobian_engine.assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_coupling);
      }

      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        prestage_engine.assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
         residual_engine.assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
        jacobian_engine.assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv)
      {
        prestage_engine.assembleVVolumePostSkeleton(eg,lfsv);
         residual_engine.assembleVVolumePostSkeleton(eg,lfsv);
        jacobian_engine.assembleVVolumePostSkeleton(eg,lfsv);
      }

      //! @}

    private:

      //! Reference to the wrapping local assembler object which
      //! constructed this engine
      const LocalAssembler & la;

      PreStageEngine & prestage_engine;
      ResidualEngine & residual_engine;
      JacobianEngine & jacobian_engine;

    }; // End of class OneStepJacobianResidualAssemblerEngine

    /**
       \brief The local assembler engine for DUNE grids which
       assembles the residual vector

       \tparam LA The local assembler

    */
    template<typename OSLA>
    class OneStepExplicitLocalJacobianResidualAssemblerEngine
      : public LocalAssemblerEngineBase
    {
    public:
      //! The type of the wrapping local assembler
      typedef OSLA OneStepLocalAssembler;

      //! The type for real numbers
      typedef typename OSLA::Real Real;

      typedef OSLA LocalAssembler;
      typedef typename LocalAssembler::LocalPreStageAssemblerEngine PreStageEngine;
      typedef typename LocalAssembler::LocalResidualAssemblerEngine ResidualEngine;
      typedef typename LocalAssembler::LocalAssemblerDT1::LocalJacobianAssemblerEngine JacobianEngine;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      OneStepExplicitLocalJacobianResidualAssemblerEngine
      (LocalAssembler & local_assembler_)
        : la(local_assembler_),
          prestage_engine(0),
          jacobian_engine(0)
      {}

      void setLocalPreStageEngine(PreStageEngine & prestage_engine_)
      {
        prestage_engine = & prestage_engine_;
      }

      void setLocalJacobianEngine(JacobianEngine & jacobian_engine_)
      {
        jacobian_engine = & jacobian_engine_;
      }

      //! Query methods for the global grid assembler
      //! @{
      bool requireSkeleton() const
      { return prestage_engine->requireSkeleton() || 
          jacobian_engine->requireSkeleton(); }
      bool requireSkeletonTwoSided() const
      { return prestage_engine->requireSkeletonTwoSided() || 
          jacobian_engine->requireSkeletonTwoSided(); }
      bool requireUVVolume() const
      { return prestage_engine->requireUVVolume() || 
          jacobian_engine->requireUVVolume(); }
      bool requireVVolume() const
      { return prestage_engine->requireVVolume() || 
          jacobian_engine->requireVVolume(); }
      bool requireUVSkeleton() const
      { return prestage_engine->requireUVSkeleton() || 
          jacobian_engine->requireUVSkeleton(); }
      bool requireVSkeleton() const
      { return prestage_engine->requireVSkeleton() || 
          jacobian_engine->requireVSkeleton(); }
      bool requireUVBoundary() const
      { return prestage_engine->requireUVBoundary() || 
          jacobian_engine->requireUVBoundary(); }
      bool requireVBoundary() const
      { return prestage_engine->requireVBoundary() || 
          jacobian_engine->requireVBoundary(); }
      bool requireUVVolumePostSkeleton() const
      { return prestage_engine->requireUVVolumePostSkeleton() || 
          jacobian_engine->requireUVVolumePostSkeleton(); }
      bool requireVVolumePostSkeleton() const
      { return prestage_engine->requireVVolumePostSkeleton() || 
          jacobian_engine->requireVVolumePostSkeleton(); }

      //! @}

      //! Public access to the wrapping local assembler
      const LocalAssembler & localAssembler() const { return la; }

      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG, typename LFSU, typename LFSV>
      void onBindLFSUV(const EG & eg, const LFSU & lfsu, const LFSV & lfsv){
        prestage_engine->onBindLFSUV(eg,lfsu,lfsv);
        jacobian_engine->onBindLFSUV(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void onBindLFSV(const EG & eg, const LFSV & lfsv){
        prestage_engine->onBindLFSV(eg,lfsv);
        jacobian_engine->onBindLFSV(eg,lfsv);
      }

      template<typename IG, typename LFSU, typename LFSV>
      void onBindLFSUVInside(const IG & ig, const LFSU & lfsu, const LFSV & lfsv){
        prestage_engine->onBindLFSUVInside(ig,lfsu,lfsv);
        jacobian_engine->onBindLFSUVInside(ig,lfsu,lfsv);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void onBindLFSUVOutside(const IG & ig,
                              const LFSU_S & lfsus, const LFSV_S & lfsvs,
                              const LFSU_N & lfsun, const LFSV_N & lfsvn)
      {
        prestage_engine->onBindLFSUVOutside(ig,lfsus,lfsvs,lfsun,lfsvn);
        jacobian_engine->onBindLFSUVOutside(ig,lfsus,lfsvs,lfsun,lfsvn);
      }

      template<typename IG, typename LFSV>
      void onBindLFSVInside(const IG & ig, const LFSV & lfsv){
        prestage_engine->onBindLFSVInside(ig,lfsv);
        jacobian_engine->onBindLFSVInside(ig,lfsv);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void onBindLFSVOutside(const IG & ig,
                             const LFSV_S & lfsvs,
                             const LFSV_N & lfsvn)
      {
        prestage_engine->onBindLFSVOutside(ig,lfsvs,lfsvn);
        jacobian_engine->onBindLFSVOutside(ig,lfsvs,lfsvn);
      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG, typename LFSV>
      void onUnbindLFSV(const EG & eg, const LFSV & lfsv){
        prestage_engine->onUnbindLFSV(eg,lfsv);
        jacobian_engine->onUnbindLFSV(eg,lfsv);
      }

      template<typename IG, typename LFSV>
      void onUnbindLFSVInside(const IG & ig, const LFSV & lfsv){
        prestage_engine->onUnbindLFSVInside(ig,lfsv);
        jacobian_engine->onUnbindLFSVInside(ig,lfsv);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void onUnbindLFSVOutside(const IG & ig,
                               const LFSV_S & lfsvs,
                               const LFSV_N & lfsvn)
      {
        prestage_engine->onUnbindLFSVOutside(ig,lfsvs,lfsvn);
        jacobian_engine->onUnbindLFSVOutside(ig,lfsvs,lfsvn);
      }

      template<typename EG, typename LFSU, typename LFSV>
      void onUnbindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
      {
        prestage_engine->onUnbindLFSUV(eg,lfsu,lfsv);
        jacobian_engine->onUnbindLFSUV(eg,lfsu,lfsv);
      }

      template<typename IG, typename LFSU, typename LFSV>
      void onUnbindLFSUVInside(const IG& ig, const LFSU& lfsu, const LFSV& lfsv)
      {
        prestage_engine->onUnbindLFSUVInside(ig,lfsu,lfsv);
        jacobian_engine->onUnbindLFSUVInside(ig,lfsu,lfsv);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void onUnbindLFSUVOutside(const IG& ig,
                                const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                                const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
      {
        prestage_engine->onUnbindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
        jacobian_engine->onUnbindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
      }


      //! @}

      //! Methods for loading of the local function's coefficients
      //! @{
      template<typename LFSU>
      void loadCoefficientsLFSUInside(const LFSU & lfsu_s){
        prestage_engine->loadCoefficientsLFSUInside(lfsu_s);
        jacobian_engine->loadCoefficientsLFSUInside(lfsu_s);
      }
      template<typename LFSU>
      void loadCoefficientsLFSUOutside(const LFSU & lfsu_n){
        prestage_engine->loadCoefficientsLFSUOutside(lfsu_n);
        jacobian_engine->loadCoefficientsLFSUOutside(lfsu_n);
      }
      template<typename LFSU>
      void loadCoefficientsLFSUCoupling(const LFSU & lfsu_c)
      {
        prestage_engine->loadCoefficientsLFSUCoupling(lfsu_c);
        jacobian_engine->loadCoefficientsLFSUCoupling(lfsu_c);
      }
      //! @}

      //! Notifier functions, called immediately before and after assembling
      //! @{

      void preAssembly(){
        prestage_engine->preAssembly();
        jacobian_engine->preAssembly();
      }

      void postAssembly(){
        prestage_engine->postAssembly();
        jacobian_engine->postAssembly();
      }

      //! @}

      //! Assembling methods
      //! @{

      template<typename EG>
      bool assembleCell(const EG & eg)
      {
        const bool abort_a = prestage_engine->assembleCell(eg);
        const bool abort_c = jacobian_engine->assembleCell(eg);
        return abort_a && abort_c;
      }

      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        prestage_engine->assembleUVVolume(eg,lfsu,lfsv);
        la.setWeight(-1.0);
        jacobian_engine->assembleUVVolume(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void assembleVVolume(const EG & eg, const LFSV & lfsv)
      {
        prestage_engine->assembleVVolume(eg,lfsv);
        la.setWeight(-1.0);
        jacobian_engine->assembleVVolume(eg,lfsv);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void assembleUVSkeleton(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                              const LFSU_N & lfsu_n, const LFSV_N & lfsv_n)
      {
        prestage_engine->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
        la.setWeight(-1.0);
        jacobian_engine->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void assembleVSkeleton(const IG & ig, const LFSV_S & lfsv_s, const LFSV_N & lfsv_n)
      {
        prestage_engine->assembleVSkeleton(ig,lfsv_s,lfsv_n);
        la.setWeight(-1.0);
        jacobian_engine->assembleVSkeleton(ig,lfsv_s,lfsv_n);
      }

      template<typename IG, typename LFSU, typename LFSV>
      void assembleUVBoundary(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s)
      {
        prestage_engine->assembleUVBoundary(ig,lfsu_s,lfsv_s);
        la.setWeight(-1.0);
        jacobian_engine->assembleUVBoundary(ig,lfsu_s,lfsv_s);
      }

      template<typename IG, typename LFSV>
      void assembleVBoundary(const IG & ig, const LFSV & lfsv_s)
      {
        prestage_engine->assembleVBoundary(ig,lfsv_s);
        la.setWeight(-1.0);
        jacobian_engine->assembleVBoundary(ig,lfsv_s);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N, 
               typename LFSU_C, typename LFSV_C>
      void assembleUVEnrichedCoupling(const IG & ig,
                                      const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                                      const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                                      const LFSU_C & lfsu_coupling, const LFSV_C & lfsv_coupling)
      {
        prestage_engine->assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_coupling,lfsv_coupling);
        la.setWeight(-1.0);
        jacobian_engine->assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_coupling,lfsv_coupling);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N, typename LFSV_C>
      void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSV_S & lfsv_s,
                                            const LFSV_N & lfsv_n,
                                            const LFSV_C & lfsv_coupling)
      {
        prestage_engine->assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_coupling);
        la.setWeight(-1.0);
        jacobian_engine->assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_coupling);
      }

      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        prestage_engine->assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
        la.setWeight(-1.0);
        jacobian_engine->assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv)
      {
        prestage_engine->assembleVVolumePostSkeleton(eg,lfsv);
        la.setWeight(-1.0);
        jacobian_engine->assembleVVolumePostSkeleton(eg,lfsv);
      }

      //! @}

    private:

      //! Reference to the wrapping local assembler object which
      //! constructed this engine
      LocalAssembler & la;

      PreStageEngine * prestage_engine;
      JacobianEngine * jacobian_engine;

    }; // End of class OneStepJacobianResidualAssemblerEngine


  };
};
#endif
