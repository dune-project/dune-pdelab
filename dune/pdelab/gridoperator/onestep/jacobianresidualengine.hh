#ifndef DUNE_PDELAB_GRIDOPERATOR_ONESTEP_JACOBIANRESIDUALENGINE_HH
#define DUNE_PDELAB_GRIDOPERATOR_ONESTEP_JACOBIANRESIDUALENGINE_HH

#include <cassert>

#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>

namespace Dune{
  namespace PDELab{

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

      template<typename TrialConstraintsContainer, typename TestConstraintsContainer>
      bool needsConstraintsCaching(const TrialConstraintsContainer& cu, const TestConstraintsContainer& cv) const
      {
        return cu.containsNonDirichletConstraints() or cv.containsNonDirichletConstraints();
      }

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
          prestage_engine(nullptr),
          jacobian_engine(nullptr)
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
      {
        return prestage_engine->requireSkeleton() or jacobian_engine->requireSkeleton();
      }
      bool requireSkeletonTwoSided() const
      {
        return prestage_engine->requireSkeletonTwoSided() or jacobian_engine->requireSkeletonTwoSided();
      }
      bool requireUVVolume() const
      {
        return prestage_engine->requireUVVolume() or jacobian_engine->requireUVVolume();
      }
      bool requireVVolume() const
      {
        return prestage_engine->requireVVolume() or jacobian_engine->requireVVolume();
      }
      bool requireUVSkeleton() const
      {
        return prestage_engine->requireUVSkeleton() or jacobian_engine->requireUVSkeleton();
      }
      bool requireVSkeleton() const
      {
        return prestage_engine->requireVSkeleton() or jacobian_engine->requireVSkeleton();
      }
      bool requireUVBoundary() const
      {
        return prestage_engine->requireUVBoundary() or jacobian_engine->requireUVBoundary();
      }
      bool requireVBoundary() const
      {
        return prestage_engine->requireVBoundary() or jacobian_engine->requireVBoundary();
      }
      bool requireUVVolumePostSkeleton() const
      {
        return prestage_engine->requireUVVolumePostSkeleton() or jacobian_engine->requireUVVolumePostSkeleton();
      }
      bool requireVVolumePostSkeleton() const
      {
        return prestage_engine->requireVVolumePostSkeleton() or jacobian_engine->requireVVolumePostSkeleton();
      }

      //! @}

      //! Public access to the wrapping local assembler
      const LocalAssembler & localAssembler() const
      {
        return la;
      }

      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG, typename LFSU, typename LFSV>
      void onBindLFSUV(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        if (prestage_engine->requireUVVolume() or prestage_engine->requireVVolume())
          prestage_engine->onBindLFSUV(eg,lfsu,lfsv);
        if (jacobian_engine->requireUVVolume() or jacobian_engine->requireVVolume())
          jacobian_engine->onBindLFSUV(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void onBindLFSV(const EG & eg, const LFSV & lfsv)
      {
        if (prestage_engine->requireUVVolume() or prestage_engine->requireVVolume())
          prestage_engine->onBindLFSV(eg,lfsv);
        if (jacobian_engine->requireUVVolume() or jacobian_engine->requireVVolume())
          jacobian_engine->onBindLFSV(eg,lfsv);
      }

      template<typename IG, typename LFSU, typename LFSV>
      void onBindLFSUVInside(const IG & ig, const LFSU & lfsu, const LFSV & lfsv)
      {
        if (prestage_engine->requireUVVolume() or prestage_engine->requireVVolume())
          prestage_engine->onBindLFSUVInside(ig,lfsu,lfsv);
        if (jacobian_engine->requireUVVolume() or jacobian_engine->requireVVolume())
          jacobian_engine->onBindLFSUVInside(ig,lfsu,lfsv);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void onBindLFSUVOutside(const IG & ig,
                              const LFSU_S & lfsus, const LFSV_S & lfsvs,
                              const LFSU_N & lfsun, const LFSV_N & lfsvn)
      {
        if (prestage_engine->requireSkeleton())
          prestage_engine->onBindLFSUVOutside(ig,lfsus,lfsvs,lfsun,lfsvn);
        if (jacobian_engine->requireSkeleton())
          jacobian_engine->onBindLFSUVOutside(ig,lfsus,lfsvs,lfsun,lfsvn);
      }

      template<typename IG, typename LFSV>
      void onBindLFSVInside(const IG & ig, const LFSV & lfsv)
      {
        if (prestage_engine->requireUVVolume() or prestage_engine->requireVVolume())
          prestage_engine->onBindLFSVInside(ig,lfsv);
        if (jacobian_engine->requireUVVolume() or jacobian_engine->requireVVolume())
          jacobian_engine->onBindLFSVInside(ig,lfsv);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void onBindLFSVOutside(const IG & ig,
                             const LFSV_S & lfsvs,
                             const LFSV_N & lfsvn)
      {
        if (prestage_engine->requireSkeleton())
          prestage_engine->onBindLFSVOutside(ig,lfsvs,lfsvn);
        if (jacobian_engine->requireSkeleton())
          jacobian_engine->onBindLFSVOutside(ig,lfsvs,lfsvn);
      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG, typename LFSV>
      void onUnbindLFSV(const EG & eg, const LFSV & lfsv)
      {
        if (prestage_engine->requireUVVolume() or prestage_engine->requireVVolume())
          prestage_engine->onUnbindLFSV(eg,lfsv);
        if (jacobian_engine->requireUVVolume() or jacobian_engine->requireVVolume())
          jacobian_engine->onUnbindLFSV(eg,lfsv);
      }

      template<typename IG, typename LFSV>
      void onUnbindLFSVInside(const IG & ig, const LFSV & lfsv)
      {
        if (prestage_engine->requireUVVolume() or prestage_engine->requireVVolume())
          prestage_engine->onUnbindLFSVInside(ig,lfsv);
        if (jacobian_engine->requireUVVolume() or jacobian_engine->requireVVolume())
          jacobian_engine->onUnbindLFSVInside(ig,lfsv);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void onUnbindLFSVOutside(const IG & ig,
                               const LFSV_S & lfsvs,
                               const LFSV_N & lfsvn)
      {
        if (prestage_engine->requireSkeleton())
          prestage_engine->onUnbindLFSVOutside(ig,lfsvs,lfsvn);
        if (jacobian_engine->requireSkeleton())
          jacobian_engine->onUnbindLFSVOutside(ig,lfsvs,lfsvn);
      }

      template<typename EG, typename LFSU, typename LFSV>
      void onUnbindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
      {
        if (prestage_engine->requireUVVolume() or prestage_engine->requireVVolume())
          prestage_engine->onUnbindLFSUV(eg,lfsu,lfsv);
        if (jacobian_engine->requireUVVolume() or jacobian_engine->requireVVolume())
          jacobian_engine->onUnbindLFSUV(eg,lfsu,lfsv);
      }

      template<typename IG, typename LFSU, typename LFSV>
      void onUnbindLFSUVInside(const IG& ig, const LFSU& lfsu, const LFSV& lfsv)
      {
        if (prestage_engine->requireUVVolume() or prestage_engine->requireVVolume())
          prestage_engine->onUnbindLFSUVInside(ig,lfsu,lfsv);
        if (jacobian_engine->requireUVVolume() or jacobian_engine->requireVVolume())
          jacobian_engine->onUnbindLFSUVInside(ig,lfsu,lfsv);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void onUnbindLFSUVOutside(const IG& ig,
                                const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                                const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
      {
        if (prestage_engine->requireSkeleton())
          prestage_engine->onUnbindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
        if (jacobian_engine->requireSkeleton())
          jacobian_engine->onUnbindLFSUVOutside(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
      }


      //! @}

      //! Methods for loading of the local function's coefficients
      //! @{
      template<typename LFSU>
      void loadCoefficientsLFSUInside(const LFSU & lfsu_s)
      {
        if (prestage_engine->requireUVVolume())
          prestage_engine->loadCoefficientsLFSUInside(lfsu_s);
        if (jacobian_engine->requireUVVolume())
          jacobian_engine->loadCoefficientsLFSUInside(lfsu_s);
      }

      template<typename LFSU>
      void loadCoefficientsLFSUOutside(const LFSU & lfsu_n)
      {
        if (prestage_engine->requireUVSkeleton())
          prestage_engine->loadCoefficientsLFSUOutside(lfsu_n);
        if (jacobian_engine->requireUVSkeleton())
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

      void preAssembly()
      {
        prestage_engine->preAssembly();
        jacobian_engine->preAssembly();
      }

      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
      {
        prestage_engine->postAssembly(gfsu,gfsv);
        jacobian_engine->postAssembly(gfsu,gfsv);
      }

      //! @}

      //! Assembling methods
      //! @{

      template<typename EG>
      bool assembleCell(const EG & eg)
      {
        const bool abort_a = prestage_engine->assembleCell(eg);
        const bool abort_c = jacobian_engine->assembleCell(eg);
        return abort_a and abort_c;
      }

      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        if (prestage_engine->requireUVVolume())
          prestage_engine->assembleUVVolume(eg,lfsu,lfsv);
        la.setWeight(-1.0);
        prestage_engine->setTimeInLastStage();
        if (jacobian_engine->requireUVVolume())
          jacobian_engine->assembleUVVolume(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void assembleVVolume(const EG & eg, const LFSV & lfsv)
      {
        if (prestage_engine->requireVVolume())
          prestage_engine->assembleVVolume(eg,lfsv);
        la.setWeight(-1.0);
        prestage_engine->setTimeInLastStage();
        if (jacobian_engine->requireVVolume())
          jacobian_engine->assembleVVolume(eg,lfsv);
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void assembleUVSkeleton(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                              const LFSU_N & lfsu_n, const LFSV_N & lfsv_n)
      {
        if (prestage_engine->requireUVSkeleton())
          prestage_engine->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
        la.setWeight(-1.0);
        if (jacobian_engine->requireUVSkeleton())
          jacobian_engine->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void assembleVSkeleton(const IG & ig, const LFSV_S & lfsv_s, const LFSV_N & lfsv_n)
      {
        if (prestage_engine->requireVSkeleton())
          prestage_engine->assembleVSkeleton(ig,lfsv_s,lfsv_n);
        la.setWeight(-1.0);
        if (jacobian_engine->requireVSkeleton())
          jacobian_engine->assembleVSkeleton(ig,lfsv_s,lfsv_n);
      }

      template<typename IG, typename LFSU, typename LFSV>
      void assembleUVBoundary(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s)
      {
        if (prestage_engine->requireUVBoundary())
          prestage_engine->assembleUVBoundary(ig,lfsu_s,lfsv_s);
        la.setWeight(-1.0);
        if (jacobian_engine->requireUVBoundary())
          jacobian_engine->assembleUVBoundary(ig,lfsu_s,lfsv_s);
      }

      template<typename IG, typename LFSV>
      void assembleVBoundary(const IG & ig, const LFSV & lfsv_s)
      {
        if (prestage_engine->requireVBoundary())
          prestage_engine->assembleVBoundary(ig,lfsv_s);
        la.setWeight(-1.0);
        if (jacobian_engine->requireVBoundary())
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
        if (prestage_engine->requireUVVolumePostSkeleton())
          prestage_engine->assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
        la.setWeight(-1.0);
        if (jacobian_engine->requireUVVolumePostSkeleton())
          jacobian_engine->assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
      }

      template<typename EG, typename LFSV>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv)
      {
        if (prestage_engine->requireVVolumePostSkeleton())
          prestage_engine->assembleVVolumePostSkeleton(eg,lfsv);
        la.setWeight(-1.0);
        if (jacobian_engine->requireVVolumePostSkeleton())
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

  }
}
#endif // DUNE_PDELAB_GRIDOPERATOR_ONESTEP_JACOBIANRESIDUALENGINE_HH
