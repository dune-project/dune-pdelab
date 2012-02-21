#ifndef DUNE_PDELAB_DEFAULT_RESIDUALENGINE_HH
#define DUNE_PDELAB_DEFAULT_RESIDUALENGINE_HH

#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include <dune/pdelab/constraints/constraints.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for DUNE grids which
       assembles the residual vector

       \tparam LA The local assembler

    */
    template<typename LA>
    class DefaultLocalResidualAssemblerEngine
      : public LocalAssemblerEngineBase
    {
    public:
      //! The type of the wrapping local assembler
      typedef LA LocalAssembler;

      //! The type of the local operator
      typedef typename LA::LocalOperator LOP;

      //! The type of the residual vector
      typedef typename LA::Traits::Residual Residual;
      typedef typename Residual::ElementType ResidualElement;

      //! The type of the solution vector
      typedef typename LA::Traits::Solution Solution;
      typedef typename Solution::ElementType SolutionElement;

      //! The local function spaces
      typedef typename LA::LFSU LFSU;
      typedef typename LA::LFSUCache LFSUCache;
      typedef typename LA::LFSV LFSV;
      typedef typename LA::LFSVCache LFSVCache;

      typedef typename Solution::template ConstLocalView<LFSUCache> SolutionView;
      typedef typename Residual::template LocalView<LFSVCache> ResidualView;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      DefaultLocalResidualAssemblerEngine(const LocalAssembler & local_assembler_)
        : local_assembler(local_assembler_), lop(local_assembler_.lop),
          rl_view(rl,1.0),
          rn_view(rn,1.0)
      {}

      //! Query methods for the global grid assembler
      //! @{
      bool requireSkeleton() const
      { return ( local_assembler.doAlphaSkeleton() || local_assembler.doLambdaSkeleton() ); }
      bool requireSkeletonTwoSided() const
      { return local_assembler.doSkeletonTwoSided(); }
      bool requireUVVolume() const
      { return local_assembler.doAlphaVolume(); }
      bool requireVVolume() const
      { return local_assembler.doLambdaVolume(); }
      bool requireUVSkeleton() const
      { return local_assembler.doAlphaSkeleton(); }
      bool requireVSkeleton() const
      { return local_assembler.doLambdaSkeleton(); }
      bool requireUVBoundary() const
      { return local_assembler.doAlphaBoundary(); }
      bool requireVBoundary() const
      { return local_assembler.doLambdaBoundary(); }
      bool requireUVVolumePostSkeleton() const
      { return local_assembler.doAlphaVolumePostSkeleton(); }
      bool requireVVolumePostSkeleton() const
      { return local_assembler.doLambdaVolumePostSkeleton(); }
      //! @}

      //! Public access to the wrapping local assembler
      const LocalAssembler & localAssembler() const { return local_assembler; }

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setResidual(Residual & residual_){
        global_rl_view.attach(residual_);
        global_rn_view.attach(residual_);
      }

      //! Set current solution vector. Should be called prior to
      //! assembling.
      void setSolution(const Solution & solution_){
        global_sl_view.attach(solution_);
        global_sn_view.attach(solution_);
      }

      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG>
      void onBindLFSUV(const EG & eg, const LFSUCache & lfsu_cache, const LFSVCache & lfsv_cache){
        global_sl_view.bind(lfsu_cache);
        xl.resize(lfsu_cache.size());
      }

      template<typename EG>
      void onBindLFSV(const EG & eg, const LFSVCache & lfsv_cache){
        global_rl_view.bind(lfsv_cache);
        rl.assign(lfsv_cache.size(),0.0);
      }

      template<typename IG>

      void onBindLFSUVInside(const IG & ig, const LFSUCache & lfsu_cache, const LFSVCache & lfsv_cache){
        global_sl_view.bind(lfsu_cache);
        xl.resize(lfsu_cache.size());
      }

      template<typename IG>
      void onBindLFSUVOutside(const IG & ig,
                              const LFSUCache & lfsu_s_cache, const LFSVCache & lfsv_s_cache,
                              const LFSUCache & lfsu_n_cache, const LFSVCache & lfsv_n_cache)
      {
        global_sn_view.bind(lfsu_n_cache);
        xn.resize(lfsu_n_cache.size());
      }

      template<typename IG>
      void onBindLFSVInside(const IG & ig, const LFSVCache & lfsv_cache){
        global_rl_view.bind(lfsv_cache);
        rl.assign(lfsv_cache.size(),0.0);
      }

      template<typename IG>
      void onBindLFSVOutside(const IG & ig,
                             const LFSVCache & lfsv_s_cache,
                             const LFSVCache & lfsv_n_cache)
      {
        global_rn_view.bind(lfsv_n_cache);
        rn.assign(lfsv_n_cache.size(),0.0);
      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG>
      void onUnbindLFSV(const EG & eg, const LFSVCache & lfsv_cache){
        global_rl_view.add(rl);
        global_rl_view.commit();
      }

      template<typename IG>
      void onUnbindLFSVInside(const IG & ig, const LFSVCache & lfsv_cache){
        global_rl_view.add(rl);
        global_rl_view.commit();
      }

      template<typename IG>
      void onUnbindLFSVOutside(const IG & ig,
                               const LFSVCache & lfsv_s_cache,
                               const LFSVCache & lfsv_n_cache)
      {
        global_rn_view.add(rn);
        global_rn_view.commit();
      }
      //! @}

      //! Methods for loading of the local function's coefficients
      //! @{
      void loadCoefficientsLFSUInside(const LFSUCache & lfsu_s_cache){
        global_sl_view.read(xl);
      }
      void loadCoefficientsLFSUOutside(const LFSUCache & lfsu_n_cache){
        global_sn_view.read(xn);
      }
      void loadCoefficientsLFSUCoupling(const LFSUCache & lfsu_c_cache)
      {DUNE_THROW(Dune::NotImplemented,"No coupling lfsu available for ");}
      //! @}

      //! Notifier functions, called immediately before and after assembling
      //! @{

      void postAssembly(){
        /*
        if(local_assembler.doConstraintsPostProcessing){
          Dune::PDELab::constrain_residual(*(local_assembler.pconstraintsv),global_rl_view.global_container());
        }
        */
      }

      //! @}

      //! Assembling methods
      //! @{

      /** Assemble on a given cell without function spaces.

          \return If true, the assembling for this cell is assumed to
          be complete and the assembler continues with the next grid
          cell.
       */
      template<typename EG>
      bool assembleCell(const EG & eg)
      {
        return LocalAssembler::isNonOverlapping && eg.entity().partitionType() != Dune::InteriorEntity;
      }

      template<typename EG>
      void assembleUVVolume(const EG & eg, const LFSUCache & lfsu_cache, const LFSVCache & lfsv_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolume>::
          alpha_volume(lop,eg,lfsu_cache.localFunctionSpace(),xl,lfsv_cache.localFunctionSpace(),rl_view);
      }

      template<typename EG>
      void assembleVVolume(const EG & eg, const LFSVCache & lfsv_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaVolume>::
          lambda_volume(lop,eg,lfsv_cache.localFunctionSpace(),rl_view);
      }

      template<typename IG>
      void assembleUVSkeleton(const IG & ig, const LFSUCache & lfsu_s_cache, const LFSVCache & lfsv_s_cache,
                              const LFSUCache & lfsu_n_cache, const LFSVCache & lfsv_n_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        rn_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
          alpha_skeleton(lop,ig,
                         lfsu_s_cache.localFunctionSpace(),xl,lfsv_s_cache.localFunctionSpace(),
                         lfsu_n_cache.localFunctionSpace(),xn,lfsv_n_cache.localFunctionSpace(),
                         rl_view,rn_view);
      }

      template<typename IG>
      void assembleVSkeleton(const IG & ig, const LFSVCache & lfsv_s_cache, const LFSVCache & lfsv_n_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        rn_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaSkeleton>::
          lambda_skeleton(lop, ig, lfsv_s_cache.localFunctionSpace(), lfsv_n_cache.localFunctionSpace(), rl_view, rn_view);
      }

      template<typename IG>
      void assembleUVBoundary(const IG & ig, const LFSUCache & lfsu_s_cache, const LFSVCache & lfsv_s_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          alpha_boundary(lop,ig,lfsu_s_cache.localFunctionSpace(),xl,lfsv_s_cache.localFunctionSpace(),rl_view);
      }

      template<typename IG>
      void assembleVBoundary(const IG & ig, const LFSVCache & lfsv_s_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaBoundary>::
          lambda_boundary(lop,ig,lfsv_s_cache.localFunctionSpace(),rl_view);
      }

      template<typename IG>
      static void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSUCache & lfsu_s_cache, const LFSVCache & lfsv_s_cache,
                                             const LFSUCache & lfsu_n_cache, const LFSVCache & lfsv_n_cache,
                                             const LFSUCache & lfsu_coupling_cache, const LFSVCache & lfsv_coupling_cache)
      {DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");}

      template<typename IG>
      static void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSVCache & lfsv_s_cache,
                                            const LFSVCache & lfsv_n_cache,
                                            const LFSVCache & lfsv_coupling_cache)
      {DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");}

      template<typename EG>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSUCache & lfsu_cache, const LFSVCache & lfsv_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolumePostSkeleton>::
          alpha_volume_post_skeleton(lop,eg,lfsu_cache.localFunctionSpace(),xl,lfsv_cache.localFunctionSpace(),rl_view);
      }

      template<typename EG>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSVCache & lfsv_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaVolumePostSkeleton>::
          lambda_volume_post_skeleton(lop,eg,lfsv_cache.localFunctionSpace(),rl_view);
      }

      //! @}

    private:
      //! Reference to the wrapping local assembler object which
      //! constructed this engine
      const LocalAssembler & local_assembler;

      //! Reference to the local operator
      const LOP & lop;

      //! Pointer to the current residual vector in which to assemble
      ResidualView global_rl_view;
      ResidualView global_rn_view;

      //! Pointer to the current residual vector in which to assemble
      SolutionView global_sl_view;
      SolutionView global_sn_view;

      //! The local vectors and matrices as required for assembling
      //! @{
      typedef Dune::PDELab::TrialSpaceTag LocalTrialSpaceTag;
      typedef Dune::PDELab::TestSpaceTag LocalTestSpaceTag;

      typedef Dune::PDELab::LocalVector<SolutionElement, LocalTrialSpaceTag> SolutionVector;
      typedef Dune::PDELab::LocalVector<ResidualElement, LocalTestSpaceTag> ResidualVector;

      //! Inside local coefficients
      SolutionVector xl;
      //! Outside local coefficients
      SolutionVector xn;
      //! Inside local residual
      ResidualVector rl;
      //! Outside local residual
      ResidualVector rn;
      //! Inside local residual weighted view
      typename ResidualVector::WeightedAccumulationView rl_view;
      //! Outside local residual weighted view
      typename ResidualVector::WeightedAccumulationView rn_view;
      //! @}

    }; // End of class DefaultLocalResidualAssemblerEngine

  }
}
#endif
