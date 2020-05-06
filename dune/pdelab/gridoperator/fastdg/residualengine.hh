#ifndef DUNE_PDELAB_GRIDOPERATOR_FASTDG_RESIDUALENGINE_HH
#define DUNE_PDELAB_GRIDOPERATOR_FASTDG_RESIDUALENGINE_HH

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/localoperator/callswitch.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The fast DG local assembler engine for DUNE grids which
       assembles the residual vector

       \tparam LA The local assembler

    */
    template<typename LA>
    class FastDGLocalResidualAssemblerEngine
      : public LocalAssemblerEngineBase
    {
    public:

      template<typename TrialConstraintsContainer, typename TestConstraintsContainer>
      bool needsConstraintsCaching(const TrialConstraintsContainer& cu, const TestConstraintsContainer& cv) const
      {
        return false;
      }

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
      typedef typename LFSU::Traits::GridFunctionSpace GFSU;
      typedef typename LA::LFSV LFSV;
      typedef typename LA::LFSVCache LFSVCache;
      typedef typename LFSV::Traits::GridFunctionSpace GFSV;

      typedef typename Solution::template ConstAliasedLocalView<LFSUCache> SolutionView;
      typedef typename Residual::template AliasedLocalView<LFSVCache> ResidualView;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      FastDGLocalResidualAssemblerEngine(const LocalAssembler & local_assembler_)
        : local_assembler(local_assembler_),
          lop(local_assembler_.localOperator())
          //rl_view(rl,1.0),
          //rn_view(rn,1.0)
      {}

      //! copy contructor
      /**
       * \note This does not create an exact copy.  Instead it copies the
       *       global views, such that they point to the same global vector.
       *       Local matrices/vectors are constructed freshly without copying
       *       the content.  Views into the local matrices/vectors are
       *       constructed freshly so they reference the local matrices/vector
       *       in the new object, and are given unit weight.  This essentially
       *       creates an engine object that is not currently bound to any
       *       entity, but otherwise behaves like the object it was contructed
       *       from.
       * \note This constructor is needed to implement splitting constructors
       *       in derived classes.
       */
      FastDGLocalResidualAssemblerEngine
      (const FastDGLocalResidualAssemblerEngine &other) :
        local_assembler(other.local_assembler), lop(other.lop),
        global_rl_view(other.global_rl_view),
        global_rn_view(other.global_rn_view),
        global_sl_view(other.global_sl_view),
        global_sn_view(other.global_sn_view)
        //rl_view(rl,1.0),
        //rn_view(rn,1.0)
      { }

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
      const LocalAssembler & localAssembler() const
      {
        return local_assembler;
      }

      //! Trial space constraints
      const typename LocalAssembler::Traits::TrialGridFunctionSpaceConstraints& trialConstraints() const
      {
        return localAssembler().trialConstraints();
      }

      //! Test space constraints
      const typename LocalAssembler::Traits::TestGridFunctionSpaceConstraints& testConstraints() const
      {
        return localAssembler().testConstraints();
      }

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setResidual(Residual & residual_)
      {
        global_rl_view.attach(residual_);
        global_rn_view.attach(residual_);
      }

      //! Set current solution vector. Should be called prior to
      //! assembling.
      void setSolution(const Solution & solution_)
      {
        global_sl_view.attach(solution_);
        global_sn_view.attach(solution_);
      }

      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG, typename LFSUC, typename LFSVC>
      void onBindLFSUV(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        global_sl_view.bind(lfsu_cache);
        //xl.resize(lfsu_cache.size());
      }

      template<typename EG, typename LFSVC>
      void onBindLFSV(const EG & eg, const LFSVC & lfsv_cache)
      {
        global_rl_view.bind(lfsv_cache);
        //rl.assign(lfsv_cache.size(),0.0);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void onBindLFSUVInside(const IG & ig, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        global_sl_view.bind(lfsu_cache);
        //xl.resize(lfsu_cache.size());
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void onBindLFSUVOutside(const IG & ig,
                              const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                              const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache)
      {
        global_sn_view.bind(lfsu_n_cache);
        //xn.resize(lfsu_n_cache.size());
      }

      template<typename IG, typename LFSVC>
      void onBindLFSVInside(const IG & ig, const LFSVC & lfsv_cache)
      {
        global_rl_view.bind(lfsv_cache);
        //rl.assign(lfsv_cache.size(),0.0);
      }

      template<typename IG, typename LFSVC>
      void onBindLFSVOutside(const IG & ig,
                             const LFSVC & lfsv_s_cache,
                             const LFSVC & lfsv_n_cache)
      {
        global_rn_view.bind(lfsv_n_cache);
        //rn.assign(lfsv_n_cache.size(),0.0);
      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG, typename LFSVC>
      void onUnbindLFSV(const EG & eg, const LFSVC & lfsv_cache)
      {
        //global_rl_view.add(rl);
        global_rl_view.unbind();
        global_rl_view.commit();
      }

      template<typename IG, typename LFSVC>
      void onUnbindLFSVInside(const IG & ig, const LFSVC & lfsv_cache)
      {
        //global_rl_view.add(rl);
        global_rl_view.commit();
        global_rl_view.unnbind();
      }

      template<typename IG, typename LFSVC>
      void onUnbindLFSVOutside(const IG & ig,
                               const LFSVC & lfsv_s_cache,
                               const LFSVC & lfsv_n_cache)
      {
        //global_rn_view.add(rn);
        global_rn_view.commit();
        global_rn_view.unbind();
      }
      //! @}

      //! Methods for loading of the local function's coefficients
      //! @{
      template<typename LFSUC>
      void loadCoefficientsLFSUInside(const LFSUC & lfsu_s_cache)
      {
        //global_sl_view.read(xl);
      }
      template<typename LFSUC>
      void loadCoefficientsLFSUOutside(const LFSUC & lfsu_n_cache)
      {
        //global_sn_view.read(xn);
      }
      template<typename LFSUC>
      void loadCoefficientsLFSUCoupling(const LFSUC & lfsu_c_cache)
      {
        DUNE_THROW(Dune::NotImplemented,"No coupling lfsu available for ");
      }
      //! @}

      //! Notifier functions, called immediately before and after assembling
      //! @{

      void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
      {
        if(local_assembler.doPostProcessing())
          Dune::PDELab::constrain_residual(local_assembler.testConstraints(),
                                           global_rl_view.container());

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

      template<typename EG, typename LFSUC, typename LFSVC>
      void assembleUVVolume(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        //rl_view.setWeight(local_assembler.weight());
        global_rl_view.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolume>::
          alpha_volume(lop,eg,lfsu_cache.localFunctionSpace(),global_sl_view,lfsv_cache.localFunctionSpace(),global_rl_view);
      }

      template<typename EG, typename LFSVC>
      void assembleVVolume(const EG & eg, const LFSVC & lfsv_cache)
      {
        //rl_view.setWeight(local_assembler.weight());
        global_rl_view.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaVolume>::
          lambda_volume(lop,eg,lfsv_cache.localFunctionSpace(),global_rl_view);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void assembleUVSkeleton(const IG & ig, const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                              const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache)
      {
        //rl_view.setWeight(local_assembler.weight());
        //rn_view.setWeight(local_assembler.weight());
        global_rl_view.setWeight(local_assembler.weight());
        global_rn_view.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
          alpha_skeleton(lop,ig,
                         lfsu_s_cache.localFunctionSpace(),global_sl_view,lfsv_s_cache.localFunctionSpace(),
                         lfsu_n_cache.localFunctionSpace(),global_sn_view,lfsv_n_cache.localFunctionSpace(),
                         global_rl_view,global_rn_view);
      }

      template<typename IG, typename LFSVC>
      void assembleVSkeleton(const IG & ig, const LFSVC & lfsv_s_cache, const LFSVC & lfsv_n_cache)
      {
        //rl_view.setWeight(local_assembler.weight());
        //rn_view.setWeight(local_assembler.weight());
        global_rl_view.setWeight(local_assembler.weight());
        global_rn_view.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaSkeleton>::
          lambda_skeleton(lop, ig, lfsv_s_cache.localFunctionSpace(), lfsv_n_cache.localFunctionSpace(), global_rl_view, global_rn_view);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void assembleUVBoundary(const IG & ig, const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache)
      {
        //rl_view.setWeight(local_assembler.weight());
        global_rl_view.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          alpha_boundary(lop,ig,lfsu_s_cache.localFunctionSpace(),global_sl_view,lfsv_s_cache.localFunctionSpace(),global_rl_view);
      }

      template<typename IG, typename LFSVC>
      void assembleVBoundary(const IG & ig, const LFSVC & lfsv_s_cache)
      {
        //rl_view.setWeight(local_assembler.weight());
        global_rl_view.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaBoundary>::
          lambda_boundary(lop,ig,lfsv_s_cache.localFunctionSpace(),global_rl_view);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      static void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                                             const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache,
                                             const LFSUC & lfsu_coupling_cache, const LFSVC & lfsv_coupling_cache)
      {
        DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");
      }

      template<typename IG, typename LFSVC>
      static void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSVC & lfsv_s_cache,
                                            const LFSVC & lfsv_n_cache,
                                            const LFSVC & lfsv_coupling_cache)
      {
        DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");
      }

      template<typename EG, typename LFSUC, typename LFSVC>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        //rl_view.setWeight(local_assembler.weight());
        global_rl_view.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolumePostSkeleton>::
          alpha_volume_post_skeleton(lop,eg,lfsu_cache.localFunctionSpace(),global_sl_view,lfsv_cache.localFunctionSpace(),global_rl_view);
      }

      template<typename EG, typename LFSVC>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSVC & lfsv_cache)
      {
        //rl_view.setWeight(local_assembler.weight());
        global_rl_view.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaVolumePostSkeleton>::
          lambda_volume_post_skeleton(lop,eg,lfsv_cache.localFunctionSpace(),global_rl_view);
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

      //! Inside local residual weighted view
      //typename ResidualVector::WeightedAccumulationView rl_view;
      //! Outside local residual weighted view
      //typename ResidualVector::WeightedAccumulationView rn_view;
      //! @}

    }; // End of class FastDGLocalResidualAssemblerEngine

  }
}
#endif
