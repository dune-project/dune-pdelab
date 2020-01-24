#ifndef DUNE_PDELAB_GRIDOPERATOR_FASTDG_JACOBIANAPPLYENGINE_HH
#define DUNE_PDELAB_GRIDOPERATOR_FASTDG_JACOBIANAPPLYENGINE_HH

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/localoperator/callswitch.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The fast DG local assembler engine for DUNE grids which
       assembles the local application of the Jacobian

       \tparam LA The local assembler

    */
    template<typename LA>
    class FastDGLocalJacobianApplyAssemblerEngine
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

      //! Wheter the local operator is linear
      static constexpr bool isLinear = LOP::isLinear;

      //! The type of the result vector
      typedef typename LA::Traits::Range Range;
      typedef typename Range::ElementType RangeElement;

      //! The type of the solution vector
      typedef typename LA::Traits::Domain Domain;
      typedef typename Domain::ElementType DomainElement;

      //! The local function spaces
      typedef typename LA::LFSU LFSU;
      typedef typename LA::LFSUCache LFSUCache;
      typedef typename LFSU::Traits::GridFunctionSpace GFSU;
      typedef typename LA::LFSV LFSV;
      typedef typename LA::LFSVCache LFSVCache;
      typedef typename LFSV::Traits::GridFunctionSpace GFSV;

      typedef typename Domain::template ConstAliasedLocalView<LFSUCache> DomainView;
      typedef typename Range::template AliasedLocalView<LFSVCache> RangeView;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      FastDGLocalJacobianApplyAssemblerEngine(const LocalAssembler & local_assembler_)
        : local_assembler(local_assembler_),
          lop(local_assembler_.localOperator())
      {}

      //! Query methods for the global grid assembler
      //! @{
      bool requireSkeleton() const
      { return local_assembler.doAlphaSkeleton(); }
      bool requireSkeletonTwoSided() const
      { return local_assembler.doSkeletonTwoSided(); }
      bool requireUVVolume() const
      { return local_assembler.doAlphaVolume(); }
      bool requireUVSkeleton() const
      { return local_assembler.doAlphaSkeleton(); }
      bool requireUVBoundary() const
      { return local_assembler.doAlphaBoundary(); }
      bool requireUVVolumePostSkeleton() const
      { return local_assembler.doAlphaVolumePostSkeleton(); }
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

      //! Set current solution vector. Should be called prior to
      //! assembling.
      void setSolution(const Domain & solution_)
      {
        if (isLinear)
          DUNE_THROW(Dune::Exception, "In the linear case the jacobian apply does not depend on the current solution and this method should never be called.");
        global_solution_view_inside.attach(solution_);
        global_solution_view_outside.attach(solution_);
      }

      //! Set current update vector. Should be called prior to
      //! assembling.
      void setUpdate(const Domain & update_)
      {
        global_update_view_inside.attach(update_);
        global_update_view_outside.attach(update_);
      }

      //! Set current result vector. Should be called prior to
      //! assembling.
      void setResult(Range & result_)
      {
        global_result_view_inside.attach(result_);
        global_result_view_outside.attach(result_);
      }

      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG, typename LFSUC, typename LFSVC>
      void onBindLFSUV(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        if (not isLinear)
          global_solution_view_inside.bind(lfsu_cache);
        global_update_view_inside.bind(lfsu_cache);
      }

      template<typename EG, typename LFSVC>
      void onBindLFSV(const EG & eg, const LFSVC & lfsv_cache)
      {
        global_result_view_inside.bind(lfsv_cache);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void onBindLFSUVInside(const IG & ig, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        if (not isLinear)
          global_solution_view_inside.bind(lfsu_cache);
        global_update_view_inside.bind(lfsu_cache);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void onBindLFSUVOutside(const IG & ig,
                              const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                              const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache)
      {
        if (not isLinear)
          global_solution_view_outside.bind(lfsu_n_cache);
        global_update_view_outside.bind(lfsu_n_cache);
      }

      template<typename IG, typename LFSVC>
      void onBindLFSVInside(const IG & ig, const LFSVC & lfsv_cache)
      {
        global_result_view_inside.bind(lfsv_cache);
      }

      template<typename IG, typename LFSVC>
      void onBindLFSVOutside(const IG & ig,
                             const LFSVC & lfsv_s_cache,
                             const LFSVC & lfsv_n_cache)
      {
        global_result_view_outside.bind(lfsv_n_cache);
      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG, typename LFSVC>
      void onUnbindLFSV(const EG & eg, const LFSVC & lfsv_cache)
      {
        global_result_view_inside.commit();
        global_result_view_inside.unbind();
      }

      template<typename IG, typename LFSVC>
      void onUnbindLFSVInside(const IG & ig, const LFSVC & lfsv_cache)
      {
        global_result_view_inside.commit();
        global_result_view_inside.unbind();
      }

      template<typename IG, typename LFSVC>
      void onUnbindLFSVOutside(const IG & ig,
                               const LFSVC & lfsv_s_cache,
                               const LFSVC & lfsv_n_cache)
      {
        global_result_view_outside.commit();
        global_result_view_outside.unbind();
      }
      //! @}

      //! Methods for loading of the local function's coefficients
      //! @{
      template<typename LFSUC>
      void loadCoefficientsLFSUInside(const LFSUC & lfsu_s_cache)
      {}

      template<typename LFSUC>
      void loadCoefficientsLFSUOutside(const LFSUC & lfsu_n_cache)
      {}

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
                                           global_result_view_inside.container());
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
        global_result_view_inside.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolume>::
          jacobian_apply_volume(lop,eg,lfsu_cache.localFunctionSpace(),global_solution_view_inside,global_update_view_inside,lfsv_cache.localFunctionSpace(),global_result_view_inside);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void assembleUVSkeleton(const IG & ig, const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                              const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache)
      {
        global_result_view_inside.setWeight(local_assembler.weight());
        global_result_view_outside.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
          jacobian_apply_skeleton(lop,ig,
                                  lfsu_s_cache.localFunctionSpace(),global_solution_view_inside,global_update_view_inside,lfsv_s_cache.localFunctionSpace(),
                                  lfsu_n_cache.localFunctionSpace(),global_solution_view_outside,global_update_view_outside,lfsv_n_cache.localFunctionSpace(),
                         global_result_view_inside,global_result_view_outside);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void assembleUVBoundary(const IG & ig, const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache)
      {
        global_result_view_inside.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          jacobian_apply_boundary(lop,ig,lfsu_s_cache.localFunctionSpace(),global_solution_view_inside,global_update_view_inside,lfsv_s_cache.localFunctionSpace(),global_result_view_inside);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      static void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                                             const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache,
                                             const LFSUC & lfsu_coupling_cache, const LFSVC & lfsv_coupling_cache)
      {
        DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");
      }

      template<typename EG, typename LFSUC, typename LFSVC>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        global_result_view_inside.setWeight(local_assembler.weight());
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolumePostSkeleton>::
          jacobian_apply_volume_post_skeleton(lop,eg,lfsu_cache.localFunctionSpace(),global_solution_view_inside,global_update_view_inside,lfsv_cache.localFunctionSpace(),global_result_view_inside);
      }

      //! @}

    private:
      //! Reference to the wrapping local assembler object which
      //! constructed this engine
      const LocalAssembler & local_assembler;

      //! Reference to the local operator
      const LOP & lop;

      //! Pointer to the current solution vector
      DomainView global_solution_view_inside;
      DomainView global_solution_view_outside;

      //! Pointer to the current update vector
      DomainView global_update_view_inside;
      DomainView global_update_view_outside;

      //! Pointer to the current result vector in which to assemble
      RangeView global_result_view_inside;
      RangeView global_result_view_outside;

      //! The local vectors and matrices as required for assembling
      //! @{
      typedef Dune::PDELab::TrialSpaceTag LocalTrialSpaceTag;
      typedef Dune::PDELab::TestSpaceTag LocalTestSpaceTag;

      typedef Dune::PDELab::LocalVector<DomainElement, LocalTrialSpaceTag> DomainVector;
      typedef Dune::PDELab::LocalVector<RangeElement, LocalTestSpaceTag> RangeVector;

      //! @}

    }; // End of class FastDGLocalJacobianAssemblerEngine

  }
}
#endif
