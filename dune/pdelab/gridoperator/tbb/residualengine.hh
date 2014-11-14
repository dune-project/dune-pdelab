#ifndef DUNE_PDELAB_TBB_RESIDUALENGINE_HH
#define DUNE_PDELAB_TBB_RESIDUALENGINE_HH

#include <cstddef>
#include <memory>
#include <mutex>

#include <boost/utility.hpp>

#include <tbb/tbb_stddef.h>

#include <dune/common/shared_ptr.hh>

#include <dune/pdelab/backend/common/threadedvectorview.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/gridoperator/default/residualengine.hh>
#include <dune/pdelab/localoperator/callswitch.hh>

namespace Dune{
  namespace PDELab{

    //! \brief simple wrapper about DefaultLocalResidualAssemblerEngine with
    //!        splitting support
    /**
     * \note We don't simply want to copy the default engine or add a
     *       splitting contructor to it since we want to be able to
     *       differentiate between engines with and without tbb support.
     *       I.e. even when coloring is employed in a manner suffiently for
     *       some engines, it may not be sufficient for other engines
     *       (e.g. pattern).
     */
    template<typename LA>
    class ColoredTBBLocalResidualAssemblerEngine :
      private boost::base_from_member<
        const std::shared_ptr<typename LA::LocalOperator> >,
      public DefaultLocalResidualAssemblerEngine<LA>
    {
      typedef typename LA::LocalOperator LOP;
      typedef boost::base_from_member<const std::shared_ptr<LOP> >
        LOPStorage;
      typedef DefaultLocalResidualAssemblerEngine<LA> Base;
    public:
      //! \brief Constructor
      /**
       * \param [in] local_assembler_ The local assembler object which
       *                              creates this engine
       */
      ColoredTBBLocalResidualAssemblerEngine(const LA &local_assembler_) :
        LOPStorage(stackobject_to_shared_ptr(local_assembler_.lop)),
        Base(local_assembler_)
      { }
      ColoredTBBLocalResidualAssemblerEngine
      ( ColoredTBBLocalResidualAssemblerEngine &other, tbb::split) :
        LOPStorage(std::make_shared<LOP>(*other.LOPStorage::member,
                                         tbb::split())),
        Base(other, *LOPStorage::member)
      { }
      void join(ColoredTBBLocalResidualAssemblerEngine &other)
      {
        LOPStorage::member->join(*other.LOPStorage::member);
      }
    };

    /**
       \brief The local assembler engine for DUNE grids which
       assembles the residual vector

       This version is for uses with threaded assemblers.  It collects the
       to-be-updated entries in a buffer and writes them once the commit limit
       is reached, obtaining a lock on the vector to prevent races.

       \tparam LA The local assembler

    */
    template<typename LA>
    class BatchedTBBLocalResidualAssemblerEngine
      : public LocalAssemblerEngineBase
    {
      typedef typename LA::Mutex Mutex;
    public:

      static const bool needs_constraints_caching = true;

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
      typedef typename GFSU::Traits::GridView GridView;

      typedef typename Solution::template ConstLocalView<LFSUCache> SolutionView;
      typedef ThreadedVectorView<Residual, LFSVCache, Mutex> ResidualView;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      BatchedTBBLocalResidualAssemblerEngine
        (const LocalAssembler & local_assembler_)
        : local_assembler(local_assembler_),
          lop(stackobject_to_shared_ptr(local_assembler_.lop)),
          rl_view(rl,1.0),
          rn_view(rn,1.0)
      {}

      BatchedTBBLocalResidualAssemblerEngine
        (BatchedTBBLocalResidualAssemblerEngine &other, tbb::split)
        : local_assembler(other.local_assembler),
          lop(std::make_shared<LOP>(*other.lop, tbb::split())),
          global_rl_view(make_shared<typename ResidualView::Buffer>
                         (*other.global_rl_view.buffer())),
          global_rn_view(global_rl_view.buffer()),
          global_sl_view(other.global_sl_view),
          global_sn_view(other.global_sn_view),
          rl_view(rl,1.0),
          rn_view(rn,1.0)
      { }

      void join(BatchedTBBLocalResidualAssemblerEngine &other)
      {
        other.global_rl_view.commit();
        other.global_rn_view.commit();
        lop->join(*other.lop);
      }

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
      void setResidual(Residual & residual_, Mutex &mutex)
      {
        auto buffer =
          make_shared<typename ResidualView::Buffer>(residual_, mutex);
        global_rl_view.attach(buffer);
        global_rn_view.attach(buffer);
      }

      //! Set current solution vector. Should be called prior to
      //! assembling.
      void setSolution(const Solution & solution_){
        global_sl_view.attach(solution_);
        global_sn_view.attach(solution_);
      }

      //! Set current solution vector. Should be called prior to
      //! assembling.
      void setAutocommit(std::size_t autocommit)
      {
        global_rl_view.buffer()->set_autocommit(autocommit);
      }

      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG, typename LFSUC, typename LFSVC>
      void onBindLFSUV(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache){
        global_sl_view.bind(lfsu_cache);
        xl.resize(lfsu_cache.size());
      }

      template<typename EG, typename LFSVC>
      void onBindLFSV(const EG & eg, const LFSVC & lfsv_cache){
        global_rl_view.bind(lfsv_cache);
        rl.assign(lfsv_cache.size(),0.0);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void onBindLFSUVInside(const IG & ig, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache){
        global_sl_view.bind(lfsu_cache);
        xl.resize(lfsu_cache.size());
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void onBindLFSUVOutside(const IG & ig,
                              const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                              const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache)
      {
        global_sn_view.bind(lfsu_n_cache);
        xn.resize(lfsu_n_cache.size());
      }

      template<typename IG, typename LFSVC>
      void onBindLFSVInside(const IG & ig, const LFSVC & lfsv_cache){
        global_rl_view.bind(lfsv_cache);
        rl.assign(lfsv_cache.size(),0.0);
      }

      template<typename IG, typename LFSVC>
      void onBindLFSVOutside(const IG & ig,
                             const LFSVC & lfsv_s_cache,
                             const LFSVC & lfsv_n_cache)
      {
        global_rn_view.bind(lfsv_n_cache);
        rn.assign(lfsv_n_cache.size(),0.0);
      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG, typename LFSVC>
      void onUnbindLFSV(const EG & eg, const LFSVC & lfsv_cache){
        global_rl_view.add(rl);
      }

      template<typename IG, typename LFSVC>
      void onUnbindLFSVInside(const IG & ig, const LFSVC & lfsv_cache){
        global_rl_view.add(rl);
      }

      template<typename IG, typename LFSVC>
      void onUnbindLFSVOutside(const IG & ig,
                               const LFSVC & lfsv_s_cache,
                               const LFSVC & lfsv_n_cache)
      {
        global_rn_view.add(rn);
      }
      //! @}

      //! Methods for loading of the local function's coefficients
      //! @{
      template<typename LFSUC>
      void loadCoefficientsLFSUInside(const LFSUC & lfsu_s_cache){
        global_sl_view.read(xl);
      }
      template<typename LFSUC>
      void loadCoefficientsLFSUOutside(const LFSUC & lfsu_n_cache){
        global_sn_view.read(xn);
      }
      template<typename LFSUC>
      void loadCoefficientsLFSUCoupling(const LFSUC & lfsu_c_cache)
      {DUNE_THROW(Dune::NotImplemented,"No coupling lfsu available for ");}
      //! @}

      //! Notifier functions, called immediately before and after assembling
      //! @{

      void postAssembly(const GFSU& gfsu, const GFSV& gfsv){
        global_rl_view.commit();
        global_rn_view.commit();

        if(local_assembler.doPostProcessing)
          {
            Dune::PDELab::constrain_residual(*(local_assembler.pconstraintsv),global_rl_view.container());
          }
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
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolume>::
          alpha_volume(*lop, eg, lfsu_cache.localFunctionSpace(), xl,
                       lfsv_cache.localFunctionSpace(), rl_view);
      }

      template<typename EG, typename LFSVC>
      void assembleVVolume(const EG & eg, const LFSVC & lfsv_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaVolume>::
          lambda_volume(*lop, eg, lfsv_cache.localFunctionSpace(), rl_view);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void assembleUVSkeleton(const IG & ig, const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                              const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        rn_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
          alpha_skeleton(*lop, ig,
                         lfsu_s_cache.localFunctionSpace(),xl,lfsv_s_cache.localFunctionSpace(),
                         lfsu_n_cache.localFunctionSpace(),xn,lfsv_n_cache.localFunctionSpace(),
                         rl_view,rn_view);
      }

      template<typename IG, typename LFSVC>
      void assembleVSkeleton(const IG & ig, const LFSVC & lfsv_s_cache, const LFSVC & lfsv_n_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        rn_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaSkeleton>::
          lambda_skeleton(*lop, ig, lfsv_s_cache.localFunctionSpace(),
                          lfsv_n_cache.localFunctionSpace(), rl_view, rn_view);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void assembleUVBoundary(const IG & ig, const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          alpha_boundary(*lop, ig, lfsu_s_cache.localFunctionSpace(), xl,
                         lfsv_s_cache.localFunctionSpace(), rl_view);
      }

      template<typename IG, typename LFSVC>
      void assembleVBoundary(const IG & ig, const LFSVC & lfsv_s_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaBoundary>::
          lambda_boundary(*lop, ig, lfsv_s_cache.localFunctionSpace(),
                          rl_view);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      static void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                                             const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache,
                                             const LFSUC & lfsu_coupling_cache, const LFSVC & lfsv_coupling_cache)
      {DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");}

      template<typename IG, typename LFSVC>
      static void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSVC & lfsv_s_cache,
                                            const LFSVC & lfsv_n_cache,
                                            const LFSVC & lfsv_coupling_cache)
      {DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");}

      template<typename EG, typename LFSUC, typename LFSVC>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolumePostSkeleton>::
          alpha_volume_post_skeleton(*lop, eg, lfsu_cache.localFunctionSpace(),
                                     xl, lfsv_cache.localFunctionSpace(),
                                     rl_view);
      }

      template<typename EG, typename LFSVC>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSVC & lfsv_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaVolumePostSkeleton>::
          lambda_volume_post_skeleton(*lop, eg,
                                      lfsv_cache.localFunctionSpace(),
                                      rl_view);
      }

      //! @}

    private:
      //! Reference to the wrapping local assembler object which
      //! constructed this engine
      const LocalAssembler & local_assembler;

      //! Reference to the local operator
      const shared_ptr<LOP> lop;

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

    }; // End of class BatchedTBBLocalResidualAssemblerEngine

    /**
       \brief The local assembler engine for DUNE grids which
       assembles the residual vector

       \tparam LA The local assembler

    */
    template<typename LA>
    class TBBLocalResidualAssemblerEngine
      : public LocalAssemblerEngineBase
    {
    public:

      static const bool needs_constraints_caching = true;

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
      typedef typename GFSU::Traits::GridView GridView;

      typedef typename Solution::template ConstLocalView<LFSUCache> SolutionView;
      typedef typename Residual::template LocalView<LFSVCache> ResidualView;

      //! Lock manager used for updating the residual
      typedef typename LA::LockManager LockManager;
      typedef typename LockManager::value_type Mutex;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      TBBLocalResidualAssemblerEngine(const LocalAssembler & local_assembler_)
        : local_assembler(local_assembler_),
          lop(stackobject_to_shared_ptr(local_assembler_.lop)),
          rl_view(rl,1.0),
          rn_view(rn,1.0)
      {}

      TBBLocalResidualAssemblerEngine(TBBLocalResidualAssemblerEngine &other,
                                      tbb::split)
        : local_assembler(other.local_assembler),
          lop(std::make_shared<LOP>(*other.lop, tbb::split())),
          lockmgr(other.lockmgr),
          global_rl_view(other.global_rl_view),
          global_rn_view(other.global_rn_view),
          global_sl_view(other.global_sl_view),
          global_sn_view(other.global_sn_view),
          rl_view(rl,1.0),
          rn_view(rn,1.0)
      {}

      void join(TBBLocalResidualAssemblerEngine &other)
      {
        lop->join(*other.lop);
      }

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

      //! Set Lock manager.  Should be called prior to assembling.
      void setLockManager(LockManager & lockmgr_)
      {
        lockmgr = &lockmgr_;
      }

      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG, typename LFSUC, typename LFSVC>
      void onBindLFSUV(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache){
        global_sl_view.bind(lfsu_cache);
        xl.resize(lfsu_cache.size());
      }

      template<typename EG, typename LFSVC>
      void onBindLFSV(const EG & eg, const LFSVC & lfsv_cache){
        global_rl_view.bind(lfsv_cache);
        rl.assign(lfsv_cache.size(),0.0);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void onBindLFSUVInside(const IG & ig, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache){
        global_sl_view.bind(lfsu_cache);
        xl.resize(lfsu_cache.size());
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void onBindLFSUVOutside(const IG & ig,
                              const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                              const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache)
      {
        global_sn_view.bind(lfsu_n_cache);
        xn.resize(lfsu_n_cache.size());
      }

      template<typename IG, typename LFSVC>
      void onBindLFSVInside(const IG & ig, const LFSVC & lfsv_cache){
        global_rl_view.bind(lfsv_cache);
        rl.assign(lfsv_cache.size(),0.0);
      }

      template<typename IG, typename LFSVC>
      void onBindLFSVOutside(const IG & ig,
                             const LFSVC & lfsv_s_cache,
                             const LFSVC & lfsv_n_cache)
      {
        global_rn_view.bind(lfsv_n_cache);
        rn.assign(lfsv_n_cache.size(),0.0);
      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG, typename LFSVC>
      void onUnbindLFSV(const EG & eg, const LFSVC & lfsv_cache){
        std::lock_guard<Mutex> guard((*lockmgr)[eg.entity()]);
        global_rl_view.add(rl);
        global_rl_view.commit();
      }

      template<typename IG, typename LFSVC>
      void onUnbindLFSVInside(const IG & ig, const LFSVC & lfsv_cache){
        std::lock_guard<Mutex> guard((*lockmgr)[*ig.inside()]);
        global_rl_view.add(rl);
        global_rl_view.commit();
      }

      template<typename IG, typename LFSVC>
      void onUnbindLFSVOutside(const IG & ig,
                               const LFSVC & lfsv_s_cache,
                               const LFSVC & lfsv_n_cache)
      {
        std::lock_guard<Mutex> guard((*lockmgr)[*ig.outside()]);
        global_rn_view.add(rn);
        global_rn_view.commit();
      }
      //! @}

      //! Methods for loading of the local function's coefficients
      //! @{
      template<typename LFSUC>
      void loadCoefficientsLFSUInside(const LFSUC & lfsu_s_cache){
        global_sl_view.read(xl);
      }
      template<typename LFSUC>
      void loadCoefficientsLFSUOutside(const LFSUC & lfsu_n_cache){
        global_sn_view.read(xn);
      }
      template<typename LFSUC>
      void loadCoefficientsLFSUCoupling(const LFSUC & lfsu_c_cache)
      {DUNE_THROW(Dune::NotImplemented,"No coupling lfsu available for ");}
      //! @}

      //! Notifier functions, called immediately before and after assembling
      //! @{

      void postAssembly(const GFSU& gfsu, const GFSV& gfsv){
        if(local_assembler.doPostProcessing)
          {
            Dune::PDELab::constrain_residual(*(local_assembler.pconstraintsv),global_rl_view.container());
          }
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
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolume>::
          alpha_volume(*lop, eg, lfsu_cache.localFunctionSpace(), xl,
                       lfsv_cache.localFunctionSpace(), rl_view);
      }

      template<typename EG, typename LFSVC>
      void assembleVVolume(const EG & eg, const LFSVC & lfsv_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaVolume>::
          lambda_volume(*lop, eg, lfsv_cache.localFunctionSpace(), rl_view);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void assembleUVSkeleton(const IG & ig, const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                              const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        rn_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
          alpha_skeleton(*lop, ig,
                         lfsu_s_cache.localFunctionSpace(),xl,lfsv_s_cache.localFunctionSpace(),
                         lfsu_n_cache.localFunctionSpace(),xn,lfsv_n_cache.localFunctionSpace(),
                         rl_view,rn_view);
      }

      template<typename IG, typename LFSVC>
      void assembleVSkeleton(const IG & ig, const LFSVC & lfsv_s_cache, const LFSVC & lfsv_n_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        rn_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaSkeleton>::
          lambda_skeleton(*lop, ig, lfsv_s_cache.localFunctionSpace(),
                          lfsv_n_cache.localFunctionSpace(), rl_view, rn_view);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void assembleUVBoundary(const IG & ig, const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          alpha_boundary(*lop, ig, lfsu_s_cache.localFunctionSpace(), xl,
                         lfsv_s_cache.localFunctionSpace(), rl_view);
      }

      template<typename IG, typename LFSVC>
      void assembleVBoundary(const IG & ig, const LFSVC & lfsv_s_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaBoundary>::
          lambda_boundary(*lop, ig, lfsv_s_cache.localFunctionSpace(),
                          rl_view);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      static void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                                             const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache,
                                             const LFSUC & lfsu_coupling_cache, const LFSVC & lfsv_coupling_cache)
      {DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");}

      template<typename IG, typename LFSVC>
      static void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSVC & lfsv_s_cache,
                                            const LFSVC & lfsv_n_cache,
                                            const LFSVC & lfsv_coupling_cache)
      {DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");}

      template<typename EG, typename LFSUC, typename LFSVC>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolumePostSkeleton>::
          alpha_volume_post_skeleton(*lop, eg, lfsu_cache.localFunctionSpace(),
                                     xl, lfsv_cache.localFunctionSpace(),
                                     rl_view);
      }

      template<typename EG, typename LFSVC>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSVC & lfsv_cache)
      {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaVolumePostSkeleton>::
          lambda_volume_post_skeleton(*lop, eg,
                                      lfsv_cache.localFunctionSpace(),
                                      rl_view);
      }

      //! @}

    private:
      //! Reference to the wrapping local assembler object which
      //! constructed this engine
      const LocalAssembler & local_assembler;

      //! Reference to the local operator
      const shared_ptr<LOP> lop;

      LockManager *lockmgr;

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
#endif // DUNE_PDELAB_TBB_RESIDUALENGINE_HH
