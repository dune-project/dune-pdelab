#ifndef DUNE_PDELAB_GRIDOPERATOR_FASTDG_PATTERNENGINE_HH
#define DUNE_PDELAB_GRIDOPERATOR_FASTDG_PATTERNENGINE_HH

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/localoperator/callswitch.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The fast DG local assembler engine for DUNE grids which creates
       the matrix pattern

       \tparam LA The local assembler

    */
    template<typename LA>
    class FastDGLocalPatternAssemblerEngine
      : public LocalAssemblerEngineBase
    {
    public:

      template<typename TrialConstraintsContainer, typename TestConstraintsContainer>
      bool needsConstraintsCaching(const TrialConstraintsContainer& cu, const TestConstraintsContainer& cv) const
      {
        return cu.containsNonDirichletConstraints() || cv.containsNonDirichletConstraints();
      }

      //! The type of the wrapping local assembler
      typedef LA LocalAssembler;

      //! The type of the local operator
      typedef typename LA::LocalOperator LOP;

      //! The local function spaces
      typedef typename LA::LFSU LFSU;
      typedef typename LA::LFSUCache LFSUCache;
      typedef typename LFSU::Traits::GridFunctionSpace GFSU;
      typedef typename LA::LFSV LFSV;
      typedef typename LA::LFSVCache LFSVCache;
      typedef typename LFSV::Traits::GridFunctionSpace GFSV;

      //! helper classes
      typedef typename LA::Traits::BorderDOFExchanger BorderDOFExchanger;
      typedef typename BorderDOFExchanger::BorderPattern BorderPattern;

      //! The type of the solution vector
      typedef typename LA::Traits::MatrixPattern Pattern;

      typedef Dune::PDELab::LocalSparsityPattern LocalPattern;

      typedef std::size_t size_type;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      FastDGLocalPatternAssemblerEngine(const LocalAssembler & local_assembler_,
                                         shared_ptr<typename LA::Traits::BorderDOFExchanger> border_dof_exchanger)
        : local_assembler(local_assembler_)
        , lop(local_assembler.localOperator())
        , pattern(nullptr)
        , _border_dof_exchanger(border_dof_exchanger)
      {}

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
      void setPattern(Pattern & pattern_)
      {
        pattern = &pattern_;
      }

      //! Query methods for the global grid assembler
      //! @{

      bool requireSkeleton() const
      {
        return local_assembler.doPatternSkeleton();
      }

      bool requireUVVolume() const
      {
        return local_assembler.doPatternVolume();
      }

      bool requireUVSkeleton() const
      {
        return local_assembler.doPatternSkeleton();
      }

      bool requireUVBoundary() const
      {
        return local_assembler.doPatternBoundary();
      }

      bool requireUVVolumePostSkeleton() const
      {
        return local_assembler.doPatternVolumePostSkeleton();
      }

      //! @}

      //! @}

      template<typename LFSVC, typename LFSUC>
      void add_border_pattern(std::true_type, const LFSVC& lfsv_cache, const LFSUC& lfsu_cache,
                              const LocalPattern& p)
      {
        if (local_assembler.reconstructBorderEntries() &&
            !communicationCache().initialized())
          {
            communicationCache().addEntries(lfsv_cache,lfsu_cache,p);
          }
      }

      template<typename LFSVC, typename LFSUC>
      void add_border_pattern(std::false_type, const LFSVC& lfsv_cache, const LFSUC& lfsu_cache,
                              const LocalPattern& p)
      {}

      template<typename LFSVC, typename LFSUC>
      void add_pattern(const LFSVC& lfsv_cache, const LFSUC& lfsu_cache, const LocalPattern& p)
      {
        for (size_type k=0; k<p.size(); ++k)
          local_assembler.add_entry(*pattern,
                                    lfsv_cache,p[k].i(),
                                    lfsu_cache,p[k].j()
                                    );

        add_border_pattern(std::integral_constant<bool,LocalAssembler::isNonOverlapping>(),
                           lfsv_cache,
                           lfsu_cache,
                           p);
      }


      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG, typename LFSUC, typename LFSVC>
      void onUnbindLFSUV(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        add_pattern(lfsv_cache,lfsu_cache,localpattern);
        localpattern.clear();
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void onUnbindLFSUVOutside(const IG& ig,
                                const LFSUC& lfsu_s_cache, const LFSVC& lfsv_s_cache,
                                const LFSUC& lfsu_n_cache, const LFSVC& lfsv_n_cache)
      {
        add_pattern(lfsv_s_cache,lfsu_n_cache,localpattern_sn);
        localpattern_sn.clear();
        add_pattern(lfsv_n_cache,lfsu_s_cache,localpattern_ns);
        localpattern_ns.clear();
      }

      //! @}

      //! Assembling methods
      //! @{

      template<typename EG, typename LFSUC, typename LFSVC>
      void assembleUVVolume(const EG & eg, const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternVolume>::
          pattern_volume(lop,lfsu_cache.localFunctionSpace(),lfsv_cache.localFunctionSpace(),localpattern);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void assembleUVSkeleton(const IG & ig, const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                              const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternSkeleton>::
          pattern_skeleton(lop,
                           lfsu_s_cache.localFunctionSpace(),lfsv_s_cache.localFunctionSpace(),
                           lfsu_n_cache.localFunctionSpace(),lfsv_n_cache.localFunctionSpace(),
                           localpattern_sn, localpattern_ns);
      }

      template<typename IG, typename LFSUC, typename LFSVC>
      void assembleUVBoundary(const IG & ig, const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternBoundary>::
          pattern_boundary(lop,lfsu_s_cache.localFunctionSpace(),lfsv_s_cache.localFunctionSpace(),localpattern);
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
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternVolumePostSkeleton>::
          pattern_volume_post_skeleton(lop,lfsu_cache.localFunctionSpace(),lfsv_cache.localFunctionSpace(),localpattern);
      }


      void postAssembly(const GFSU& gfsu, const GFSV& gfsv){
        post_border_pattern_assembly(std::integral_constant<bool,LocalAssembler::isNonOverlapping>(),
                                     gfsu,
                                     gfsv);
      }

      void post_border_pattern_assembly(std::true_type, const GFSU& gfsu, const GFSV& gfsv)
      {
        if(local_assembler.doPostProcessing() &&
           local_assembler.reconstructBorderEntries())
          {
            communicationCache().finishInitialization();

            typename LA::Traits::BorderDOFExchanger::template PatternExtender<Pattern>
              data_handle(*_border_dof_exchanger,gfsu,gfsv,*pattern);
            gfsv.gridView().communicate(data_handle,
                                        InteriorBorder_InteriorBorder_Interface,
                                        ForwardCommunication);
          }
      }

      void post_border_pattern_assembly(std::false_type, const GFSU& gfsu, const GFSV& gfsv)
      {}

      //! @}


    private:

      typename LA::Traits::BorderDOFExchanger::CommunicationCache&
      communicationCache()
      {
        return _border_dof_exchanger->communicationCache();
      }

      const typename LA::Traits::BorderDOFExchanger::CommunicationCache&
      communicationCache() const
      {
        return _border_dof_exchanger->communicationCache();
      }

      //! Reference to the wrapping local assembler object
      const LocalAssembler & local_assembler;

      //! Reference to the local operator
      const LOP & lop;

      //! Pointer to the current matrix pattern container
      Pattern * pattern;

      //! Local pattern object used in assembling
      LocalPattern localpattern;
      LocalPattern localpattern_sn, localpattern_ns;

      BorderPattern _border_pattern;

      shared_ptr<BorderDOFExchanger> _border_dof_exchanger;

    }; // End of class FastDGLocalPatternAssemblerEngine

  }
}
#endif
