#ifndef DUNE_PDELAB_DEFAULT_PATTERNENGINE_HH
#define DUNE_PDELAB_DEFAULT_PATTERNENGINE_HH

#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include <dune/pdelab/gridoperator/common/borderdofexchanger.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for DUNE grids which creates
       the matrix pattern

       \tparam LA The local assembler

    */
    template<typename LA>
    class DefaultLocalPatternAssemblerEngine
      : public LocalAssemblerEngineBase
    {
    public:

      static const bool needs_constraints_caching = true;

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


      //! The type of the solution vector
      typedef typename LA::Traits::MatrixPattern Pattern;
      typedef typename LA::Traits::BorderPattern BorderPattern;

      typedef Dune::PDELab::LocalSparsityPattern LocalPattern;

      typedef std::size_t size_type;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      DefaultLocalPatternAssemblerEngine(const LocalAssembler & local_assembler_,
                                         shared_ptr<typename LA::Traits::BorderDOFExchanger> border_dof_exchanger)
        : local_assembler(local_assembler_)
        , lop(local_assembler.lop)
        , pattern(nullptr)
        , _border_dof_exchanger(border_dof_exchanger)
      {}

      //! Public access to the wrapping local assembler
      const LocalAssembler & localAssembler() const
      {
        return local_assembler;
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

      void add_border_pattern(true_type, const LFSVCache& lfsv_cache, const LFSUCache& lfsu_cache, const LocalPattern& p)
      {
        if (local_assembler.reconstructBorderEntries() &&
            !communicationCache().initialized())
          {
            for (auto& entry : p)
              {
                const auto& di = lfsv_cache.dof_index(entry.i());
                const auto& dj = lfsu_cache.dof_index(entry.j());

                size_type row_geometry_index = LFSV::Traits::GridFunctionSpace::Ordering::Traits::DOFIndexAccessor::geometryType(di);
                size_type row_entity_index = LFSV::Traits::GridFunctionSpace::Ordering::Traits::DOFIndexAccessor::entityIndex(di);

                size_type col_geometry_index = LFSU::Traits::GridFunctionSpace::Ordering::Traits::DOFIndexAccessor::geometryType(dj);
                size_type col_entity_index = LFSU::Traits::GridFunctionSpace::Ordering::Traits::DOFIndexAccessor::entityIndex(dj);
                // We are only interested in connections between two border entities.
                if (!communicationCache().isBorderEntity(row_geometry_index,row_entity_index) ||
                    !communicationCache().isBorderEntity(col_geometry_index,col_entity_index))
                  continue;

                communicationCache().addEntry(di,col_geometry_index,col_entity_index,dj.treeIndex());
              }
          }
      }

      void add_border_pattern(false_type, const LFSVCache& lfsv_cache, const LFSUCache& lfsu_cache, const LocalPattern& p)
      {}

      void add_pattern(const LFSVCache& lfsv_cache, const LFSUCache& lfsu_cache, const LocalPattern& p)
      {
        for (size_type k=0; k<p.size(); ++k)
          local_assembler.add_entry(*pattern,
                                    lfsv_cache,p[k].i(),
                                    lfsu_cache,p[k].j()
                                    );

        add_border_pattern(integral_constant<bool,LocalAssembler::isNonOverlapping>(),
                           lfsv_cache,
                           lfsu_cache,
                           p);
      }


      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG>
      void onUnbindLFSUV(const EG & eg, const LFSUCache & lfsu_cache, const LFSVCache & lfsv_cache)
      {
        add_pattern(lfsv_cache,lfsu_cache,localpattern);
        localpattern.clear();
      }

      template<typename IG>
      void onUnbindLFSUVOutside(const IG& ig,
                                const LFSUCache& lfsu_s_cache, const LFSVCache& lfsv_s_cache,
                                const LFSUCache& lfsu_n_cache, const LFSVCache& lfsv_n_cache)
      {
        add_pattern(lfsv_s_cache,lfsu_n_cache,localpattern_sn);
        localpattern_sn.clear();
        add_pattern(lfsv_n_cache,lfsu_s_cache,localpattern_ns);
        localpattern_ns.clear();
      }

      //! @}

      //! Assembling methods
      //! @{

      template<typename EG>
      void assembleUVVolume(const EG & eg, const LFSUCache & lfsu_cache, const LFSVCache & lfsv_cache)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternVolume>::
          pattern_volume(lop,lfsu_cache.localFunctionSpace(),lfsv_cache.localFunctionSpace(),localpattern);
      }

      template<typename IG>
      void assembleUVSkeleton(const IG & ig, const LFSUCache & lfsu_s_cache, const LFSVCache & lfsv_s_cache,
                              const LFSUCache & lfsu_n_cache, const LFSVCache & lfsv_n_cache)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternSkeleton>::
          pattern_skeleton(lop,
                           lfsu_s_cache.localFunctionSpace(),lfsv_s_cache.localFunctionSpace(),
                           lfsu_n_cache.localFunctionSpace(),lfsv_n_cache.localFunctionSpace(),
                           localpattern_sn, localpattern_ns);
      }

      template<typename IG>
      void assembleUVBoundary(const IG & ig, const LFSUCache & lfsu_s_cache, const LFSVCache & lfsv_s_cache)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternBoundary>::
          pattern_boundary(lop,lfsu_s_cache,lfsv_s_cache,localpattern);
      }

      template<typename IG>
      static void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSUCache & lfsu_s_cache, const LFSVCache & lfsv_s_cache,
                                             const LFSUCache & lfsu_n_cache, const LFSVCache & lfsv_n_cache,
                                             const LFSUCache & lfsu_coupling_cache, const LFSVCache & lfsv_coupling_cache)
      {
        DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");
      }

      template<typename IG>
      static void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSVCache & lfsv_s_cache,
                                            const LFSVCache & lfsv_n_cache,
                                            const LFSVCache & lfsv_coupling_cache)
      {
        DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");
      }

      template<typename EG>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSUCache & lfsu_cache, const LFSVCache & lfsv_cache)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternVolumePostSkeleton>::
          pattern_volume_post_skeleton(lop,lfsu_cache.localFunctionSpace(),lfsv_cache.localFunctionSpace(),localpattern);
      }


      void postAssembly(const GFSU& gfsu, const GFSV& gfsv){
        post_border_pattern_assembly(integral_constant<bool,LocalAssembler::isNonOverlapping>(),
                                     gfsu,
                                     gfsv);
      }

      void post_border_pattern_assembly(true_type, const GFSU& gfsu, const GFSV& gfsv)
      {
        if(local_assembler.doPostProcessing &&
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

      void post_border_pattern_assembly(false_type, const GFSU& gfsu, const GFSV& gfsv)
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

      shared_ptr<typename LA::Traits::BorderDOFExchanger> _border_dof_exchanger;

    }; // End of class DefaultLocalPatternAssemblerEngine

  }
}
#endif
