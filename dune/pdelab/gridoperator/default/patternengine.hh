#ifndef DUNE_PDELAB_DEFAULT_PATTERNENGINE_HH
#define DUNE_PDELAB_DEFAULT_PATTERNENGINE_HH

#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>

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
      //! The type of the wrapping local assembler
      typedef LA LocalAssembler;

      //! The type of the local operator
      typedef typename LA::LocalOperator LOP;

      //! The local function spaces
      typedef typename LA::LFSU LFSU;
      typedef typename LA::LFSV LFSV;

      //! The type of the solution vector
      typedef typename LA::Traits::MatrixPattern Pattern;

      typedef Dune::PDELab::LocalSparsityPattern LocalPattern;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      DefaultLocalPatternAssemblerEngine(const LocalAssembler & local_assembler_)
        : local_assembler(local_assembler_), lop(local_assembler_.lop),
          invalid_pattern(static_cast<Pattern*>(0)), pattern(invalid_pattern)
      {}

      //! Public access to the wrapping local assembler
      const LocalAssembler & localAssembler() const { return local_assembler; }

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setPattern(Pattern & pattern_){
        pattern = &pattern_;
      }

      //! Query methods for the global grid assembler
      //! @{
      bool requireSkeleton() const
      { return local_assembler.doPatternSkeleton(); }

      bool requireUVVolume() const
      { return local_assembler.doPatternVolume(); }
      bool requireUVSkeleton() const
      { return local_assembler.doPatternSkeleton(); }
      bool requireUVBoundary() const
      { return local_assembler.doPatternBoundary(); }
      bool requireUVVolumePostSkeleton() const
      { return local_assembler.doPatternVolumePostSkeleton(); }
      //! @}

      //! @}


      void add_pattern(const LFSV& lfsv, const LFSU& lfsu, const LocalPattern& p)
      {
        for (size_t k=0; k<p.size(); ++k)
          local_assembler.add_entry(*pattern,
                                    lfsv.globalIndex(p[k].i()),
                                    lfsu.globalIndex(p[k].j())
                                    );
      }


      //! Called when the local function space is about to be rebound or
      //! discarded
      //! @{
      template<typename EG>
      void onUnbindLFSUV(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        add_pattern(lfsv,lfsu,localpattern);
        localpattern.clear();
      }

      template<typename IG>
      void onUnbindLFSUVOutside(const IG& ig,
                                const LFSU& lfsu_s, const LFSV& lfsv_s,
                                const LFSU& lfsu_n, const LFSV& lfsv_n)
      {
        add_pattern(lfsv_s,lfsu_n,localpattern_sn);
        localpattern_sn.clear();
        add_pattern(lfsv_n,lfsu_s,localpattern_ns);
        localpattern_ns.clear();
      }

      //! @}

      //! Assembling methods
      //! @{

      template<typename EG>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternVolume>::
          pattern_volume(lop,lfsu,lfsv,localpattern);
      }

      template<typename IG>
      void assembleUVSkeleton(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s,
                              const LFSU & lfsu_n, const LFSV & lfsv_n)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternSkeleton>::
          pattern_skeleton(lop,lfsu_s,lfsv_s,lfsu_n,lfsv_n,
                           localpattern_sn, localpattern_ns);
      }

      template<typename IG>
      void assembleUVBoundary(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternBoundary>::
          pattern_boundary(lop,lfsu_s,lfsv_s,localpattern);
      }

      template<typename IG>
      static void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSU & lfsu_s, const LFSV & lfsv_s,
                                             const LFSU & lfsu_n, const LFSV & lfsv_n,
                                             const LFSU & lfsu_coupling, const LFSV & lfsv_coupling)
      {DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");}

      template<typename IG>
      static void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSV & lfsv_s,
                                            const LFSV & lfsv_n,
                                            const LFSV & lfsv_coupling)
      {DUNE_THROW(Dune::NotImplemented,"Assembling of coupling spaces is not implemented for ");}

      template<typename EG>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternVolumePostSkeleton>::
          pattern_volume_post_skeleton(lop,lfsu,lfsv,localpattern);
      }

      //! @}


    private:

      //! Reference to the wrapping local assembler object
      const LocalAssembler & local_assembler;

      //! Reference to the local operator
      const LOP & lop;

      //! Default value indicating an invalid solution pointer
      Pattern * const invalid_pattern;

      //! Pointer to the current matrix pattern container
      Pattern * pattern;

      //! Local pattern object used in assembling
      LocalPattern localpattern;
      LocalPattern localpattern_sn, localpattern_ns;
    }; // End of class DefaultLocalPatternAssemblerEngine

  };
};
#endif
