#ifndef DUNE_PDELAB_DEFAULT_PATTERNENGINE_HH
#define DUNE_PDELAB_DEFAULT_PATTERNENGINE_HH
namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for DUNE grids which creates
       the matrix pattern

       \tparam LA The local assembler

    */
    template<typename LA>
    class DefaultLocalPatternAssemblerEngine
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
      typedef typename LA::Pattern Pattern;
  
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
      const LocalAssembler & localAssembler(){ return local_assembler; }

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setPattern(Pattern & pattern_){
        pattern = &pattern_;
      }

      //! Query methods for the global grid assembler
      //! @{
      bool requireSkeleton() const 
      { return local_assembler.doPatternSkeleton(); }
      bool requireSkeletonTwoSided() const
      { return false; }
      bool requireUVVolume() const
      { return local_assembler.doPatternVolume(); }
      bool requireVVolume() const
      { return false; }
      bool requireUVSkeleton() const
      { return local_assembler.doPatternSkeleton(); }
      bool requireVSkeleton() const
      { return false; }
      bool requireUVBoundary() const
      { return local_assembler.doPatternBoundary(); }
      bool requireVBoundary() const
      { return false; }
      bool requireUVEnrichedCoupling() const
      { return false; }
      bool requireVEnrichedCoupling() const
      { return false; }
      bool requireUVVolumePostSkeleton() const
      { return local_assembler.doPatternVolumePostSkeleton(); }
      bool requireVVolumePostSkeleton() const
      { return false; }
      //! @}


      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG>
      void onBindLFSUV(const EG & eg, const LFSU & lfsu, const LFSV & lfsv){
        // Reset link container
        localpattern = LocalPattern();
      }

      template<typename EG>
      void onBindLFSV(const EG & eg, const LFSV & lfsv){}

      template<typename IG>
      void onBindLFSUVInside(const IG & ig, const LFSU & lfsu, const LFSV & lfsv){}

      template<typename IG>
      void onBindLFSUVOutside(const IG & ig, const LFSU & lfsun, const LFSV & lfsvn){
        // Reset link container
        localpattern_sn = LocalPattern();
        localpattern_ns = LocalPattern();
      }

      template<typename IG>
      void onBindLFSVInside(const IG & ig, const LFSV & lfsv){}

      template<typename IG>
      void onBindLFSVOutside(const IG & ig, const LFSV & lfsvn){}
      //! @}


      //! Called when the local function space is about to be rebound or
      //! discarded 
      //! @{
      template<typename EG>
      void onUnbindLFSUV(const EG & eg, const LFSU & lfsu, const LFSV & lfsv){
      
        for (size_t k=0; k<localpattern.size(); ++k)
          local_assembler.add_entry(*pattern,
                                    lfsv.globalIndex(localpattern[k].i()),
                                    lfsu.globalIndex(localpattern[k].j())
                                    );
      }

      template<typename EG>
      void onUnbindLFSV(const EG & eg, const LFSV & lfsv){ }

      template<typename IG>
      void onUnbindLFSUVInside(const IG & ig, const LFSU & lfsu, const LFSV & lfsv){}

      template<typename IG>
      void onUnbindLFSUVOutside(const IG & ig, const LFSU & lfsun, const LFSV & lfsvn){}

      template<typename IG>
      void onUnbindLFSVInside(const IG & ig, const LFSV & lfsv){

      }

      template<typename IG>
      void onUnbindLFSVOutside(const IG & ig, const LFSV & lfsvn){

      }
      //! @}

      //! Methods for loading of the local function's coefficients
      //! @{
      void loadCoefficientsLFSUInside(const LFSU & lfsu_s){}
      void loadCoefficientsLFSUOutside(const LFSU & lfsu_n){}
      void loadCoefficientsLFSUCoupling(const LFSU & lfsu_c){}
      //! @}

      //! Notifier functions, called immediately before and after assembling
      //! @{
      void preAssembly(){}
      void postAssembly(){}
      //! @}

      //! Assembling methods
      //! @{
      template<typename EG>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternVolume>::
          pattern_volume(lop,lfsu,lfsv,localpattern);
      }

      template<typename EG>
      void assembleVVolume(const EG & eg, const LFSV & lfsv){}

      template<typename IG>
      void assembleUVSkeleton(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s,
                              const LFSU & lfsu_n, const LFSV & lfsv_n)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternSkeleton>::
          pattern_skeleton(lop,lfsu_s,lfsv_s,lfsu_n,lfsv_n,
                           localpattern_sn, localpattern_ns);

        for (size_t k=0; k<localpattern_sn.size(); ++k)
          local_assembler.add_entry(*pattern,
                                    lfsv_s.globalIndex(localpattern_sn[k].i()),
                                    lfsu_n.globalIndex(localpattern_sn[k].j())
                                    );

        for (size_t k=0; k<localpattern_ns.size(); ++k)
          local_assembler.add_entry(*pattern,
                                    lfsv_n.globalIndex(localpattern_ns[k].i()),
                                    lfsu_s.globalIndex(localpattern_ns[k].j())
                                    );
      }

      template<typename IG>
      void assembleVSkeleton(const IG & ig, const LFSV & lfsv_s, const LFSV & lfsv_n){}

      template<typename IG>
      void assembleUVBoundary(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s)
      {
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doPatternBoundary>::
          pattern_boundary(lop,lfsu_s,lfsv_s,localpattern);
      }

      template<typename IG>
      void assembleVBoundary(const IG & ig, const LFSV & lfsv_s){}

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

      template<typename EG>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv){}

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
