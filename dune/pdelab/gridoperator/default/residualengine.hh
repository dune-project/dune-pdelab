#ifndef DUNE_PDELAB_DEFAULT_RESIDUALENGINE_HH
#define DUNE_PDELAB_DEFAULT_RESIDUALENGINE_HH

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for DUNE grids which
       assembles the residual vector

       \tparam LA The local assembler

    */
    template<typename LA>
    class DefaultLocalResidualAssemblerEngine
    {
    public:
      //! The type of the wrapping local assembler
      typedef LA LocalAssembler;

      //! The type of the local operator
      typedef typename LA::LocalOperator LOP;

      //! The type of the residual vector
      typedef typename LA::Residual Residual;
      typedef typename Residual::ElementType ResidualElement;

      //! The type of the solution vector
      typedef typename LA::Solution Solution;
      typedef typename Solution::ElementType SolutionElement;

      //! The local function spaces
      typedef typename LA::LFSU LFSU;
      typedef typename LA::LFSV LFSV;

      /**
         \brief Constructor 

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      DefaultLocalResidualAssemblerEngine(const LocalAssembler & local_assembler_)
        : local_assembler(local_assembler_), lop(local_assembler_.lop), 
          invalid_residual(static_cast<Residual*>(0)), invalid_solution(static_cast<Solution*>(0)),
          residual(invalid_residual), 
          solution(invalid_solution)
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
      bool requireUVEnrichedCoupling() const
      { return false; }
      bool requireVEnrichedCoupling() const
      { return false; }
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
        residual = &residual_;
      }

      //! Set current solution vector. Should be called prior to
      //! assembling.
      void setSolution(const Solution & solution_){
        solution = &solution_;
      }
    
      //! Called immediately after binding of local function space in
      //! global assembler.
      //! @{
      template<typename EG>
      void onBindLFSUV(const EG & eg, const LFSU & lfsu, const LFSV & lfsv){
        xl.resize(lfsu.size());
      }

      template<typename EG>
      void onBindLFSV(const EG & eg, const LFSV & lfsv){
        rl.assign(lfsv.size(),0.0);
        rlc.assign(lfsv.size(),0.0);
      }

      template<typename IG>
      void onBindLFSUVInside(const IG & ig, const LFSU & lfsu, const LFSV & lfsv){
        xl.resize(lfsu.size());
      }

      template<typename IG>
      void onBindLFSUVOutside(const IG & ig, const LFSU & lfsun, const LFSV & lfsvn){
        xn.resize(lfsun.size());
      }

      template<typename IG>
      void onBindLFSVInside(const IG & ig, const LFSV & lfsv){
        rl.assign(lfsv.size(),0.0);
        rlc.assign(lfsv.size(),0.0);
      }

      template<typename IG>
      void onBindLFSVOutside(const IG & ig, const LFSV & lfsvn){
        rn.assign(lfsvn.size(),0.0);
        rnc.assign(lfsvn.size(),0.0);
      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded 
      //! @{
      template<typename EG>
      void onUnbindLFSUV(const EG & eg, const LFSU & lfsu, const LFSV & lfsv){}

      template<typename EG>
      void onUnbindLFSV(const EG & eg, const LFSV & lfsv){
        lfsv.vadd(rl,*residual);
      }

      template<typename IG>
      void onUnbindLFSUVInside(const IG & ig, const LFSU & lfsu, const LFSV & lfsv){}

      template<typename IG>
      void onUnbindLFSUVOutside(const IG & ig, const LFSU & lfsun, const LFSV & lfsvn){}

      template<typename IG>
      void onUnbindLFSVInside(const IG & ig, const LFSV & lfsv){
        lfsv.vadd(rl,*residual);
      }

      template<typename IG>
      void onUnbindLFSVOutside(const IG & ig, const LFSV & lfsvn){
        lfsvn.vadd(rn,*residual);
      }
      //! @}

      //! Methods for loading of the local function's coefficients
      //! @{
      void loadCoefficientsLFSUInside(const LFSU & lfsu_s){
        lfsu_s.vread(*solution,xl);
      }
      void loadCoefficientsLFSUOutside(const LFSU & lfsu_n){
        lfsu_n.vread(*solution,xn);
      }
      void loadCoefficientsLFSUCoupling(const LFSU & lfsu_c)
      {DUNE_THROW(Dune::NotImplemented,"No coupling lfsu available for ");}
      //! @}


#ifndef DOXYGEN
      // Dummy implementation for bind notifiers on coupling function
      // spaces
      void onBindLFSUCoupling(const LFSU & lfsu){}
      void onUnbindLFSUCoupling(const LFSU & lfsu){}
      void onBindLFSVCoupling(const LFSV & lfsv){}
      void onUnbindLFSVCoupling(const LFSV & lfsv){}
#endif

      //! Notifier functions, called immediately before and after assembling
      //! @{
      void preAssembly(){ }
      void postAssembly(){ 
        if(local_assembler.doConstraintsPostProcessing){
          Dune::PDELab::constrain_residual(*(local_assembler.pconstraintsv),*residual); 
        }
      }
      //! @}

      //! Assembling methods
      //! @{

      template<typename EG>
      bool assembleCell(const EG & eg)
      {
        return LocalAssembler::isNonOverlapping && eg.entity().partitionType() != Dune::InteriorEntity;
      }

      template<typename EG>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        assign(rlc,0);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolume>::
          alpha_volume(lop,eg,lfsu,xl,lfsv,rlc);
        axpy(rl,local_assembler.weight,rlc);
      }

      template<typename EG>
      void assembleVVolume(const EG & eg, const LFSV & lfsv)
      {
        assign(rlc,0);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaVolume>::
          lambda_volume(lop,eg,lfsv,rlc);
        axpy(rl,local_assembler.weight,rlc);
      }

      template<typename IG>
      void assembleUVSkeleton(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s,
                              const LFSU & lfsu_n, const LFSV & lfsv_n)
      {
        assign(rlc,0);
        assign(rnc,0);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
          alpha_skeleton(lop,ig,
                         lfsu_s,xl,lfsv_s,
                         lfsu_n,xn,lfsv_n,
                         rlc,rnc);

        axpy(rl,local_assembler.weight,rlc);
        axpy(rn,local_assembler.weight,rnc);
      }

      template<typename IG>
      void assembleVSkeleton(const IG & ig, const LFSV & lfsv_s, const LFSV & lfsv_n)
      {
        assign(rlc,0);
        assign(rnc,0);

        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaSkeleton>::
          lambda_skeleton(lop, ig, lfsv_s, lfsv_n, rlc, rnc);

        axpy(rl,local_assembler.weight,rlc);
        axpy(rn,local_assembler.weight,rnc);
      }

      template<typename IG>
      void assembleUVBoundary(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s)
      {
        assign(rlc,0);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          alpha_boundary(lop,ig,lfsu_s,xl,lfsv_s,rlc);
        axpy(rl,local_assembler.weight,rlc);
      }

      template<typename IG>
      void assembleVBoundary(const IG & ig, const LFSV & lfsv_s)
      {
        assign(rlc,0);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaBoundary>::
          lambda_boundary(lop,ig,lfsv_s,rlc);
        axpy(rl,local_assembler.weight,rlc);
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
        assign(rlc,0);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolumePostSkeleton>::
          alpha_volume_post_skeleton(lop,eg,lfsu,xl,lfsv,rlc);
        axpy(rl,local_assembler.weight,rlc);
      }

      template<typename EG>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv)
      {
        assign(rlc,0);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doLambdaVolumePostSkeleton>::
          lambda_volume_post_skeleton(lop,eg,lfsv,rlc);
        axpy(rl,local_assembler.weight,rlc);
      }

      //! @}

    private:
      //! Reference to the wrapping local assembler object which
      //! constructed this engine
      const LocalAssembler & local_assembler;

      //! Reference to the local operator
      const LOP & lop;

      //! Default value indicating an invalid residual pointer
      Residual * const invalid_residual;

      //! Default value indicating an invalid solution pointer
      Solution * const invalid_solution;

      //! Pointer to the current residual vector in which to assemble
      Residual * residual;

      //! Pointer to the current residual vector in which to assemble
      const Solution * solution;

      //! The local vectors and matrices as required for assembling
      //! @{
      typedef Dune::PDELab::TrialSpaceTag LocalTrialSpaceTag;
      typedef Dune::PDELab::TestSpaceTag LocalTestSpaceTag;
      //! Inside local coefficients
      Dune::PDELab::LocalVector<SolutionElement, LocalTrialSpaceTag> xl;
      //! Outside local coefficients
      Dune::PDELab::LocalVector<SolutionElement, LocalTrialSpaceTag> xn;
      //! Inside local residual
      Dune::PDELab::LocalVector<ResidualElement, LocalTestSpaceTag> rl;
      //! Outside local residual
      Dune::PDELab::LocalVector<ResidualElement, LocalTestSpaceTag> rn;
      //! Inside local residual copy
      Dune::PDELab::LocalVector<ResidualElement, LocalTestSpaceTag> rlc;
      //! Outside local residual copy
      Dune::PDELab::LocalVector<ResidualElement, LocalTestSpaceTag> rnc;
      //! @}

      template<typename VC, typename VF>
      void axpy(VC & y, const VF & a, const VC & x){
        assert(y.size() == x.size());
        for(unsigned int i=0; i<y.size(); ++i)
          y[i] += a * x[i];
      }

      template<typename VC, typename VF>
      void assign(VC & y, const VF & a){
        for(unsigned int i=0; i<y.size(); ++i)
          y[i] = a;
      }

    }; // End of class DefaultLocalResidualAssemblerEngine

  };
};
#endif
