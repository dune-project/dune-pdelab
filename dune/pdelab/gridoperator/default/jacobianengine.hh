#ifndef DUNE_PDELAB_DEFAULT_JACOBIANENGINE_HH
#define DUNE_PDELAB_DEFAULT_JACOBIANENGINE_HH

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for DUNE grids which
       assembles the jacobian matrix

       \tparam LA The local assembler

    */
    template<typename LA>
    class DefaultLocalJacobianAssemblerEngine
    {
    public:
      //! The type of the wrapping local assembler
      typedef LA LocalAssembler;

      //! The type of the local operator
      typedef typename LA::LocalOperator LOP;

      //! The local function spaces
      typedef typename LA::LFSU LFSU;
      typedef typename LA::LFSV LFSV;

      //! The type of the jacobian matrix
      typedef typename LA::Jacobian Jacobian;
      typedef typename Jacobian::ElementType JacobianElement;

      //! The type of the solution vector
      typedef typename LA::Solution Solution;
      typedef typename Solution::ElementType SolutionElement;
  
      /**
         \brief Constructor 

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      DefaultLocalJacobianAssemblerEngine(const LocalAssembler & local_assembler_)
        : local_assembler(local_assembler_), lop(local_assembler_.lop), 
          invalid_jacobian(static_cast<Jacobian*>(0)),
          invalid_solution(static_cast<Solution*>(0)),
          jacobian(invalid_jacobian), 
          solution(invalid_solution)
      {}
  
      //! Query methods for the global grid assembler
      //! @{
      bool requireSkeleton() const 
      { return local_assembler.doAlphaSkeleton(); }
      bool requireSkeletonTwoSided() const
      { return local_assembler.doSkeletonTwoSided(); }
      bool requireUVVolume() const
      { return local_assembler.doAlphaVolume(); }
      bool requireVVolume() const
      { return false; }
      bool requireUVSkeleton() const
      { return local_assembler.doAlphaSkeleton(); }
      bool requireVSkeleton() const
      { return false; }
      bool requireUVBoundary() const
      { return local_assembler.doAlphaBoundary(); }
      bool requireVBoundary() const
      { return false; }
      bool requireUVEnrichedCoupling() const
      { return false; }
      bool requireVEnrichedCoupling() const
      { return false; }
      bool requireUVVolumePostSkeleton() const
      { return local_assembler.doAlphaVolumePostSkeleton(); }
      bool requireVVolumePostSkeleton() const
      { return false; }
      //! @}

      //! Public access to the wrapping local assembler
      const LocalAssembler & localAssembler(){ return local_assembler; }

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setJacobian(Jacobian & jacobian_){
        jacobian = &jacobian_;
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
        al.assign(lfsv.size() ,lfsu.size(),0.0);
        alc.assign(lfsv.size() ,lfsu.size(),0.0);
      }

      template<typename EG>
      void onBindLFSV(const EG & eg, const LFSV & lfsv){

      }

      template<typename IG>
      void onBindLFSUVInside(const IG & ig, const LFSU & lfsu, const LFSV & lfsv){
      }

      template<typename IG>
      void onBindLFSUVOutside(const IG & ig, const LFSU & lfsun, const LFSV & lfsvn){
        xn.resize(lfsun.size());
      }

      template<typename IG>
      void onBindLFSVInside(const IG & ig, const LFSV & lfsv){

      }

      template<typename IG>
      void onBindLFSVOutside(const IG & ig, const LFSV & lfsvn){

      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded 
      //! @{
      template<typename EG>
      void onUnbindLFSUV(const EG & eg, const LFSU & lfsu, const LFSV & lfsv){
        local_assembler.etadd(lfsv,lfsu,al,*jacobian);
      }

      template<typename EG>
      void onUnbindLFSV(const EG & eg, const LFSV & lfsv){}

      template<typename IG>
      void onUnbindLFSUVInside(const IG & ig, const LFSU & lfsu, const LFSV & lfsv){}

      template<typename IG>
      void onUnbindLFSUVOutside(const IG & ig, const LFSU & lfsun, const LFSV & lfsvn){}

      template<typename IG>
      void onUnbindLFSVInside(const IG & ig, const LFSV & lfsv){}

      template<typename IG>
      void onUnbindLFSVOutside(const IG & ig, const LFSV & lfsvn){}
      //! @}

      //! Methods for loading of the local function's coefficients
      //! @{
      void loadCoefficientsLFSUInside(const LFSU & lfsu){
        lfsu.vread(*solution,xl);
      }
      void loadCoefficientsLFSUOutside(const LFSU & lfsun){
        lfsun.vread(*solution,xn);
      }
      void loadCoefficientsLFSUCoupling(const LFSU & lfsu_c)
      {DUNE_THROW(Dune::NotImplemented,"No coupling lfsu available for ");}
      //! @}

#ifndef DOXYGEN
      // Dummy implementation for bind notifiers on coupling function
      // spaces
      void onBindLFSUVCoupling(const LFSU & lfsu, const LFSV & lfsv){}
      void onUnbindLFSUVCoupling(const LFSU & lfsu, const LFSV & lfsv){}
      void onBindLFSVCoupling(const LFSV & lfsv){}
      void onUnbindLFSVCoupling(const LFSV & lfsv){}
#endif

      //! Notifier functions, called immediately before and after assembling
      //! @{
      void preAssembly(){}
      void postAssembly(){ 
        if(local_assembler.doConstraintsPostProcessing){
          typedef typename LocalAssembler::Base::Traits::TestConstraintsType::const_iterator 
            global_row_iterator;       
          for (global_row_iterator cit=(local_assembler.pconstraintsv)->begin(); 
               cit!=(local_assembler.pconstraintsv)->end(); ++cit)
            local_assembler.set_trivial_row(cit->first,cit->second,*jacobian);
        }
      }
      //! @}

      //! Assembling methods
      //! @{
      template<typename EG>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        assign(alc,0.0);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolume>::
          jacobian_volume(lop,eg,lfsu,xl,lfsv,alc);
        update(al,local_assembler.weight,alc);
      }

      template<typename EG>
      void assembleVVolume(const EG & eg, const LFSV & lfsv)
      {}

      template<typename IG>
      void assembleUVSkeleton(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s,
                              const LFSU & lfsu_n, const LFSV & lfsv_n)
      {
        assign(alc,0.0);

        al_sn.assign(lfsv_s.size() ,lfsu_n.size(),0.0);
        al_ns.assign(lfsv_n.size(),lfsu_s.size() ,0.0);
        al_nn.assign(lfsv_n.size(),lfsu_n.size(),0.0);

        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
          jacobian_skeleton(lop,ig,lfsu_s,xl,lfsv_s,lfsu_n,xn,lfsv_n,alc,al_sn,al_ns,al_nn);

        update(al,local_assembler.weight,alc);

        al_sn *= local_assembler.weight;
        al_ns *= local_assembler.weight;
        al_nn *= local_assembler.weight;

        local_assembler.etadd(lfsv_s,lfsu_n,al_sn,*jacobian);
        local_assembler.etadd(lfsv_n,lfsu_s,al_ns,*jacobian);
        local_assembler.etadd(lfsv_n,lfsu_n,al_nn,*jacobian);
      }

      template<typename IG>
      void assembleVSkeleton(const IG & ig, const LFSV & lfsv_s, const LFSV & lfsv_n){}

      template<typename IG>
      void assembleUVBoundary(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s)
      {
        assign(alc,0.0);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          jacobian_boundary(lop,ig,lfsu_s,xl,lfsv_s,alc);                              
        update(al,local_assembler.weight,alc);
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
        assign(alc,0.0);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolumePostSkeleton>::
          jacobian_volume_post_skeleton(lop,eg,lfsu,xl,lfsv,alc);
        update(al,local_assembler.weight,alc);
      }

      template<typename EG>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv){}

      //! @}

    private:
      //! Reference to the wrapping local assembler object which
      //! constructed this engine
      const LocalAssembler & local_assembler;

      //! Reference to the local operator
      const LOP & lop;

      //! Default value indicating an invalid residual pointer
      Jacobian * const invalid_jacobian;

      //! Default value indicating an invalid solution pointer
      Solution * const invalid_solution;

      //! Pointer to the current residual vector in which to assemble
      Jacobian * jacobian;

      //! Pointer to the current residual vector in which to assemble
      const Solution * solution;

      //! The local vectors and matrices as required for assembling
      //! @{
      typedef Dune::PDELab::TrialSpaceTag LocalTrialSpaceTag;
      typedef Dune::PDELab::TestSpaceTag LocalTestSpaceTag;
      Dune::PDELab::LocalVector<SolutionElement, LocalTrialSpaceTag> xl;
      Dune::PDELab::LocalVector<SolutionElement, LocalTrialSpaceTag> xn;
      Dune::PDELab::LocalMatrix<JacobianElement> al;
      Dune::PDELab::LocalMatrix<JacobianElement> al_sn;
      Dune::PDELab::LocalMatrix<JacobianElement> al_ns;
      Dune::PDELab::LocalMatrix<JacobianElement> al_nn;

      Dune::PDELab::LocalMatrix<JacobianElement> alc;
      //! @}

      template<typename VC, typename VF>
      void update(VC & A, const VF & a, const VC & B){
        assert(A.nrows() == B.nrows());
        assert(A.ncols() == B.ncols());
        const unsigned int R = A.nrows();
        const unsigned int C = A.ncols();
        for(unsigned int r=0; r<R; ++r)
          for(unsigned int c=0; c<C; ++c)
            A(r,c) += a * B(r,c);
      }

      template<typename VC, typename VF>
      void assign(VC & A, const VF & a){
        const unsigned int R = A.nrows();
        const unsigned int C = A.ncols();
        for(unsigned int r=0; r<R; ++r)
          for(unsigned int c=0; c<C; ++c)
            A(r,c) = a;
      }

    
    }; // End of class DefaultLocalJacobianAssemblerEngine

  };
};
#endif
