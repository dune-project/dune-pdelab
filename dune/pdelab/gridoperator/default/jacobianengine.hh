#ifndef DUNE_PDELAB_DEFAULT_JACOBIANENGINE_HH
#define DUNE_PDELAB_DEFAULT_JACOBIANENGINE_HH

#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for DUNE grids which
       assembles the jacobian matrix

       \tparam LA The local assembler

    */
    template<typename LA>
    class DefaultLocalJacobianAssemblerEngine
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
          solution(invalid_solution),
          al_view(al,1.0),
          al_sn_view(al_sn,1.0),
          al_ns_view(al_ns,1.0),
          al_nn_view(al_nn,1.0)
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
      }

      //! @}

      //! Called when the local function space is about to be rebound or
      //! discarded 
      //! @{
      template<typename EG>
      void onUnbindLFSUV(const EG & eg, const LFSU & lfsu, const LFSV & lfsv){
        local_assembler.etadd(lfsv,lfsu,al,*jacobian);
      }

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

      //! Notifier functions, called immediately before and after assembling
      //! @{
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
      bool assembleCell(const EG & eg)
      {
        return LocalAssembler::isNonOverlapping && eg.entity().partitionType() != Dune::InteriorEntity;
      }

      template<typename EG>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        al_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolume>::
          jacobian_volume(lop,eg,lfsu,xl,lfsv,al_view);
      }

      template<typename IG>
      void assembleUVSkeleton(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s,
                              const LFSU & lfsu_n, const LFSV & lfsv_n)
      {
        al_sn.assign(lfsv_s.size() ,lfsu_n.size(),0.0);
        al_ns.assign(lfsv_n.size(),lfsu_s.size() ,0.0);
        al_nn.assign(lfsv_n.size(),lfsu_n.size(),0.0);

        al_view.setWeight(local_assembler.weight);
        al_sn_view.setWeight(local_assembler.weight);
        al_ns_view.setWeight(local_assembler.weight);
        al_nn_view.setWeight(local_assembler.weight);

        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
          jacobian_skeleton(lop,ig,lfsu_s,xl,lfsv_s,lfsu_n,xn,lfsv_n,al_view,al_sn_view,al_ns_view,al_nn_view);

        local_assembler.etadd(lfsv_s,lfsu_n,al_sn,*jacobian);
        local_assembler.etadd(lfsv_n,lfsu_s,al_ns,*jacobian);
        local_assembler.etadd(lfsv_n,lfsu_n,al_nn,*jacobian);
      }

      template<typename IG>
      void assembleUVBoundary(const IG & ig, const LFSU & lfsu_s, const LFSV & lfsv_s)
      {
        al_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          jacobian_boundary(lop,ig,lfsu_s,xl,lfsv_s,al_view);
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
        al_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolumePostSkeleton>::
          jacobian_volume_post_skeleton(lop,eg,lfsu,xl,lfsv,al_view);
      }

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

      typedef Dune::PDELab::LocalVector<SolutionElement, LocalTrialSpaceTag> SolutionVector;
      typedef Dune::PDELab::LocalMatrix<JacobianElement> JacobianMatrix;

      SolutionVector xl;
      SolutionVector xn;

      JacobianMatrix al;
      JacobianMatrix al_sn;
      JacobianMatrix al_ns;
      JacobianMatrix al_nn;

      typename JacobianMatrix::WeightedAccumulationView al_view;
      typename JacobianMatrix::WeightedAccumulationView al_sn_view;
      typename JacobianMatrix::WeightedAccumulationView al_ns_view;
      typename JacobianMatrix::WeightedAccumulationView al_nn_view;

      //! @}

    
    }; // End of class DefaultLocalJacobianAssemblerEngine

  };
};
#endif
