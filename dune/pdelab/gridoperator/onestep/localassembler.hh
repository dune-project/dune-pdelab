#ifndef DUNE_PDELAB_ONESTEP_LOCAL_ASSEMBLER_HH
#define DUNE_PDELAB_ONESTEP_LOCAL_ASSEMBLER_HH

#include <dune/pdelab/gridoperator/onestep/residualengine.hh>
#include <dune/pdelab/gridoperator/onestep/patternengine.hh>
#include <dune/pdelab/gridoperator/onestep/jacobianengine.hh>
#include <dune/pdelab/gridoperator/onestep/prestageengine.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/common/typetree.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler for one step methods

       \tparam LA0 The local assembler for the temporal derivative term of order zero
       \tparam LA1 The local assembler for the temporal derivative term of order one
    */
    template<typename LA0, typename LA1>
    class OneStepLocalAssembler 
      : public Dune::PDELab::LocalAssemblerBase< 
      typename LA0::Traits::MatrixBackend, 
      typename LA0::Traits::TrialConstraintsType, 
      typename LA0::Traits::TestConstraintsType>
    {
    public:

      //! The traits class from the 
      typedef typename LA0::Traits Traits;

      //! The types of the local assemblers of order one and zero
      typedef LA0 LocalAssemblerDT0;
      typedef LA1 LocalAssemblerDT1;

      //! The base class
      typedef Dune::PDELab::LocalAssemblerBase
      < typename LA0::Traits::MatrixBackend, 
        typename LA0::Traits::TrialConstraintsType, 
        typename LA0::Traits::TestConstraintsType> Base;

      //! The local assembler engines
      //! @{
      typedef OneStepLocalPatternAssemblerEngine<OneStepLocalAssembler> LocalPatternAssemblerEngine;
      typedef OneStepLocalPreStageAssemblerEngine<OneStepLocalAssembler> LocalPreStageAssemblerEngine;
      typedef OneStepLocalResidualAssemblerEngine<OneStepLocalAssembler> LocalResidualAssemblerEngine;
      typedef OneStepLocalJacobianAssemblerEngine<OneStepLocalAssembler> LocalJacobianAssemblerEngine;

      friend class OneStepLocalPatternAssemblerEngine<OneStepLocalAssembler>;
      friend class OneStepLocalPreStageAssemblerEngine<OneStepLocalAssembler>;
      friend class OneStepLocalResidualAssemblerEngine<OneStepLocalAssembler>;
      friend class OneStepLocalJacobianAssemblerEngine<OneStepLocalAssembler>;
      //! @}

      void static_checks(){
        dune_static_assert((is_same<typename LA0::Pattern,typename LA1::Pattern>::value),
                           "Received two local assemblers which are non-compatible "
                           "due to different matrix pattern types");
        dune_static_assert((is_same<typename LA0::Jacobian,typename LA1::Jacobian>::value),
                           "Received two local assemblers which are non-compatible "
                           "due to different jacobian types");
        dune_static_assert((is_same<typename LA0::Solution,typename LA1::Solution>::value),
                           "Received two local assemblers which are non-compatible "
                           "due to different solution vector types");
        dune_static_assert((is_same<typename LA0::Residual,typename LA1::Residual>::value),
                           "Received two local assemblers which are non-compatible "
                           "due to different residual vector types");
        dune_static_assert((is_same<typename LA0::Real,typename LA1::Real>::value),
                           "Received two local assemblers which are non-compatible "
                           "due to different real number types");
      }

      //! The local operators type for real numbers e.g. time
      typedef typename LA0::Real Real;

      //! The residual representation type
      typedef typename LA0::Residual Residual;

      //! The solution representation type
      typedef typename LA0::Solution Solution;

      //! The jacobian representation type
      typedef typename LA0::Jacobian Jacobian;

      //! The matrix pattern representation type
      typedef typename LA0::Pattern Pattern;

      //! The type of the one step parameter object
      typedef Dune::PDELab::TimeSteppingParameterInterface<Real> OneStepParameters;

      //! Constructor with empty constraints
      OneStepLocalAssembler (LA0 & la0_, LA1 & la1_, Residual & const_residual_) 
        : la0(la0_), la1(la1_), const_residual(const_residual_), 
          time(0.0), stage(0),
          pattern_engine(*this), prestage_engine(*this), residual_engine(*this), jacobian_engine(*this)
      { static_checks(); }

      //! Notifies the local assembler about the current time of
      //! assembling. Should be called before assembling if the local
      //! operator has time dependencies.
      void preStep(Real time_, Real dt_, int stages_){
        time = time_;
        dt = dt_;
        la0.preStep(time_,dt_, stages_);
        la1.preStep(time_,dt_, stages_);
      }

      //! Set the one step method parameters
      void setMethod(const OneStepParameters & method_){
        osp_method = & method_;
      }

      //! Set the current stage of the one step scheme
      void setStage(int stage_){
        stage = stage_;
      }

      //! Access time at given stage
      Real timeAtStage(int stage_){
        return time+osp_method->d(stage_)*dt;
      }

      //! Access methods which provid "ready to use" engines
      //! @{

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalPatternAssemblerEngine & localPatternAssemblerEngine
      (Pattern & p)
      {
        pattern_engine.setPattern(p);
        return pattern_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalPreStageAssemblerEngine & localPreStageAssemblerEngine
      (const std::vector<Solution*> & x)
      {
        prestage_engine.setSolutions(x);
        prestage_engine.setConstResidual(const_residual);
        return prestage_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalResidualAssemblerEngine & localResidualAssemblerEngine
      (Residual & r, const Solution & x)
      {
        residual_engine.setSolution(x);
        residual_engine.setResidual(r);
        return residual_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalJacobianAssemblerEngine & localJacobianAssemblerEngine
      (Jacobian & a, const Solution & x)
      {
        jacobian_engine.setSolution(x);
        jacobian_engine.setJacobian(a);
        return jacobian_engine;
      }

      //! @}

    private:

      //! The local assemblers for the temporal derivative of order
      //! one and zero
      //! @{
      LA0 & la0;
      LA1 & la1;
      //! @}

      //! The one step parameter object containing the generalized
      //! butcher tableau parameters
      const OneStepParameters * osp_method;

      //! The constant part of the residual
      Residual & const_residual;
      
      //! The current time of assembling
      Real time, dt;

      //! The current stage of the one step scheme
      int stage;

      //! The engine member objects
      //! @{
      LocalPatternAssemblerEngine  pattern_engine;
      LocalPreStageAssemblerEngine prestage_engine;
      LocalResidualAssemblerEngine residual_engine;
      LocalJacobianAssemblerEngine jacobian_engine;
      //! @}
    };

  };
};
#endif
