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
    template<typename GO, typename LA0, typename LA1>
    class OneStepLocalAssembler
      : public Dune::PDELab::LocalAssemblerBase<
      typename GO::Traits::MatrixBackend,
      typename GO::Traits::TrialGridFunctionSpaceConstraints,
      typename GO::Traits::TestGridFunctionSpaceConstraints>
    {
    public:

      //! The types of the local assemblers of order one and zero
      typedef LA0 LocalAssemblerDT0;
      typedef LA1 LocalAssemblerDT1;

      typedef Dune::PDELab::LocalAssemblerTraits<GO> Traits;

      //! The base class
      typedef Dune::PDELab::LocalAssemblerBase
      <typename GO::Traits::MatrixBackend,
       typename GO::Traits::TrialGridFunctionSpaceConstraints,
       typename GO::Traits::TestGridFunctionSpaceConstraints> Base;

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
        dune_static_assert((is_same<typename LA0::Traits::MatrixBackend::Pattern,
                            typename LA1::Traits::MatrixBackend::Pattern>::value),
                           "Received two local assemblers which are non-compatible "
                           "due to different matrix pattern types");
        dune_static_assert((is_same<typename LA0::Traits::Jacobian,
                            typename LA1::Traits::Jacobian>::value),
                           "Received two local assemblers which are non-compatible "
                           "due to different jacobian types");
        dune_static_assert((is_same<typename LA0::Traits::Solution,
                            typename LA1::Traits::Solution>::value),
                           "Received two local assemblers which are non-compatible "
                           "due to different solution vector types");
        dune_static_assert((is_same<typename LA0::Traits::Residual,
                            typename LA1::Traits::Residual>::value),
                           "Received two local assemblers which are non-compatible "
                           "due to different residual vector types");
      }

      //! The local operators type for real numbers e.g. time
      typedef typename Traits::RangeField Real;

      // //! The residual representation type
      // typedef typename LA0::GridOperator::Traits::Range Residual;

      // //! The solution representation type
      // typedef typename LA0::GridOperator::Traits::Domain Solution;

      // //! The jacobian representation type
      // typedef typename LA0::GridOperator::Traits::Jacobian Jacobian;

      // //! The matrix pattern representation type
      // typedef typename LA0::GridOperator::Traits::MatrixBackend::Pattern Pattern;

      //! The type of the one step parameter object
      typedef Dune::PDELab::TimeSteppingParameterInterface<Real> OneStepParameters;

      //! Constructor with empty constraints
      OneStepLocalAssembler (LA0 & la0_, LA1 & la1_, typename Traits::Residual & const_residual_)
        : Base(la0_.trialConstraints(),la0_.testConstraints()),
          la0(la0_), la1(la1_), const_residual(const_residual_),
          time(0.0), use_mass_dt(false), stage(0),
          pattern_engine(*this), prestage_engine(*this), residual_engine(*this), jacobian_engine(*this)
      { static_checks(); }

      //! Notifies the local assembler about the current time of
      //! assembling. Should be called before assembling if the local
      //! operator has time dependencies.
      void preStep(Real time_, Real dt_, int stages_){
        time = time_;
        dt = dt_;
        if(use_mass_dt){
          dt_factor0 = 1.0;
          dt_factor1 = 1.0 / dt;
        }
        else{
          dt_factor0 = dt;
          dt_factor1 = 1.0;
        }

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

      //! Determines whether the time step size is multiplied to the
      //! mass term (first order time derivative) or the elliptic term
      //! (zero-th order time derivative).
      void divideMassTermByDeltaT(bool v){
        use_mass_dt = v;
      }

      //! Access time at given stage
      Real timeAtStage(int stage_){
        return time+osp_method->d(stage_)*dt;
      }

      //! Access time at given stage
      Real timeAtStage(){
        return time+osp_method->d(stage)*dt;
      }

      //! Access methods which provid "ready to use" engines
      //! @{

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalPatternAssemblerEngine & localPatternAssemblerEngine
      (typename Traits::MatrixPattern & p)
      {
        pattern_engine.setPattern(p);
        return pattern_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalPreStageAssemblerEngine & localPreStageAssemblerEngine
      (const std::vector<typename Traits::Solution*> & x)
      {
        prestage_engine.setSolutions(x);
        prestage_engine.setConstResidual(const_residual);
        return prestage_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalResidualAssemblerEngine & localResidualAssemblerEngine
      (typename Traits::Residual & r, const typename Traits::Solution & x)
      {
        residual_engine.setSolution(x);
        residual_engine.setResidual(r);
        return residual_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalJacobianAssemblerEngine & localJacobianAssemblerEngine
      (typename Traits::Jacobian & a, const typename Traits::Solution & x)
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
      typename Traits::Residual & const_residual;

      //! The current time of assembling
      Real time;

      //! The time step size
      Real dt;

      /** The time step factors for assembling. Depending on the value
       of \a use_mass_dt, it will hold:

       dt_factor0 = dt and dt_factor1 = 1.0

       or

       dt_factor0 = 1.0 and dt_factor1 = 1.0 / dt
      */
      Real dt_factor0, dt_factor1;

      //! Determines, whether the time step size will be multiplied
      //! with the time derivative term of first of zero-th order.
      bool use_mass_dt;

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
