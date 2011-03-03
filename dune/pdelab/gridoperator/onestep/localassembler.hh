#ifndef DUNE_UDG_LOCAL_ASSEMBLER_HH
#define DUNE_UDG_LOCAL_ASSEMBLER_HH

#include <dune/udg/pdelab/assembler/residualengine.hh>
#include <dune/udg/pdelab/assembler/patternengine.hh>
#include <dune/udg/pdelab/assembler/jacobianengine.hh>
#include <dune/pdelab/gridoperatorspace/assemblerutilities.hh>
#include <dune/pdelab/common/typetree.hh>

namespace Dune{
  namespace UDG{

    /**
       \brief The local assembler for one step methods

       \tparam LA0 The local assembler for the temporal derivative term of order zero
       \tparam LA1 The local assembler for the temporal derivative term of order one
    */
    template<typename LA0, typename LA1>
    class OneStepLocalAssembler 
      : public Dune::PDELab::LocalAssemblerBase<
      typename LA0::Traits::MatrixBackendType,
      typename LA0::Traits::TrialConstraintsType,
      typename LA0::Traits::TestConstraintsType  >
    {
    public:

      //! The traits class from the 
      typedef typename LA0::Traits Traits;

      //! The types of the local assemblers of order one and zero
      typedef LA0 LocalAssemblerDT0;
      typedef LA1 LocalAssemblerDT1;

      //! The local operators type for real numbers e.g. time
      typedef typename LA1::Real Real;

      //! The type of the one step parameter object
      Dune::PDELab::TimeSteppingParameterInterface<Real> OneStepParameters;

      //! The local assembler engines
      //! @{
      typedef OneStepLocalPatternAssemblerEngine<OneStepLocalAssembler> LocalPatternAssemblerEngine;
      typedef OneStepLocalResidualAssemblerEngine<OneStepLocalAssembler> LocalResidualAssemblerEngine;
      typedef OneStepLocalJacobianAssemblerEngine<OneStepLocalAssembler> LocalJacobianAssemblerEngine;

      friend class OneStepLocalPatternAssemblerEngine<OneStepLocalAssembler>;
      friend class OneStepLocalResidualAssemblerEngine<OneStepLocalAssembler>;
      friend class OneStepLocalJacobianAssemblerEngine<OneStepLocalAssembler>;
      //! @}

      void static_checks(){
        Dune::dune_static_assert(Dune::is_same<LA0::Pattern,LA1::Pattern>::true,
                                 "Received two local assemblers which are non-compatible "
                                 "due to different matrix pattern types");
        Dune::dune_static_assert(Dune::is_same<LA0::Jacobian,LA1::Jacobian>::true,
                                 "Received two local assemblers which are non-compatible "
                                 "due to different jacobian types");
        Dune::dune_static_assert(Dune::is_same<LA0::Solution,LA1::Solution>::true,
                                 "Received two local assemblers which are non-compatible "
                                 "due to different solution vector types");
        Dune::dune_static_assert(Dune::is_same<LA0::Residual,LA1::Residual>::true,
                                 "Received two local assemblers which are non-compatible "
                                 "due to different residual vector types");
        Dune::dune_static_assert(Dune::is_same<LA0::Real,LA1::Real>::true,
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


      //! Constructor with empty constraints
      OneStepLocalAssembler (LA0 & la0_, LA1 & la1_, OneStepParameters & method_, 
                             Residual & const_residual_) 
        : la0(la0_), la1(la1_), method(method_), const_residual(const_residual_), 
          time(0.0), weight(1.0), stage(0)
          pattern_engine(*this), residual_engine(*this), jacobian_engine(*this)
      { static_checks(); }

      //! Constructor for non trivial constraints
      OneStepLocalAssembler (LA0 & la0_, LA1 & la1_, OneStepParameters & method_, 
                             Residual & const_residual_, const CU& cu_, const CV& cv_) 
        : Base(cu_, cv_), 
          la0(la0_), la1(la1_), method(method_), const_residual(const_residual_), 
          time(0.0), weight(1.0),
          pattern_engine(*this), residual_engine(*this), jacobian_engine(*this)
      { static_checks(); }

      //! Notifies the local assembler about the current time of
      //! assembling. Should be called before assembling if the local
      //! operator has time dependencies.
      void setTime(Real time_){
        lop.setTime(time_);
        time = time_;
      }

      //! Notifies the assembler about the current weight of assembling.
      void setWeight(RangeField weight_){
        weight = weight_;
      }

      //! Set the current stage of the one step scheme
      void setStage(int stage_){
        stage = stage_;
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
      LocalResidualAssemblerEngine & localResidualAssemblerEngine
      (Residual & r, const Solution & x)
      {
        residual_engine.setResidual(r);
        residual_engine.setConstResidual(const_residual);
        residual_engine.setSolution(x);
        return residual_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalJacobianAssemblerEngine & localJacobianAssemblerEngine
      (Jacobian & a, const Solution & x)
      {
        jacobian_engine.setJacobian(a);
        jacobian_engine.setSolution(x);
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
      const OneStepParameters & method;

      //! The constant part of the residual
      Residual & const_residual;
      
      //! The current time of assembling
      Real time;

      //! The current weight of assembling
      RangeField weight;

      //! The current stage of the one step scheme
      int stage;

      //! The engine member objects
      //! @{
      LocalPatternAssemblerEngine  pattern_engine;
      LocalResidualAssemblerEngine residual_engine;
      LocalJacobianAssemblerEngine jacobian_engine;
      //! @}
    };

  };
};
#endif
