#ifndef DUNE_PDELAB_ONESTEP_OPERATOR_HH
#define DUNE_PDELAB_ONESTEP_OPERATOR_HH

#include <dune/pdelab/gridoperator/onestep/localassembler.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief A standard grid operator implementation suitable for the
       combination with a one step time stepping method.

       
       \tparam GFSU GridFunctionSpace for ansatz functions
       \tparam GFSV GridFunctionSpace for test functions
       \tparam MB The matrix backend to be used for representation of the jacobian
       \tparam DF The domain field type of the operator
       \tparam RF The range field type of the operator
       \tparam ST The type of the sub triangulation
       \tparam nonoverlapping_mode Switch for nonoverlapping grids
       \tparam CU   Constraints maps for the individual dofs (trial space)
       \tparam CV   Constraints maps for the individual dofs (test space)

    */
    template<typename GO0, typename GO1>
    class OneStepGridOperator
    {
    public:

      //! The matrix backend to be used for the jacobian matrix
      typedef typename GO0::MatrixBackend MatrixBackend;

      //! The sparsity pattern container for the jacobian matrix
      typedef typename GO0::Pattern Pattern;

      //! The domain and range types of the operator
      //! @{
      typedef typename GO0::DomainField DomainField;
      typedef typename GO0::RangeField RangeField;
      typedef typename GO0::Jacobian Jacobian;
      typedef typename GO0::Range Range;
      typedef typename GO0::Domain Domain;
      //! @}

      template <typename MFT>
      struct MatrixContainer{
        typedef Jacobian Type;
      };

      //! The global UDG assembler type
      typedef typename GO0::Assembler Assembler;

      //! The local assembler types of the subordinate grid operators
      //! @{
      typedef typename GO0::LocalAssembler LocalAssemblerDT0;
      typedef typename GO1::LocalAssembler LocalAssemblerDT1;
      //! @}

      //! The local UDG assembler type
      typedef OneStepLocalAssembler<LocalAssemblerDT0,LocalAssemblerDT1> LocalAssembler;

      //! The type of the one step method parameters
      typedef typename LocalAssembler::OneStepParameters OneStepParameters;

      //! Constructor for non trivial constraints
      OneStepGridOperator(const GO0 & go0_, const GO1 & go1_)
        : global_assembler(go0_.assembler()), 
          la0(go0_.localAssembler()), la1(go1_.localAssembler()),
          const_residual(go0_.testGridFunctionSpace()),
          local_assembler(la0,la1, const_residual)
      {}

      //! Get the trial grid function space
      const GFSU& trialGridFunctionSpace() const
      {
        return global_assembler.trialGridFunctionSpace();
      }

      //! Get the test grid function space
      const GFSV& testGridFunctionSpace() const
      {
        return global_assembler.testGridFunctionSpace();
      }

      //! Get dimension of space u
      typename GFSU::Traits::SizeType globalSizeU () const
      {
        return trialGridFunctionSpace().globalSize();
      }

      //! Get dimension of space v
      typename GFSV::Traits::SizeType globalSizeV () const
      {
        return testGridFunctionSpace().globalSize();
      }

      Assembler & assembler(){ return global_assembler; }

      LocalAssembler & localAssembler(){ return local_assembler; }

      //! Fill pattern of jacobian matrix
      void fill_pattern(Pattern & p) const {
        typedef typename LocalAssembler::LocalPatternAssemblerEngine PatternEngine;
        PatternEngine & pattern_engine = local_assembler.localPatternAssemblerEngine(p);
        global_assembler.assemble(pattern_engine);
      }

      //! Assemble constant part of residual
      void preStage(unsigned int stage, const std::vector<Domain*> & x){
        typedef typename LocalAssembler::LocalPreStageAssemblerEngine PreStageEngine;
        local_assembler.setStage(stage);
        PreStageEngine & prestage_engine = local_assembler.localPreStageAssemblerEngine(stage);
        global_assembler.assemble(prestage_engine);
      }

      //! Assemble residual
      void residual(const Domain & x, Range & r) const {
        typedef typename LocalAssembler::LocalResidualAssemblerEngine ResidualEngine;
        ResidualEngine & residual_engine = local_assembler.localResidualAssemblerEngine(r,x);
        global_assembler.assemble(residual_engine);
      }

      //! Assemble jacobian
      void jacobian(const Domain & x, Jacobian & a) const {
        typedef typename LocalAssembler::LocalJacobianAssemblerEngine JacobianEngine;
        JacobianEngine & jacobian_engine = local_assembler.localJacobianAssemblerEngine(a,x);
        global_assembler.assemble(jacobian_engine);
      }
      
      //! parametrize assembler with a time-stepping method
      void preStep (const TimeSteppingParameterInterface<TReal>& method_, TReal time_, TReal dt_)
      {
        local_assembler.setMethod(method_);
        local_assembler.setTime(time_);
        la0.preStep(time,dt,method->s());
        la1.preStep(time,dt,method->s());
      }

      //! to be called after step is completed
      void postStep ()
      {
        la0.postStep();
        la1.postStep();
      }

      //! to be called after stage is completed
      void postStage ()
      {
        la0.postStage();
        la1.postStage();
      }

      //! to be called once before each stage
      TReal suggestTimestep (TReal dt) const
      {
        TReal suggested_dt = la.suggestTimestep(dt);
        if (gfsu.gridview().comm().size()>1)
          suggested_dt =  gfsu.gridview().comm().min(suggested_dt);
        return suggested_dt;
      }

    private:
      Assembler global_assembler;
      LocalAssemblerDT0 & la0;
      LocalAssemblerDT1 & la1;
      Range const_residual;
      mutable LocalAssembler local_assembler;
    };

  };
};
#endif
