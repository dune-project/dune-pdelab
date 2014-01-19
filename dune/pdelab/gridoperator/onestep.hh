
#ifndef DUNE_PDELAB_ONESTEP_OPERATOR_HH
#define DUNE_PDELAB_ONESTEP_OPERATOR_HH

#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/gridoperator/onestep/localassembler.hh>
#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
#include <dune/pdelab/constraints/common/constraints.hh>

namespace Dune{
  namespace PDELab{

    template<typename GO0, typename GO1, bool implicit = true>
    class OneStepGridOperator
    {
    public:

      //! The sparsity pattern container for the jacobian matrix
      typedef typename GO0::Pattern Pattern;

      //! The global UDG assembler type
      typedef typename GO0::Traits::Assembler Assembler;

      //! The local assembler types of the subordinate grid operators
      //! @{
      typedef typename GO0::Traits::LocalAssembler LocalAssemblerDT0;
      typedef typename GO1::Traits::LocalAssembler LocalAssemblerDT1;
      //! @}

      //! The local assembler type
      typedef OneStepLocalAssembler<OneStepGridOperator,LocalAssemblerDT0,LocalAssemblerDT1> LocalAssembler;

      //! The BorderDOFExchanger
      typedef typename GO0::BorderDOFExchanger BorderDOFExchanger;

      //! The grid operator traits
      typedef Dune::PDELab::GridOperatorTraits
      <typename GO0::Traits::TrialGridFunctionSpace,
       typename GO0::Traits::TestGridFunctionSpace,
       typename GO0::Traits::MatrixBackend,
       typename GO0::Traits::DomainField,
       typename GO0::Traits::RangeField,
       typename GO0::Traits::JacobianField,
       typename GO0::Traits::TrialGridFunctionSpaceConstraints,
       typename GO0::Traits::TestGridFunctionSpaceConstraints,
       Assembler,
       LocalAssembler> Traits;

      //! The io types of the operator
      //! @{
      typedef typename Traits::Domain Domain;
      typedef typename Traits::Range Range;
      typedef typename Traits::Jacobian Jacobian;
      //! @}

      template <typename MFT>
      struct MatrixContainer{
        typedef Jacobian Type;
      };

      //! The type for real number e.g. time
      typedef typename LocalAssembler::Real Real;

      //! The type of the one step method parameters
      typedef typename LocalAssembler::OneStepParameters OneStepParameters;

      //! Constructor for non trivial constraints
      OneStepGridOperator(GO0 & go0_, GO1 & go1_)
        : global_assembler(go0_.assembler()),
          go0(go0_), go1(go1_),
          la0(go0_.localAssembler()), la1(go1_.localAssembler()),
          const_residual( go0_.testGridFunctionSpace() ),
          local_assembler(la0,la1, const_residual)
      {
        GO0::setupGridOperators(Dune::tie(go0_,go1_));
        if(!implicit)
          local_assembler.setDTAssemblingMode(LocalAssembler::DoNotAssembleDT);
      }

      //! Determines whether the time step size is multiplied to the
      //! mass term (first order time derivative) or the elliptic term
      //! (zero-th order time derivative).
      void divideMassTermByDeltaT()
      {
        if(!implicit)
          DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");
        local_assembler.setDTAssemblingMode(LocalAssembler::DivideOperator1ByDT);
      }
      void multiplySpatialTermByDeltaT()
      {
        if(!implicit)
          DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");
        local_assembler.setDTAssemblingMode(LocalAssembler::MultiplyOperator0ByDT);
      }

      //! Get the trial grid function space
      const typename Traits::TrialGridFunctionSpace& trialGridFunctionSpace() const
      {
        return global_assembler.trialGridFunctionSpace();
      }

      //! Get the test grid function space
      const typename Traits::TestGridFunctionSpace& testGridFunctionSpace() const
      {
        return global_assembler.testGridFunctionSpace();
      }

      //! Get dimension of space u
      typename Traits::TrialGridFunctionSpace::Traits::SizeType globalSizeU () const
      {
        return trialGridFunctionSpace().globalSize();
      }

      //! Get dimension of space v
      typename Traits::TestGridFunctionSpace::Traits::SizeType globalSizeV () const
      {
        return testGridFunctionSpace().globalSize();
      }

      Assembler & assembler() const { return global_assembler; }

      LocalAssembler & localAssembler() const { return local_assembler; }

      //! Fill pattern of jacobian matrix
      void fill_pattern(Pattern & p) const {
        if(implicit){
          typedef typename LocalAssembler::LocalPatternAssemblerEngine PatternEngine;
          PatternEngine & pattern_engine = local_assembler.localPatternAssemblerEngine(p);
          global_assembler.assemble(pattern_engine);
        } else {
          typedef typename LocalAssembler::LocalExplicitPatternAssemblerEngine PatternEngine;
          PatternEngine & pattern_engine = local_assembler.localExplicitPatternAssemblerEngine(p);
          global_assembler.assemble(pattern_engine);
        }
      }

      //! Assemble constant part of residual
      void preStage(unsigned int stage, const std::vector<Domain*> & x){
        if(!implicit){DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");}

        typedef typename LocalAssembler::LocalPreStageAssemblerEngine PreStageEngine;
        local_assembler.setStage(stage);
        PreStageEngine & prestage_engine = local_assembler.localPreStageAssemblerEngine(x);
        global_assembler.assemble(prestage_engine);
        //Dune::printvector(std::cout,const_residual.base(),"const residual","row",4,9,1);
      }

      //! Assemble residual
      void residual(const Domain & x, Range & r) const {
        if(!implicit){DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");}

        typedef typename LocalAssembler::LocalResidualAssemblerEngine ResidualEngine;
        ResidualEngine & residual_engine = local_assembler.localResidualAssemblerEngine(r,x);
        global_assembler.assemble(residual_engine);
        //Dune::printvector(std::cout,r.base(),"residual","row",4,9,1);
      }

      //! Assemble jacobian
      void jacobian(const Domain & x, Jacobian & a) const {
        if(!implicit){DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");}

        typedef typename LocalAssembler::LocalJacobianAssemblerEngine JacobianEngine;
        JacobianEngine & jacobian_engine = local_assembler.localJacobianAssemblerEngine(a,x);
        global_assembler.assemble(jacobian_engine);
        //printmatrix(std::cout,a.base(),"global stiffness matrix","row",9,1);
      }

      //! Assemble jacobian and residual simultaneously for explicit treatment
      void explicit_jacobian_residual(unsigned int stage, const std::vector<Domain*> & x,
                                      Jacobian & a, Range & r1, Range & r0)
      {
        if(implicit){DUNE_THROW(Dune::Exception,"This function should not be called in implicit mode");}

        local_assembler.setStage(stage);

        typedef typename LocalAssembler::LocalExplicitJacobianResidualAssemblerEngine
          ExplicitJacobianResidualEngine;

        ExplicitJacobianResidualEngine & jacobian_residual_engine
          = local_assembler.localExplicitJacobianResidualAssemblerEngine(a,r0,r1,x);

        global_assembler.assemble(jacobian_residual_engine);
      }

      //! Interpolate constrained values from given function f
      template<typename F, typename X>
      void interpolate (unsigned stage, const X& xold, F& f, X& x) const
      {
        // Set time in boundary value function
        f.setTime(local_assembler.timeAtStage(stage));

        go0.localAssembler().setTime(local_assembler.timeAtStage(stage));

        // Interpolate
        go0.interpolate(xold,f,x);

        // Copy non-constrained dofs from old time step
        Dune::PDELab::copy_nonconstrained_dofs(local_assembler.trialConstraints(),xold,x);
      }

      //! set time stepping method
      void setMethod (const TimeSteppingParameterInterface<Real>& method_)
      {
        local_assembler.setMethod(method_);
      }

      //! parametrize assembler with a time-stepping method
      void preStep (const TimeSteppingParameterInterface<Real>& method_, Real time_, Real dt_)
      {
        local_assembler.setMethod(method_);
        local_assembler.preStep(time_,dt_,method_.s());
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
      Real suggestTimestep (Real dt) const
      {
        Real suggested_dt = std::min(la0.suggestTimestep(dt),la1.suggestTimestep(dt));
        if (trialGridFunctionSpace().gridView().comm().size()>1)
          suggested_dt =  trialGridFunctionSpace().gridView().comm().min(suggested_dt);
        return suggested_dt;
      }

      void update()
      {
        go0.update();
        go1.update();
        const_residual = Range(go0.testGridFunctionSpace());
      }

      const typename Traits::MatrixBackend& matrixBackend() const
      {
        return go0.matrixBackend();
      }

    private:
      Assembler & global_assembler;
      GO0 & go0;
      GO1 & go1;
      LocalAssemblerDT0 & la0;
      LocalAssemblerDT1 & la1;
      Range const_residual;
      mutable LocalAssembler local_assembler;
    };

  }
}
#endif
