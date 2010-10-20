// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_MULTISTEP_GRIDOPERATORSPACE_HH
#define DUNE_PDELAB_MULTISTEP_GRIDOPERATORSPACE_HH

#include <cmath>
#include <cstddef>

#include <dune/common/shared_ptr.hh>
#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>

#include <dune/pdelab/gridfunctionspace/constraints.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include <dune/pdelab/localoperator/weightedsum.hh>
#include <dune/pdelab/multistep/parameter.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup MultiStepMethods Multi-Step Methods
    //! \ingroup PDELab
    //! \{

    template<typename TReal, typename LOPs>
    class MultiStepGridOperatorSpaceBase
    {
    protected:
      typedef typename ForEachType<AddRefTypeEvaluator, LOPs>::Type LOPRefs;

      WeightedSumLocalOperator<TReal, LOPs> sumLOP;

      MultiStepGridOperatorSpaceBase(const LOPRefs& lops)
        : sumLOP(lops)
      { }
    };

    ////////////////////////////////////////////////////////////////////////
    //
    //  Operator
    //

    //! Generic assembler for multi-step nth-order problems
    /**
     * \tparam TReal               type to represent time values (and
     *                             coefficients of time-stepping schemes)
     * \tparam R                   type that stores a residual vector
     * \tparam GFSU                GridFunctionSpace for ansatz functions
     * \tparam GFSV                GridFunctionSpace for test functions
     * \tparam LOPs                Tuple of LocalOperators.  Element n is the
     *                             local operator for r_n.
     * \tparam CU                  assembled constraints for the space U
     * \tparam CV                  assembled constraints for the space V
     * \tparam B                   linear algebra backend
     * \tparam nonoverlapping_mode switch to assemble for nonoverlapping grids
     *
     * The local operators should be derived from
     * InstationaryLocalOperatorDefaultMethods or comply to its interface by
     * other means.
     */
    template<typename TReal,
             typename R,
             typename GFSU,
             typename GFSV,
             typename LOPs,
             typename CU=EmptyTransformation,
             typename CV=EmptyTransformation,
             typename B=StdVectorFlatMatrixBackend,
             bool nonoverlapping_mode=false>
    class MultiStepGridOperatorSpace
      : private MultiStepGridOperatorSpaceBase<TReal, LOPs>,
        public GridOperatorSpace<GFSU, GFSV,
                                 WeightedSumLocalOperator<TReal, LOPs>,
                                 CU, CV, B, nonoverlapping_mode>
    {
      typedef MultiStepGridOperatorSpaceBase<TReal, LOPs> Base;
      typedef GridOperatorSpace<GFSU, GFSV,
                                WeightedSumLocalOperator<TReal, LOPs>,
                                CU, CV, B, nonoverlapping_mode> SGOS;

      typedef typename Base::LOPRefs LOPRefs;

    public:
      // export order of this problem
      static const std::size_t order = tuple_size<LOPs>::value - 1;

    private:
      using Base::sumLOP;
      // parameter class
      typedef MultiStepParameterInterface<TReal, order> Parameters;
      const Parameters *parameters;
      // constant part of the residual
      shared_ptr<R> r0;
      // time at end of step
      TReal tn;
      // current time step size
      TReal dt;

      // set the weights for the given step within the multi-step scheme
      // parameters and dt must have been setup before use.
      void setWeights(std::size_t step) {
        for(std::size_t i = 0; i <= order; ++i)
          sumLOP.setWeight(parameters->alpha(step,i)/std::pow(dt,int(i)), i);
      }

    public:
      //! construct
      /**
       * \param parameters_ Parameters for the multi-step scheme.
       * \param gfsu        Trial GridFunctionsSpace.
       * \param gfsv        Test GridFunctionsSpace.
       * \param lops        Tuple of references to the local operators.
       *
       * This constructor uses empty constraints.
       */
      MultiStepGridOperatorSpace
      ( const Parameters& parameters_,
        const GFSU& gfsu, const GFSV& gfsv,
        const LOPRefs& lops)
        : Base(lops), SGOS(gfsu, gfsv, sumLOP),
          parameters(&parameters_),
          r0(), tn(0), dt(1)
      {
        // This may not be the time step that is used in the end, but it is
        // sufficient for fill_pattern() to work
        setWeights(0);
      }

      //! construct
      /**
       * \param parameters_ Parameters for the multi-step scheme.
       * \param gfsu        Trial GridFunctionsSpace.
       * \param gfsv        Test GridFunctionsSpace.
       * \param lops        Tuple of references to the local operators.
       * \param cu          Constraints on the trial GridFunctionsSpace.
       * \param cv          Constraints on the test GridFunctionsSpace.
       */
      MultiStepGridOperatorSpace
      ( const Parameters& parameters_,
        const GFSU& gfsu, const CU& cu,
        const GFSV& gfsv, const CV& cv,
        const LOPRefs& lops)
        : Base(lops), SGOS(gfsu, cu, gfsv, cv, sumLOP),
          parameters(&parameters_),
          r0(), tn(0), dt(1)
      {
        // This may not be the time step that is used in the end, but it is
        // sufficient for fill_pattern() to work
        setWeights(0);
      }

      //! construct
      /**
       * \param gfsu Trial GridFunctionsSpace.
       * \param gfsv Test GridFunctionsSpace.
       * \param lops Tuple of references to the local operators.
       *
       * This constructor uses empty constraints and does not initialize the
       * parameters.
       */
      MultiStepGridOperatorSpace
      ( const GFSU& gfsu, const GFSV& gfsv,
        const LOPRefs& lops)
        : Base(lops), SGOS(gfsu, gfsv, sumLOP),
          parameters(0),
          r0(), tn(0), dt(1)
      { }

      //! construct
      /**
       * \param gfsu Trial GridFunctionsSpace.
       * \param gfsv Test GridFunctionsSpace.
       * \param lops Tuple of references to the local operators.
       * \param cu   Constraints on the trial GridFunctionsSpace.
       * \param cv   Constraints on the test GridFunctionsSpace.
       *
       * This constructor does not initialize the parameters.
       */
      MultiStepGridOperatorSpace
      ( const GFSU& gfsu, const CU& cu,
        const GFSV& gfsv, const CV& cv,
        const LOPRefs& lops)
        : Base(lops), SGOS(gfsu, cu, gfsv, cv, sumLOP),
          parameters(0),
          r0(), tn(0), dt(1)
      { }

      //! Set the methods to use
      void setMethod(const Parameters& parameters_) {
        parameters = &parameters_;
        // This may not be the time step that is used in the end, but it is
        // sufficient for fill_pattern() to work
        setWeights(0);
      }

      //! prepare for doing a step
      /**
       * Set the current time, and calculate the constant part of the residual.
       *
       * \param time      Time at beginning of step.
       * \param dt_       Size of step.
       * \param oldvalues The old values.  Must support the expression
       *                  *oldvalues[i] for 0 <= i < steps.  *oldvalues[0]
       *                  should yield a reference to a vector with the values
       *                  at time time, *oldvalues[1] should yield a reference
       *                  to a vector with the values at time time-dt and
       *                  generally *oldvalues[i] should yield a reference to
       *                  a vector with the old values at time time-i*dt.
       *
       * This method calls preStep(time,dt,1) and preStage(time+dt, 1) on all
       * local operators.
       */
      template<typename OldValues>
      void preStep(TReal time, TReal dt_, const OldValues& oldvalues)
      {
        dt = dt_;
        tn = time+dt;

        // allocate constant part of residual
        r0.reset(new R(this->testGridFunctionSpace(), 0.0));

        // prepare the local operators
        sumLOP.preStep(time, dt, 1);
        sumLOP.preStage(tn, 1);

        for(unsigned step = 1; step <= parameters->steps(); ++step) {
          setWeights(step);
          sumLOP.setTime(tn-(step-1)*dt);
          // residual is additive
          SGOS::residual(*oldvalues[step-1], *r0);
        }
        // reset weights and time to the end of the step
        setWeights(0);
        sumLOP.setTime(tn);
      }

      //! to be called after step is completed
      /**
       * Invokes postStage() and postStep() on the local operators.
       */
      void postStep ()
      {
        sumLOP.postStage();
        sumLOP.postStep();
        // free constant part of residual
        r0.reset();
      }

      //! Interpolate constrained values
      /**
       * \param xold  Vector with the old values.  Used to obtain the
       *              non-constrained values.
       * \param f     Function to evaluate to obtain the contrained values.
       *              The time on the function should be set apropriately
       *              before calling this method.
       * \param x     Where to store the combination of xold and the
       *              interpolated values.
       *
       * \note xold and x may not refer to the same object.
       */
      template<typename F, typename X>
      void interpolate (const X& xold, F& f, X& x) const
      {
        // make x obey the boundary values
        Dune::PDELab::interpolate(f,this->trialGridFunctionSpace(),x);

        // copy non-constrained dofs from old time step
        copy_nonconstrained_dofs(this->testConstraints(),xold,x);
      }

      //////////////////////////////////////////////////////////////////////
      //
      // residual methods
      //

      //! generic evaluation of residual
      /**
       * \param r residual (needs to be cleared before this method is called)
       */
      template<typename X>
      void residual (const X& x, R& r) const
      {
        // copy constant part of residual
        r += *r0;

        SGOS::residual(x, r);

        // set residual to zero on constrained dofs
        constrain_residual(this->testConstraints(),r);
      }
    };

    //! \} group MultiStepMethods
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTISTEP_GRIDOPERATORSPACE_HH
