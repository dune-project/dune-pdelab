// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_MULTISTEP_GRIDOPERATORSPACE_HH
#define DUNE_PDELAB_MULTISTEP_GRIDOPERATORSPACE_HH

#include <cmath>
#include <cstddef>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/typetraits.hh>

#include <dune/pdelab/gridfunctionspace/constraints.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include <dune/pdelab/localoperator/weightedsum.hh>
#include <dune/pdelab/multistep/cache.hh>
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
          sumLOP.setWeight(parameters->alpha(step,i) / std::pow(dt, TReal(i)),
                           i);
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

    //! Caching MultiStepGridOperatorSpace
    /**
     * \tparam Step                Type used for step numbers.
     * \tparam Time                Type to represent temporal values (and
     *                             coefficients of time-stepping schemes)
     * \tparam Coeffs              Coefficient type used for residual and
     *                             ansatz vectors and the matrices.
     * \tparam GFSU                GridFunctionSpace for ansatz functions.
     * \tparam GFSV                GridFunctionSpace for test functions.
     * \tparam LOPs                Tuple of LocalOperators.  Element n is the
     *                             type of the local operator for r_n.
     * \tparam CU                  Constraints maps for the individual dofs
     *                             (trial space)
     * \tparam CV                  Constraints maps for the individual dofs
     *                             (test space)
     * \tparam B                   The matrix backend.
     * \tparam nonoverlapping_mode Indicates whether assembling is done for
     *                             overlap cells
     */
    template<class Step, class Time, class Coeffs,
             class GFSU, class GFSV, class LOPs,
             class CU=EmptyTransformation, class CV=EmptyTransformation,
             class B=StdVectorFlatMatrixBackend,
             bool nonoverlapping_mode=false>
    class CachedMultiStepGridOperatorSpace :
      public GridOperatorBase<GFSU,GFSV,CU,CV,B>
    {
    public:
      //! contruct a matrix type from the backend
      /**
       * \tparam E Type of the matrix coefficients.
       */
      template<typename E>
      struct MatrixContainer
      {
        //! \brief define Type as the Type of a Matrix of E's
        typedef typename B::template Matrix<CachedMultiStepGridOperatorSpace,
                                            E> Type;
      private:
        MatrixContainer () {}
      };

    private:
      // extract useful types
      typedef GridOperatorBase<GFSU,GFSV,CU,CV,B> Base;

      typedef typename GFSU::template VectorContainer<Coeffs>::Type
        UnknownVector;
      typedef typename GFSV::template VectorContainer<Coeffs>::Type
        ResidualVector;
      typedef typename MatrixContainer<Coeffs>::Type Matrix;

      typedef MultiStepCache<ResidualVector, UnknownVector, Matrix, Step,
                             Time> Cache;

      typedef typename ForEachType<AddRefTypeEvaluator, LOPs>::Type LOPRefs;

      // tuple of shared_ptrs to the component Stationary GOS
      template<class LOP>
      class LOPToSharedPtrGOSTypeEvaluator {
        typedef GridOperatorSpace<GFSU, GFSV,
                                  typename remove_reference<LOP>::type,
                                  EmptyTransformation,
                                  EmptyTransformation,
                                  B, nonoverlapping_mode> GOS;
      public:
        typedef shared_ptr<GOS> Type;
        static Type apply(const typename remove_reference<LOP>::type &lop,
                          const GFSU& gfsu, const GFSV& gfsv)
        { return shared_ptr<GOS>(new GOS(gfsu, gfsv, lop)); }
      };
      typedef typename ForEachType<LOPToSharedPtrGOSTypeEvaluator, LOPs>::Type
        GOSPtrTuple;

    public:
      // export order of this problem
      static const std::size_t order = tuple_size<LOPs>::value - 1;

    private:
      typedef MultiStepParameterInterface<Time, order> Parameters;

      // member variables
      const Parameters *parameters;
      shared_ptr<Cache> cache;
      LOPRefs lops;
      GOSPtrTuple gosPtrs;
      // constant part of the residual
      shared_ptr<ResidualVector> r0;
      // number of the values at the end of the step
      Step currentStep;
      // time at end of step
      Time tn;
      // current time step size
      Time dt;

      // loop over local operator references
      ForEachValue<LOPRefs> lopLoop;
      // loop over stationary grid operator pointers
      ForEachValue<GOSPtrTuple> gosLoop;

    public:
      typedef typename Base::Traits Traits;

      //! construct
      /**
       * \param parameters_ Parameters for the multi-step scheme.
       * \param gfsu        Trial GridFunctionsSpace.
       * \param gfsv        Test GridFunctionsSpace.
       * \param lops_       Tuple of references to the local operators.
       *
       * This constructor uses empty constraints.
       */
      CachedMultiStepGridOperatorSpace(const Parameters& parameters_,
                                       const GFSU& gfsu, const GFSV& gfsv,
                                       const LOPRefs& lops_) :
        Base(gfsu, gfsv),
        parameters(&parameters_), cache(new Cache), lops(lops_),
        gosPtrs(transformTuple<LOPToSharedPtrGOSTypeEvaluator>
                (lops, gfsu, gfsv)),
        r0(), currentStep(0), tn(0), dt(1),
        lopLoop(lops), gosLoop(gosPtrs)
      { }

      //! construct
      /**
       * \param parameters_ Parameters for the multi-step scheme.
       * \param gfsu        Trial GridFunctionsSpace.
       * \param gfsv        Test GridFunctionsSpace.
       * \param lops_       Tuple of references to the local operators.
       * \param cu          Constraints on the trial GridFunctionsSpace.
       * \param cv          Constraints on the test GridFunctionsSpace.
       */
      CachedMultiStepGridOperatorSpace(const Parameters& parameters_,
                                       const GFSU& gfsu, const CU& cu,
                                       const GFSV& gfsv, const CV& cv,
                                       const LOPRefs& lops_) :
        Base(gfsu, cu, gfsv, cv),
        parameters(&parameters_), cache(new Cache), lops(lops_),
        gosPtrs(transformTuple<LOPToSharedPtrGOSTypeEvaluator>
                (lops, gfsu, gfsv)),
        r0(), currentStep(0), tn(0), dt(1),
        lopLoop(lops), gosLoop(gosPtrs)
      { }

      //! construct
      /**
       * \param gfsu  Trial GridFunctionsSpace.
       * \param gfsv  Test GridFunctionsSpace.
       * \param lops_ Tuple of references to the local operators.
       *
       * This constructor uses empty constraints and does not initialize the
       * parameters.
       */
      CachedMultiStepGridOperatorSpace(const GFSU& gfsu, const GFSV& gfsv,
                                       const LOPRefs& lops_) :
        Base(gfsu, gfsv),
        parameters(0), cache(new Cache), lops(lops_),
        gosPtrs(transformTuple<LOPToSharedPtrGOSTypeEvaluator>
                (lops, gfsu, gfsv)),
        r0(), currentStep(0), tn(0), dt(1),
        lopLoop(lops), gosLoop(gosPtrs)
      { }

      //! construct
      /**
       * \param gfsu  Trial GridFunctionsSpace.
       * \param gfsv  Test GridFunctionsSpace.
       * \param lops_ Tuple of references to the local operators.
       * \param cu    Constraints on the trial GridFunctionsSpace.
       * \param cv    Constraints on the test GridFunctionsSpace.
       *
       * This constructor does not initialize the parameters.
       */
      CachedMultiStepGridOperatorSpace(const GFSU& gfsu, const CU& cu,
                                       const GFSV& gfsv, const CV& cv,
                                       const LOPRefs& lops_) :
        Base(gfsu, cu, gfsv, cv),
        parameters(0), cache(new Cache), lops(lops_),
        gosPtrs(transformTuple<LOPToSharedPtrGOSTypeEvaluator>
                (lops, gfsu, gfsv)),
        r0(), currentStep(0), tn(0), dt(1),
        lopLoop(lops), gosLoop(gosPtrs)
      { }

      //! Set the methods to use
      void setMethod(const Parameters& parameters_) {
        parameters = &parameters_;
      }

      //! get the cache object
      shared_ptr<Cache> getCache() const { return cache; }
      //! set the cache object to use
      void setCache(const shared_ptr<Cache> &cache_) { cache = cache_; }

    private:
      // LOCAL VISITORS

      // call preStep() on the local operators
      class PreStepVisitor {
        Time time;
        Time dt;

      public:
        PreStepVisitor(Time time_, Time dt_) : time(time_), dt(dt_) { }

        template<class LOP>
        void visit(LOP& lop) { lop.preStep(time, dt, 1); }
      };

      // call postStep() on the local operators
      struct PostStepVisitor {
        template<class LOP>
        void visit(LOP& lop) { lop.postStep(); }
      };

        // sumLOP.preStage(tn, 1);

      // call preStage() on the local operators
      class PreStageVisitor {
        Time tn;

      public:
        explicit PreStageVisitor(Time tn_) : tn(tn_) { }

        template<class LOP>
        void visit(LOP& lop) { lop.preStage(tn, 1); }
      };

      // call postStage() on the local operators
      struct PostStageVisitor {
        template<class LOP>
        void visit(LOP& lop) { lop.postStage(); }
      };

      // call setTime() on the local operators
      class SetTimeVisitor {
        Time time;

      public:
        explicit SetTimeVisitor(Time time_) : time(time_) { }

        template<class LOP>
        void visit(LOP& lop) { lop.setTime(time); }
      };

      // GLOBAL VISITORS

      // collect patterns of the current step
      template<class Pattern>
      class PatternVisitor {
        const CachedMultiStepGridOperatorSpace& mgos;
        Pattern &pattern;
        std::size_t currentElem;

      public:
        PatternVisitor(const CachedMultiStepGridOperatorSpace& mgos_,
                       Pattern &pattern_) :
          mgos(mgos_), pattern(pattern_), currentElem(0)
        { }

        template<class GOSPtr>
        void visit(const GOSPtr& gosPtr) {
          if(mgos.parameters->alpha(0, currentElem) != 0.0)
            gosPtr->fill_pattern(pattern);
          ++currentElem;
        }
      };

      // add up the residual valuess of the given step
      class ResidualValueVisitor {
        const CachedMultiStepGridOperatorSpace& mgos;
        ResidualVector &residual;
        unsigned backStep;
        std::size_t currentElem;

      public:
        ResidualValueVisitor(const CachedMultiStepGridOperatorSpace& mgos_,
                             unsigned backStep_, ResidualVector &residual_) :
          mgos(mgos_), residual(residual_), backStep(backStep_), currentElem(0)
        { }

        template<class GOSPtr>
        void visit(const GOSPtr& gosPtr) {
          if(mgos.parameters->alpha(backStep, currentElem) != 0.0)
            residual.axpy
              ( mgos.parameters->alpha(backStep, currentElem) /
                  std::pow(mgos.dt, Time(currentElem)),
                *mgos.getResidualValue(*gosPtr, currentElem, backStep));
          ++currentElem;
        }
      };

      // add up the Jacobians of the given step
      class JacobianVisitor {
        const CachedMultiStepGridOperatorSpace& mgos;
        Matrix &composedJacobian;
        unsigned backStep;
        std::size_t currentElem;

      public:
        JacobianVisitor(const CachedMultiStepGridOperatorSpace& mgos_,
                        unsigned backStep_, Matrix &composedJacobian_) :
          mgos(mgos_), composedJacobian(composedJacobian_),
          backStep(backStep_), currentElem(0)
        { }

        template<class GOSPtr>
        void visit(const GOSPtr& gosPtr) {
          if(mgos.parameters->alpha(backStep, currentElem) != 0.0)
            composedJacobian.axpy
              ( mgos.parameters->alpha(backStep, currentElem) /
                  std::pow(mgos.dt, Time(currentElem)),
                *mgos.getJacobian(*gosPtr, currentElem, backStep));
          ++currentElem;
        }
      };

      // add up the zero-residuals of the given step
      class ZeroResidualVisitor {
        const CachedMultiStepGridOperatorSpace& mgos;
        ResidualVector &residual;
        unsigned backStep;
        std::size_t currentElem;

      public:
        ZeroResidualVisitor(const CachedMultiStepGridOperatorSpace& mgos_,
                            unsigned backStep_, ResidualVector &residual_) :
          mgos(mgos_), residual(residual_), backStep(backStep_), currentElem(0)
        { }

        template<class GOSPtr>
        void visit(const GOSPtr& gosPtr) {
          if(mgos.parameters->alpha(backStep, currentElem) != 0.0)
            residual.axpy
              ( mgos.parameters->alpha(backStep, currentElem) /
                  std::pow(mgos.dt, Time(currentElem)),
                *mgos.getZeroResidual(*gosPtr, currentElem, backStep));
          ++currentElem;
        }
      };

      //! get a zero-residual from a particular stationary GOS
      /**
       * Does not yet support fetching from the cache.
       */
      template<typename GOS>
      shared_ptr<const ResidualVector>
      getZeroResidual(const GOS& gos, std::size_t order,
                      unsigned backStep) const
      {
        shared_ptr<ResidualVector>
          result(new ResidualVector(this->testGridFunctionSpace(), 0));

        if(cache->getPolicy()->hasPureLinearAlpha(order, backStep))
          gos.zero_residual(*result);
        else
          gos.residual(ResidualVector(this->trialGridFunctionSpace(), 0),
                       *result);

        return result;
      }

      //! get a Jacobian from a particular stationary GOS
      /**
       * Does not yet support fetching from the cache (thats what the
       * otherwise unneeded parameter order is for).
       */
      template<typename GOS>
      shared_ptr<const Matrix>
      getJacobian(const GOS& gos, std::size_t order, unsigned backStep) const {
        // AFFINE: This method is only meaningful in the affine case, so don't
        // even consider non-linearity
        if(!cache->getPolicy()->isAffine(order, currentStep-Step(backStep)))
          DUNE_THROW(InvalidStateException,
                     "CachedMultiStepGridOperatorSpace::getJacobian(): This "
                     "method is meaningful only for affine operator "
                     "R" << order);

        shared_ptr<Matrix> result(new Matrix(*this));
        *result = 0;

        gos.jacobian(ResidualVector(this->trialGridFunctionSpace(), 0),
                     *result);

        return result;
      }

      //! get a residual value from a particular stationary GOS
      /**
       * Will try to fetch the value from the cache first.  If the value has
       * to be computed, will try to store the value in the cache.
       */
      template<typename GOS>
      shared_ptr<const ResidualVector>
      getResidualValue(const GOS& gos, std::size_t order,
                       unsigned backStep) const
      {
        try {
          return cache->getResidualValue(order, currentStep-Step(backStep));
        }
        catch(const NotInCache&) {
          shared_ptr<ResidualVector>
            result(new ResidualVector(this->testGridFunctionSpace(), 0));

          gos.residual(*cache->getUnknowns(currentStep-Step(backStep)),
                       *result);
          cache->setResidualValue(order, currentStep-Step(backStep), result);
          return result;
        }
      }

      //! get the composed jacobian
      /**
       * fetch it from the cache if possible, otherwise compute it and store
       * it in teh cache along the way.
       */
      shared_ptr<const Matrix> getComposedJacobian() const {
        // AFFINE: This method is only meaningful in the affine case, so don't
        // even consider non-linearity
        if(!cache->getPolicy()->isComposedAffine(currentStep))
          DUNE_THROW(InvalidStateException,
                     "CachedMultiStepGridOperatorSpace::getComposeJacobian(): "
                     "This method is meaningful only for affine operators");

        try { return cache->getComposedJacobian(currentStep); }
        catch(const NotInCache&) {
          shared_ptr<Matrix> result(new Matrix(*this));
          *result = 0;

          { JacobianVisitor visitor(*this, 0, *result);
            gosLoop.apply(visitor); }

          // apply constraints
          typedef typename CV::const_iterator global_row_iterator;
          for(global_row_iterator cit = this->pconstraintsv->begin();
              cit != this->pconstraintsv->end(); ++cit)
            set_trivial_row(cit->first,cit->second,*result);

          cache->setComposedJacobian(currentStep, result);
          return result;
        }
      }

    public:
      //! prepare for doing a step
      /**
       * Set the current time, and calculate the constant part of the residual.
       *
       * \param step      Number of step, and step index of the values
       *                  computed during that step (i.e. the ones at the end
       *                  of that step).
       * \param startTime Time at beginning of step.
       * \param dt_       Size of step.
       *
       * This method calls preStep(time,dt,1) and preStage(time+dt, 1) on all
       * local operators.
       */
      void preStep(Step step, Time startTime, Time dt_)
      {
        currentStep = step;
        dt = dt_;
        tn = startTime+dt;

        // allocate constant part of residual
        r0.reset(new ResidualVector(this->testGridFunctionSpace(), 0.0));

        // prepare the local operators
        { PreStepVisitor visitor(tn-dt, dt); lopLoop.apply(visitor); }
        { PreStageVisitor visitor(tn); lopLoop.apply(visitor); }

        for(unsigned backStep = 1; backStep <= parameters->steps(); ++backStep)
        {
          { SetTimeVisitor visitor(tn-backStep*dt); lopLoop.apply(visitor); }
          { ResidualValueVisitor visitor(*this, backStep, *r0);
            gosLoop.apply(visitor); }
        }
        // reset time to the end of the step
        { SetTimeVisitor visitor(tn); lopLoop.apply(visitor); }

        if(cache->getPolicy()->isComposedAffine(currentStep)) {
          // AFFINE: also include the contribution of the zero-residuals of
          // the current time step
          { ZeroResidualVisitor visitor(*this, 0, *r0);
            gosLoop.apply(visitor); }

          // set residual to zero on constrained dofs
          Dune::PDELab::constrain_residual(*this->pconstraintsv,*r0);
        }
        else {
          // NON-LINEAR: nothing to be done
        }
      }

      //! to be called after step is completed
      /**
       * Invokes postStage() and postStep() on the local operators.
       */
      void postStep ()
      {
        { PostStageVisitor visitor; lopLoop.apply(visitor); }
        { PostStepVisitor visitor; lopLoop.apply(visitor); }
        // free constant part of residual
        r0.reset();
      }

      /**\brief Construct global sparsity pattern from local description

         This function can be called by the Matrix to get the sparsity pattern.
         Assumes that the pattern is initially empty.
      */
      template<typename Pattern>
      void fill_pattern (Pattern& globalpattern) const
      {
        PatternVisitor<Pattern> visitor(*this, globalpattern);
        gosLoop.apply(visitor);
      }

      //! generic evaluation of residual
      /**
       * \param r residual (needs to be cleared before this method is called)
       */
      template<typename X, typename R>
      void residual (const X& x, R& r) const
      {
        r += *r0;
        if(cache->getPolicy()->isComposedAffine(currentStep))
          // AFFINE: r0 already contains zero-residual part for step n and is
          // already constrained
          jacobian_apply(x,r);
        else {
          // NON-LINEAR: r0 contains contributions from older time-step only
          DUNE_THROW(NotImplemented,
                     "CachedMultiStepGridOperatorSpace::residual(): Not "
                     "implemented for non-linear operators");

          // set residual to zero on constrained dofs
          Dune::PDELab::constrain_residual(*this->pconstraintsv,r);
        }
      }

      //! \brief Assemble constant part of residual for operators with purely
      //!        linear alpha_*() methods
      /**
       * The result of this method is identical to calling residual() with
       * parameter \c x initialized to 0.  However, this method works only
       * when the local operator has purely linear alpha_volume(),
       * alpha_skeleton(), alpha_boundary() and alpha_volume_post_skeleton()
       * methods.  Note that this is a stronger demand than requiring the
       * operator to be affine -- affine operators may still have an affine
       * shift in their alpha_*() methods.  For this method to work any affine
       * shift must be implemented in the lambda_*() methods.
       *
       * \note This method is meaningless in the context of non-linear
       *       operators.
       *
       * \param r residual (needs to be cleared before this method is called)
       */
      template<typename R>
      void zero_residual(R& r) const {
        // AFFINE: This method is only meaningful in the affine case, so don't
        // even consider non-linearity
        if(!cache->getPolicy()->isComposedAffine(currentStep))
          DUNE_THROW(InvalidStateException,
                     "CachedMultiStepGridOperatorSpace::zero_residual(): "
                     "This method is meaningful only for affine operators");

        r += *r0;
        // AFFINE: r0 already contains zero-residual part for step n and is
        // already constrained
      }

      //! generic application of Jacobian
      template<typename X, typename Y>
      void jacobian_apply (X& x, Y& y) const
      {
        if(cache->getPolicy()->isComposedAffine(currentStep))
          // AFFINE: The jacobian already applies suitable constraints
          getComposedJacobian()->umv(x,y);
        else {
          // NON-LINEAR: r0 contains contributions from older time-step only
          DUNE_THROW(NotImplemented,
                     "CachedMultiStepGridOperatorSpace::jacobian_apply(): Not "
                     "implemented for non-linear operators");

          // set residual to zero on constrained dofs
          Dune::PDELab::constrain_residual(*this->pconstraintsv,y);
        }
      }

      //! generic assembly of Jacobian
      /**
       * \param x Where (in the space spanned by the dofs) to evaluate the Jacobian
       * \param a Jacobian (needs to be cleared before passed to this method)
       */
      template<typename X, typename A>
      void jacobian (const X& x, A& a) const
      {
        if(cache->getPolicy()->isComposedAffine(currentStep))
          // AFFINE: The jacobian already applies suitable constraints
          a += *getComposedJacobian();
        else {
          // NON-LINEAR
          DUNE_THROW(NotImplemented,
                     "CachedMultiStepGridOperatorSpace::jacobian(): Not "
                     "implemented for non-linear operators");

          // NON-LINEAR: (possiply) need to apply constraints
          typedef typename CV::const_iterator global_row_iterator;
          for(global_row_iterator cit = this->pconstraintsv->begin();
              cit != this->pconstraintsv->end(); ++cit)
            set_trivial_row(cit->first,cit->second,a);
        }
      }

    };

    //! \} group MultiStepMethods
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTISTEP_GRIDOPERATORSPACE_HH
