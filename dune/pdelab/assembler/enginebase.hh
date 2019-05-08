// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_ENGINEBASE_HH
#define DUNE_PDELAB_ASSEMBLER_ENGINEBASE_HH

#include <dune/pdelab/common/typetraits.hh>

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/flavor.hh>
#include <dune/pdelab/instationary/onestepmethods.hh>

namespace Dune {
  namespace PDELab {

    enum class DTScaling {multiply, divide};

    enum class Galerkin { disable = 0, enable = 1, automatic };

    static constexpr auto enableGalerkin = std::integral_constant<Galerkin,Galerkin::enable>{};
    static constexpr auto disableGalerkin = std::integral_constant<Galerkin,Galerkin::disable>{};
    static constexpr auto automaticGalerkin = std::integral_constant<Galerkin,Galerkin::automatic>{};




    template<bool instationary_ = false, Galerkin galerkin_ = Galerkin::automatic>
    struct DefaultEngineParametersBase
    {

      static constexpr bool instationary = instationary_;

      static constexpr bool fastDG = false;

      template<typename TrialSpace, typename TestSpace>
      static constexpr bool galerkin =
        galerkin_ == Galerkin::automatic
        ? std::is_same_v<TrialSpace,TestSpace>
        : bool(galerkin_);

    };


    template<bool instationary_ = false, Galerkin galerkin_ = Galerkin::automatic>
    struct DefaultResidualEngineParameters
      : DefaultEngineParametersBase<instationary_,galerkin_>
    {

      static constexpr bool assembleVariablePart = true;

      static constexpr bool assembleConstantPart = true;

      static constexpr bool assembleOffDiagonalSkeletonPart = true;

      static constexpr bool assembleDiagonalSkeletonPart = true;

    };


    template<bool instationary_ = false, Galerkin galerkin_ = Galerkin::automatic>
    struct DefaultJacobianEngineParameters
      : DefaultEngineParametersBase<instationary_,galerkin_>
    {

      static constexpr bool assembleVariablePart = true;

      static constexpr bool assembleConstantPart = false;

      static constexpr bool assembleOffDiagonalSkeletonPart = true;

      static constexpr bool assembleDiagonalSkeletonPart = true;

    };

    template<bool instationary_ = false, Galerkin galerkin_ = Galerkin::automatic>
    struct DefaultApplyJacobianEngineParameters
      : DefaultEngineParametersBase<instationary_,galerkin_>
    {

      static constexpr bool assembleVariablePart = true;

      static constexpr bool assembleConstantPart = false;

      static constexpr bool assembleOffDiagonalSkeletonPart = true;

      static constexpr bool assembleDiagonalSkeletonPart = true;

    };

    template<Galerkin galerkin_ = Galerkin::automatic>
    struct DefaultPatternEngineParameters
      : DefaultEngineParametersBase<false,galerkin_>
    {};


    template<bool instationary_ = false, Galerkin galerkin_ = Galerkin::automatic>
    struct DefaultGridOperatorParameters
      : DefaultEngineParametersBase<instationary_,galerkin_>
    {

      using ResidualEngineParameters      = DefaultResidualEngineParameters<instationary_,galerkin_>;
      using JacobianEngineParameters      = DefaultJacobianEngineParameters<instationary_,galerkin_>;
      using ApplyJacobianEngineParameters = DefaultApplyJacobianEngineParameters<instationary_,galerkin_>;
      using PatternEngineParameters       = DefaultPatternEngineParameters<galerkin_>;

    };




    template<bool instationary_ = false, Galerkin galerkin_ = Galerkin::automatic>
    struct FastDGResidualEngineParameters
      : DefaultResidualEngineParameters<instationary_,galerkin_>
    {
      static constexpr bool fastDG = true;
    };

    template<bool instationary_ = false, Galerkin galerkin_ = Galerkin::automatic>
    struct FastDGJacobianEngineParameters
      : DefaultJacobianEngineParameters<instationary_,galerkin_>
    {
      static constexpr bool fastDG = true;
    };

    template<bool instationary_ = false, Galerkin galerkin_ = Galerkin::automatic>
    struct FastDGApplyJacobianEngineParameters
      : DefaultApplyJacobianEngineParameters<instationary_,galerkin_>
    {
      static constexpr bool fastDG = true;
    };

    template<Galerkin galerkin_ = Galerkin::automatic>
    struct FastDGPatternEngineParameters
      : DefaultPatternEngineParameters<galerkin_>
    {
      static constexpr bool fastDG = true;
    };

    template<bool instationary_ = false, Galerkin galerkin_ = Galerkin::automatic>
    struct FastDGGridOperatorParameters
      : DefaultEngineParametersBase<instationary_,galerkin_>
    {

      static constexpr bool fastDG = true;

      using ResidualEngineParameters      = FastDGResidualEngineParameters<instationary_,galerkin_>;
      using JacobianEngineParameters      = FastDGJacobianEngineParameters<instationary_,galerkin_>;
      using ApplyJacobianEngineParameters = FastDGApplyJacobianEngineParameters<instationary_,galerkin_>;
      using PatternEngineParameters       = FastDGPatternEngineParameters<galerkin_>;

    };

    template<typename TrialSpace_,
             typename TestSpace_,
             typename TrialConstraints_,
             typename TestConstraints_,
             bool enable_flavors>
    struct LocalFunctionSpaceTypes;

    template<
      typename TrialSpace_,
      typename TestSpace_,
      typename TrialConstraints_,
      typename TestConstraints_
      >
    struct LocalFunctionSpaceTypes<TrialSpace_,TestSpace_,TrialConstraints_,TestConstraints_,true>
    {

      using TestSpace        = TestSpace_;
      using TestConstraints  = TestConstraints_;

      template<typename Context, typename Flavor>
      using TestLocalSpace   = LocalFunctionSpace<TestSpace,Flavor>;

      template<typename Context, typename Flavor>
      using TestSpaceCache   = LFSIndexCache<TestLocalSpace<Context,Flavor>,TestConstraints,Context::fastDG()>;


      using TrialSpace       = TrialSpace_;
      using TrialConstraints = TrialConstraints_;

      template<typename Context, typename Flavor>
      using TrialLocalSpace = LocalFunctionSpace<TrialSpace,Flavor>;

      template<typename Context, typename Flavor>
      using TrialSpaceCache = LFSIndexCache<TrialLocalSpace<Context,Flavor>,TrialConstraints,Context::fastDG()>;

    };

    template<
      typename TrialSpace_,
      typename TestSpace_,
      typename TrialConstraints_,
      typename TestConstraints_
      >
    struct LocalFunctionSpaceTypes<TrialSpace_,TestSpace_,TrialConstraints_,TestConstraints_,false>
    {

      using TestSpace        = TestSpace_;
      using TestConstraints  = TestConstraints_;

      template<typename Context, typename>
      using TestLocalSpace   = LocalFunctionSpace<TestSpace,Flavor::Generic>;

      template<typename Context, typename>
      using TestSpaceCache   = LFSIndexCache<TestLocalSpace<Context,Flavor::Generic>,TestConstraints,Context::fastDG()>;


      using TrialSpace       = TrialSpace_;
      using TrialConstraints = TrialConstraints_;

      template<typename Context, typename>
      using TrialLocalSpace  = LocalFunctionSpace<TrialSpace,Flavor::Generic>;

      template<typename Context, typename>
      using TrialSpaceCache  = LFSIndexCache<TrialLocalSpace<Context,Flavor::Generic>,TrialConstraints,Context::fastDG()>;

    };



    template<
      typename TrialSpace_,
      typename TestSpace_,
      typename TrialConstraints_,
      typename TestConstraints_,
      bool enable_flavors
      >
    class FunctionSpaceProvider {

    public:

      using TrialConstraints = TrialConstraints_;
      using TestConstraints  = TestConstraints_;

    protected:

      using Types = LocalFunctionSpaceTypes<
        TrialSpace_,
        TestSpace_,
        TrialConstraints_,
        TestConstraints_,
        enable_flavors
        >;

    public:

      template<typename Context>
      using TestLocalSpace   = typename Types::template TestLocalSpace<Context,
                                                                       std::conditional_t<
                                                                         enable_flavors,
                                                                         typename Context::Flavor::Test,
                                                                         Flavor::Generic
                                                                         >
                                                                       >;

      template<typename Context>
      using TestSpaceCache   = typename Types::template TestSpaceCache<Context,
                                                                       std::conditional_t<
                                                                         enable_flavors,
                                                                         typename Context::Flavor::Test,
                                                                         Flavor::Generic
                                                                         >
                                                                       >;

      template<typename Context>
      using TrialLocalSpace  = typename Types::template TrialLocalSpace<Context,
                                                                        std::conditional_t<
                                                                          enable_flavors,
                                                                          typename Context::Flavor::Trial,
                                                                          Flavor::Generic
                                                                          >
                                                                        >;

      template<typename Context>
      using TrialSpaceCache  = typename Types::template TrialSpaceCache<Context,
                                                                        std::conditional_t<
                                                                          enable_flavors,
                                                                          typename Context::Flavor::Trial,
                                                                          Flavor::Generic
                                                                          >
                                                                        >;

      static constexpr bool unconstrained()
      {
        return std::is_same_v<TestConstraints,EmptyTransformation>;
      }

      static constexpr bool constrained()
      {
        return not unconstrained();
      }

      template<typename Context>
      TestSpaceCache<Context> makeTestSpaceCache(const Context&) const
      {
        return TestSpaceCache<Context>(testConstraints());
      }

      template<typename Context>
      TrialSpaceCache<Context> makeTrialSpaceCache(const Context&) const
      {
        return TrialSpaceCache<Context>(trialConstraints());
      }

      const TrialConstraints& trialConstraints() const
      {
        if constexpr (constrained())
          return *_trial_constraints;
        else
          return _trial_constraints;
      }

      const TestConstraints& testConstraints() const
      {
        if constexpr (constrained())
          return *_test_constraints;
        else
          return _test_constraints;
      }

    protected:

      template<int i = 0>
      FunctionSpaceProvider(
        const TrialConstraints* trial_constraints,
        const TestConstraints* test_constraints,
        std::enable_if_t<i == 0 and constrained(),int> = 0
        )
        : _trial_constraints(trial_constraints)
        , _test_constraints(test_constraints)
      {}

      template<int i = 0>
      FunctionSpaceProvider(
        const TrialConstraints* trial_constraints,
        const TestConstraints* test_constraints,
        std::enable_if_t<i == 0 and unconstrained(),int> = 0
        )
      {}

      template<int i = 0>
      FunctionSpaceProvider(
        std::enable_if_t<i == 0 and unconstrained(),int> = 0
        )
      {}

    private:

      std::conditional_t<
        std::is_same_v<TrialConstraints,EmptyTransformation>,
        TrialConstraints,
        const TrialConstraints*
        > _trial_constraints;

      std::conditional_t<
        std::is_same_v<TestConstraints,EmptyTransformation>,
        TestConstraints,
        const TestConstraints*
        > _test_constraints;


    };




    template<typename Context>
    struct TimeData
      : public Context
    {
      using TimeReal = typename Context::Engine::TimeReal;
      using Context_ = Context;

      static constexpr bool instationary()
      {
        return Context_::Engine::instationary();
      }

      TimeReal time() const
      {
        return Context_::engine().time();
      }

      TimeReal dt() const
      {
        return Context_::engine().dt();
      }

      TimeReal t0() const
      {
        return Context_::engine().t0();
      }

      TimeData(Context&& ctx)
        : Context(std::move(ctx))
      {}

    };

    template<typename Context>
    auto timeData(Context&& context)
    {
      return TimeData<Context>(std::move(context));
    }


    template<typename Real, bool instationary_>
    class InstationaryEngineBase
    {

    public:

      using TimeReal = Real;

      static constexpr bool instationary()
      {
        return instationary_;
      }

      TimeReal dt() const
      {
        return _dt;
      }

      TimeReal time() const
      {
        return _time;
      }

      TimeReal t0() const
      {
        return _t0;
      }

      void setOneStepMethod(std::shared_ptr<const OneStep::Method<TimeReal>> one_step_method)
      {
        _one_step_method = one_step_method;
      }

      const OneStep::Method<TimeReal>& oneStepMethod() const
      {
        assert(_one_step_method);
        return *_one_step_method;
      }

      int startStep(TimeReal t0, TimeReal dt)
      {
        _stage = 0;
        _time = _t0 = t0;
        _dt = dt;

        _timestep_factor = 1.0;
        _timestep_time_factor = 1.0;
        switch (_scaling)
        {
        case DTScaling::divide:
          _timestep_time_factor /= _dt;
          break;
        case DTScaling::multiply:
          _timestep_factor = dt;
          break;
        default:
          DUNE_THROW(AssemblyError,"Unknown time step scaling method");
        }
        return oneStepMethod().stages();
      }

      Real weight() const
      {
        return _weight;
      }

      Real timeWeight() const
      {
        return _time_weight;
      }

      Real timestepFactor() const
      {
        return _timestep_factor;
      }

      Real timestepTimeFactor() const
      {
        return _timestep_time_factor;
      }

      int stage() const
      {
        return _stage;
      }

      void setStage(int stage)
      {
        _stage = stage;
      }

      template<typename Assembler, typename TrialVector>
      bool acceptStage(int stage, Assembler&, const TrialVector&)
      {
        if (stage == _stage)
          return false;
        assert(stage == _stage - 1);
        ++_stage;
        updateWeights();
        return true;;
      }

      void updateWeights()
      {
        if (instationary())
        {
          _weight = oneStepMethod().weight(_stage,_stage) * _timestep_factor;
          _time_weight = oneStepMethod().timeDerivativeWeight(_stage,_stage) * _timestep_time_factor;
        }
        else
        {
          _weight = 1.0;
          _time_weight = 0.0;
        }
      }

    protected:

      void setWeight(Real weight)
      {
        _weight = weight;
      }

      void setTimeWeight(Real time_weight)
      {
        _time_weight = time_weight;
      }

    private:

      std::shared_ptr<const OneStep::Method<Real>> _one_step_method = nullptr;
      int _stage = 0;
      bool _stage_accept_mode = false;
      DTScaling _scaling = DTScaling::divide;
      TimeReal _time = 0.0;
      TimeReal _t0 = 0.0;
      TimeReal _dt = 0.0;
      Real _weight = 1.0;
      Real _time_weight = 1.0;
      TimeReal _timestep_factor = 1.0;
      TimeReal _timestep_time_factor = 1.0;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_ENGINEBASE_HH
