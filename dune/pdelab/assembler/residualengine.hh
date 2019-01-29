// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_RESIDUALENGINE_HH
#define DUNE_PDELAB_ASSEMBLER_RESIDUALENGINE_HH

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/assembler/context.hh>
#include <dune/pdelab/assembler/celldata.hh>
#include <dune/pdelab/assembler/functionspacedata.hh>
#include <dune/pdelab/assembler/vectordata.hh>
#include <dune/pdelab/localoperator/guardedcalls.hh>
#include <dune/pdelab/gridfunctionspace/flavor.hh>
#include <dune/pdelab/assembler/enginebase.hh>

namespace Dune {
  namespace PDELab {

    template<
      typename TrialVector_,
      typename TestVector_,
      typename LOP,
      typename TrialConstraints_ = EmptyTransformation,
      typename TestConstraints_ = EmptyTransformation,
      bool instationary_ = false,
      Galerkin galerkin = Galerkin::automatic
      >
    class ResidualEngine
        : public InstationaryEngineBase<typename TrialVector_::value_type,instationary_>
    {

      static constexpr bool enable_flavors = not LocalOperator::disableFunctionSpaceFlavors<LOP>();

      using Types = LocalFunctionSpaceTypes<
        typename TestVector_::GridFunctionSpace,
        typename TrialVector_::GridFunctionSpace,
        TrialConstraints_,
        TestConstraints_,
        enable_flavors
        >;

      using IEB = InstationaryEngineBase<typename TrialVector_::value_type,instationary_>;

    public:

      using IEB::instationary;
      using IEB::stage;
      using IEB::setStage;
      using IEB::oneStepMethod;
      using IEB::setOneStepMethod;
      using IEB::timestepFactor;
      using IEB::timestepTimeFactor;

      using size_type        = std::size_t;

      using TrialVector      = TrialVector_;
      using TestVector       = TestVector_;
      using TrialConstraints = TrialConstraints_;
      using TestConstraints  = TestConstraints_;

      using TestSpace        = typename Types::TestSpace;

      template<typename Flavor = Flavor::Generic>
      using TestLocalSpace   = typename Types::template TestLocalSpace<Flavor>;

      template<typename Flavor = Flavor::Generic>
      using TestSpaceCache   = typename Types::template TestSpaceCache<Flavor>;


      using TrialSpace       = typename Types::TrialSpace;

      template<typename Flavor = Flavor::Generic>
      using TrialLocalSpace  = typename Types::template TrialLocalSpace<Flavor>;

      template<typename Flavor = Flavor::Generic>
      using TrialSpaceCache  = typename Types::template TrialSpaceCache<Flavor>;


      using EntitySet        = typename TestSpace::Traits::EntitySet;

      using TimeReal         = typename TrialVector::value_type;

      static constexpr bool unconstrained()
      {
        return std::is_same<TestConstraints,EmptyTransformation>::value;
      }

      static constexpr
      std::bool_constant<
        galerkin == Galerkin::automatic ?
        std::is_same<TrialSpace,TestSpace>::value
        : bool(galerkin)
        >
      isGalerkin()
      {
        return {};
      }

    private:

      LOP* _lop;

      bool _stage_accept_mode = false;

      std::vector<std::shared_ptr<TestVector>> _residuals;
      std::vector<std::shared_ptr<TestVector>> _time_residuals;

      const TrialVector* _argument;
      TestVector* _residual      = nullptr;;
      TestVector* _time_residual = nullptr;

      EmptyTransformation _empty_constraints;

      const TrialConstraints* _trial_constraints;
      const TestConstraints* _test_constraints;

    public:

      template<typename CellFlavor>
      struct Data
        : public Context::RootContext
      {

        using Flavor = CellFlavor;
        using Engine = ResidualEngine;
        using EntitySet = typename Engine::EntitySet;

        static constexpr bool skipVariablePart()
        {
          return false;
        }

        static constexpr bool skipConstantPart()
        {
          return false;
        }

        static constexpr bool skipOffDiagonalSkeletonPart()
        {
          return false;
        }

        static constexpr bool skipDiagonalSkeletonPart()
        {
          return true;
        }

        static constexpr auto isGalerkin()
        {
          return ResidualEngine::isGalerkin();
        }

        Data(Engine& engine)
          : _engine(engine)
        {}

        Engine& engine()
        {
          return _engine;
        }

        const Engine& engine() const
        {
          return _engine;
        }

      private:

        Engine& _engine;

      };


      static constexpr bool intersectionsTwoSided()
      {
        return LocalOperator::intersectionsTwoSided<LOP>();
      }

      static constexpr bool visitBoundaryIntersections()
      {
        return std::is_invocable_v<decltype(LocalOperator::boundaryIntegral()),LOP,PlaceHolder&>;
      }

      static constexpr bool visitSkeletonIntersections()
      {
        return std::is_invocable_v<decltype(LocalOperator::skeletonIntegral()),LOP,PlaceHolder&>;
      }

      static constexpr bool visitPeriodicIntersections()
      {
        return visitSkeletonIntersections();
      }

      static constexpr bool visitProcessorIntersections()
      {
        return visitSkeletonIntersections();
      }

      template<typename LOP_>
      ResidualEngine(
        const TrialVector& trial_vector,
        TestVector& test_vector,
        LOP_& lop,
        std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,std::integral_constant<Galerkin,galerkin>> = std::integral_constant<Galerkin,galerkin>{}
        )
        : _lop(&lop)
        , _argument(&trial_vector)
        , _residual(&test_vector)
        , _trial_constraints(nullptr)
        , _test_constraints(nullptr)
      {}

      ResidualEngine(const TrialVector& trial_vector, TestVector& test_vector, LOP& lop,
                     const TrialConstraints& trial_constraints, const TestConstraints& test_constraints,
                     std::integral_constant<Galerkin,galerkin> = std::integral_constant<Galerkin,galerkin>{})
        : _lop(&lop)
        , _argument(&trial_vector)
        , _residual(&test_vector)
        , _trial_constraints(&trial_constraints)
        , _test_constraints(&test_constraints)
      {}

      template<typename LOP_>
      ResidualEngine(
        LOP_& lop,
        std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,std::integral_constant<Galerkin,galerkin>> = std::integral_constant<Galerkin,galerkin>{}
        )
        : _lop(&lop)
        , _argument(nullptr)
        , _residual(nullptr)
        , _trial_constraints(nullptr)
        , _test_constraints(nullptr)
      {}

      ResidualEngine(
        LOP& lop,
        const TrialConstraints& trial_constraints,
        const TestConstraints& test_constraints,
        std::integral_constant<Galerkin,galerkin> = std::integral_constant<Galerkin,galerkin>{}
        )
        : _lop(&lop)
        , _argument(nullptr)
        , _residual(nullptr)
        , _trial_constraints(&trial_constraints)
        , _test_constraints(&test_constraints)
      {}

      void setArgument(const TrialVector& argument)
      {
        _argument = &argument;
      }

      void setResidual(TestVector& residual)
      {
        _residual = &residual;
      }

      TestVector& residual()
      {
        if (_stage_accept_mode)
        {
          assert(_residuals.size() > stage());
          assert(_residuals[stage()]);
          return *_residuals[stage()];
        }
        else
        {
          assert(_residual);
          return *_residual;
        }
      }

      template<typename h = int>
      TestVector& timeResidual()
      {
        static_assert(Std::to_true_v<h> and instationary(),"Calling timeResidual() is only allowed in instationary mode");
        if (_stage_accept_mode)
        {
          assert(_time_residuals.size() > stage());
          assert(_time_residuals[stage()]);
          return *_time_residuals[stage()];
        }
        else
        {
          assert(_time_residual);
          return *_time_residual;
        }
      }

      const TrialVector& argument() const
      {
        return *_argument;
      }

      const TestSpace& testSpace() const
      {
        return _residual->gridFunctionSpace();
      }

      const TestConstraints& testConstraints() const
      {
        return *_test_constraints;
      }

      template<typename Flavor_>
      std::enable_if_t<
        Std::to_true_type<Flavor_>::value and enable_flavors,
        TestSpaceCache<Flavor_>
        >
      makeTestSpaceCache(Flavor_) const
      {
        return TestSpaceCache<Flavor_>(testConstraints());
      }

      template<typename Flavor_>
      std::enable_if_t<
        Std::to_true_type<Flavor_>::value and not enable_flavors,
        TestSpaceCache<Flavor::Generic>
        >
      makeTestSpaceCache(Flavor_) const
      {
        return TestSpaceCache<Flavor::Generic>(testConstraints());
      }

      const TrialSpace& trialSpace() const
      {
        return _argument->gridFunctionSpace();
      }

      const TrialConstraints& trialConstraints() const
      {
        return *_trial_constraints;
      }

      template<typename Flavor_>
      std::enable_if_t<
        Std::to_true_type<Flavor_>::value and enable_flavors,
        TrialSpaceCache<Flavor_>
        >
      makeTrialSpaceCache(Flavor_) const
      {
        return TrialSpaceCache<Flavor_>(trialConstraints());
      }

      template<typename Flavor_>
      std::enable_if_t<
        Std::to_true_type<Flavor_>::value and not enable_flavors,
        TrialSpaceCache<Flavor::Generic>
        >
      makeTrialSpaceCache(Flavor_) const
      {
        return TrialSpaceCache<Flavor::Generic>(trialConstraints());
      }

      LOP& localOperator()
      {
        return *_lop;
      }

      const LOP& localOperator() const
      {
        return *_lop;
      }

      void setOneStepMethod(shared_ptr<const OneStep::Method<TimeReal>> one_step_method)
      {
        IEB::setOneStepMethod(one_step_method);
        _residuals.resize(oneStepMethod().stages());
        _time_residuals.resize(oneStepMethod().stages());
        for (auto& r : _residuals)
          r = std::make_shared<TestVector>(testSpace());
        for (auto& r : _time_residuals)
          r = std::make_shared<TestVector>(testSpace());
      }

      template<typename Assembler>
      bool acceptStage(int new_stage, Assembler& assembler, const TrialVector& solution)
      {
        if (new_stage == stage())
          return false;
        _stage_accept_mode = true;
        updateWeights();
        auto argument = _argument;
        _argument = &solution;
        auto residual = _residual;
        auto time_residual = _time_residual;
        assembler.assemble(*this);
        _residual = residual;
        _time_residual = time_residual;
        _stage_accept_mode = false;
        IEB::acceptStage(new_stage,assembler,solution);
        if (stage() == oneStepMethod().stages())
          _argument = argument;
        return true;
      }

      void updateWeights()
      {
        IEB::updateWeights();
        if (instationary() and _stage_accept_mode)
        {
          IEB::setWeight(IEB::timestepFactor());
          IEB::setTimeWeight(IEB::timestepTimeFactor());
        }
      }

      template<typename Assembler>
      auto context(const Assembler& assembler)
      {
        return
          extractContext(
            *_lop,
            intersectionDomainData(
              cellDomainData(
                outsideCell(
                  extractCellContext(
                    *_lop,
                    cellTimeResidualData<UncachedVectorView,TestVector,Flavor::Test>(
                      std::bool_constant<instationary()>{},
                      cellResidualData(
                        cachedVectorData<UncachedVectorView,TestVector,Flavor::Test,LocalViewDataMode::accumulate>(
                          cellArgumentData(
                            cachedVectorData<ConstUncachedVectorView,TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                              trialSpaceData(
                                testSpaceData(
                                  cellGridData(
                                    timeData(
                                      Data<CellFlavor::Outside<enable_flavors>>(*this)
                                      )))))))))),
                  insideCell(
                    extractCellContext(
                      *_lop,
                      cellTimeResidualData<UncachedVectorView,TestVector,Flavor::Test>(
                        std::bool_constant<instationary()>{},
                        cellResidualData(
                          cachedVectorData<UncachedVectorView,TestVector,Flavor::Test,LocalViewDataMode::accumulate>(
                            cellArgumentData(
                              cachedVectorData<ConstUncachedVectorView,TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                                trialSpaceData(
                                  testSpaceData(
                                    cellGridData(
                                      timeData(
                                        Data<CellFlavor::Inside<enable_flavors>>(*this)
                                        )))))))))))))));
      }

      template<typename Context>
      void start(Context& ctx)
      {
        invoke_if_possible(LocalOperator::start(),*_lop,ctx);
      }

      template<typename Context, typename Element, typename Index>
      bool skipCell(Context& ctx, const Element& element, Index index) const
      {
        return invoke_or(LocalOperator::skipCell(),false,*_lop,ctx,element,index);
      }

      template<typename Context, typename Element, typename Index>
      void startCell(Context& ctx, const Element& element, Index index) const
      {
        invoke_if_possible(LocalOperator::startCell(),*_lop,ctx,element,index);
      }

      template<typename Context>
      void volume(Context& ctx)
      {
        invoke_if_possible(LocalOperator::volumeIntegral(),*_lop,ctx.cellContext());
      }

      template<typename Context>
      void startIntersections(Context& ctx)
      {
        invoke_if_possible(LocalOperator::startIntersections(),*_lop,ctx.cellContext());
      }

      template<typename Context>
      void skeleton(Context& ctx)
      {
        invoke_if_possible(LocalOperator::skeletonIntegral(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void periodic(Context& ctx)
      {
        invoke_if_possible(LocalOperator::skeletonIntegral(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void boundary(Context& ctx)
      {
        invoke_if_possible(LocalOperator::boundaryIntegral(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void processor(Context& ctx)
      {
        invoke_if_possible(LocalOperator::boundaryIntegral(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void finishIntersections(Context& ctx)
      {
        invoke_if_possible(LocalOperator::finishIntersections(),*_lop,ctx.cellContext());
      }

      template<typename Context>
      void volumePostIntersections(Context& ctx)
      {
        invoke_if_possible(LocalOperator::volumeIntegralPostIntersections(),*_lop,ctx.cellContext());
      }

      template<typename Context, typename Element, typename Index>
      void finishCell(Context& ctx, const Element& element, Index index) const
      {
        invoke_if_possible(LocalOperator::finishCell(),*_lop,ctx,element,index);
      }

      template<typename Context>
      void finish(Context& ctx)
      {
        invoke_if_possible(LocalOperator::finish(),*_lop,ctx);
        if (instationary() and not _stage_accept_mode)
        {
          const auto& osm = oneStepMethod();
          for (int r = 0; r < stage(); ++r) {
            if (osm.timeDerivativeActive(stage(),r))
              _time_residual->axpy(osm.timeDerivativeWeight(stage(),r)*timestepTimeFactor(),*_time_residuals[r]);
            if (osm.active(stage(),r))
              _residual->axpy(osm.weight(stage(),r)*timestepFactor(),*_residuals[r]);
          }
        }
        constrain_residual(testConstraints(),residual());
        if constexpr(instationary())
          if (_residual != _time_residual)
            constrain_residual(testConstraints(),timeResidual());
      }

      template<typename Context>
      void result(Context& ctx)
      {}

    };

    template<typename Argument, typename Residual, typename LOP>
    ResidualEngine(
        const Argument&,
        Residual&,
        LOP&
      )
      -> ResidualEngine<
        Argument,
        Residual,
        LOP,
        EmptyTransformation,
        EmptyTransformation,
        false,
        Galerkin::automatic
        >;


    template<typename Argument, typename Residual, typename LOP, Galerkin galerkin>
    ResidualEngine(
        const Argument&,
        Residual&,
        LOP&,
        std::integral_constant<Galerkin,galerkin>
      )
      -> ResidualEngine<
        Argument,
        Residual,
        LOP,
        EmptyTransformation,
        EmptyTransformation,
        false,
        galerkin
        >;


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_RESIDUALENGINE_HH
