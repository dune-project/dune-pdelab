// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_APPLYJACOBIANENGINE_HH
#define DUNE_PDELAB_ASSEMBLER_APPLYJACOBIANENGINE_HH

#include <dune/pdelab/common/checks.hh>
#include <dune/pdelab/common/exceptions.hh>
#include <dune/pdelab/common/typetraits.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/common/uncachedvectorview.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/assembler/context.hh>
#include <dune/pdelab/assembler/celldata.hh>
#include <dune/pdelab/assembler/functionspacedata.hh>
#include <dune/pdelab/assembler/vectordata.hh>
#include <dune/pdelab/assembler/utility.hh>
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
    class ApplyJacobianEngine
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

      const TrialVector* _argument;
      const TrialVector* _linearization_point;
      TestVector* _result;

      EmptyTransformation _empty_constraints;

      TrialConstraints* _trial_constraints;
      TestConstraints* _test_constraints;

      bool _symmetric_dirichlet_constraints;

    public:

      template<typename CellFlavor>
      struct Data
        : public Context::RootContext
      {

        using Flavor = CellFlavor;
        using Engine = ApplyJacobianEngine;
        using EntitySet = typename Engine::EntitySet;

        static constexpr bool skipVariablePart()
        {
          return false;
        }

        static constexpr bool skipConstantPart()
        {
          return true;
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
          return ApplyJacobianEngine::isGalerkin();
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
        return std::is_invocable_v<decltype(LocalOperator::boundaryJacobian()),LOP,PlaceHolder&>;
      }

      static constexpr bool visitSkeletonIntersections()
      {
        return std::is_invocable_v<decltype(LocalOperator::skeletonJacobian()),LOP,PlaceHolder&>;
      }

      static constexpr bool visitPeriodicIntersections()
      {
        return visitSkeletonIntersections();
      }

      static constexpr bool visitProcessorIntersections()
      {
        return visitSkeletonIntersections();
      }

      ApplyJacobianEngine(
        const TrialVector& argument,
        const TrialVector& linearization_point,
        TestVector& result,
        LOP& lop,
        const TrialConstraints& trial_constraints,
        const TestConstraints& test_constraints,
        std::integral_constant<Galerkin,galerkin> = std::integral_constant<Galerkin,galerkin>{}
        )
        : _lop(&lop)
        , _argument(&argument)
        , _linearization_point(&linearization_point)
        , _result(&result)
        , _trial_constraints(&trial_constraints)
        , _test_constraints(&test_constraints)
        , _symmetric_dirichlet_constraints(false)
      {}

      template<typename LOP_>
      ApplyJacobianEngine(
        const TrialVector& argument,
        const TrialVector& linearization_point,
        TestVector& result,
        LOP_& lop,
        std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,std::integral_constant<Galerkin,galerkin>> = std::integral_constant<Galerkin,galerkin>{}
        )
        : _lop(&lop)
        , _argument(&argument)
        , _linearization_point(&linearization_point)
        , _result(&result)
        , _trial_constraints(&_empty_constraints)
        , _test_constraints(&_empty_constraints)
        , _symmetric_dirichlet_constraints(false)
      {}


      TestVector& result()
      {
        return *_result;
      }

      void setResult(TestVector& result)
      {
        _result = &result;
      }

      constexpr bool bindLinearizationPoint() const
      {
        return isNonLinear(localOperator());
      }

      const TrialVector& argument() const
      {
        return *_argument;
      }

      void setArgument(const TrialVector& argument)
      {
        _argument = &argument;
      }

      const TrialVector& linearizationPoint() const
      {
#if DUNE_PDELAB_ENABLE_CHECK_ASSEMBLY
        if (not isNonLinear(localOperator()))
          DUNE_THROW(AssemblyError, "Not allowed to call linearizationPoint() for linear operators");
#endif
        return *_linearization_point;
      }

      void setLinearizationPoint(const TrialVector& linearization_point)
      {
#if DUNE_PDELAB_ENABLE_CHECK_ASSEMBLY
        if (not isNonLinear(localOperator()))
          DUNE_THROW(AssemblyError, "Not allowed to call setLinearizationPoint() for linear operators");
#endif
        _linearization_point = &linearization_point;
      }

      const TestSpace& testSpace() const
      {
        return _result->gridFunctionSpace();
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

      const TestConstraints& testConstraints() const
      {
        return *_test_constraints;
      }

      const TrialSpace& trialSpace() const
      {
        return _argument->gridFunctionSpace();
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

      const TrialConstraints& trialConstraints() const
      {
        return *_trial_constraints;
      }

      LOP& localOperator()
      {
        return *_lop;
      }

      const LOP& localOperator() const
      {
        return *_lop;
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
                    cellTimeResultData(
                      std::bool_constant<instationary()>{},
                      cellResultData(
                        cachedVectorData<UncachedVectorView,TestVector,Flavor::Test,LocalViewDataMode::accumulate>(
                          cellLinearizationPointData(
                            models<Concept::PossiblyNonLinear,LOP>(),
                            cachedVectorData<ConstUncachedVectorView,TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                              cellArgumentData(
                                cachedVectorData<ConstUncachedVectorView,TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                                  trialSpaceData(
                                    testSpaceData(
                                      cellGridData(
                                        timeData(
                                          std::bool_constant<instationary()>{},
                                          Data<CellFlavor::Outside<enable_flavors>>(*this)
                                          )))))))))))),
                  insideCell(
                    extractCellContext(
                      *_lop,
                      cellTimeResultData(
                        std::bool_constant<instationary()>{},
                        cellResultData(
                          cachedVectorData<UncachedVectorView,TestVector,Flavor::Test,LocalViewDataMode::accumulate>(
                            cellLinearizationPointData(
                              models<Concept::PossiblyNonLinear,LOP>(),
                              cachedVectorData<ConstUncachedVectorView,TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                                cellArgumentData(
                                  cachedVectorData<ConstUncachedVectorView,TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                                    trialSpaceData(
                                      testSpaceData(
                                        cellGridData(
                                          timeData(
                                            std::bool_constant<instationary()>{},
                                            Data<CellFlavor::Inside<enable_flavors>>(*this)
                                            )))))))))))))))));
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
        invoke_if_possible(LocalOperator::volumeApplyJacobian(),*_lop,ctx.cellContext());
      }

      template<typename Context>
      void startIntersections(Context& ctx)
      {
        invoke_if_possible(LocalOperator::startIntersections(),*_lop,ctx.cellContext());
      }

      template<typename Context>
      void skeleton(Context& ctx)
      {
        invoke_if_possible(LocalOperator::skeletonApplyJacobian(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void periodic(Context& ctx)
      {
        invoke_if_possible(LocalOperator::skeletonApplyJacobian(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void boundary(Context& ctx)
      {
        invoke_if_possible(LocalOperator::boundaryApplyJacobian(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void processor(Context& ctx)
      {
        invoke_if_possible(LocalOperator::boundaryApplyJacobian(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void finishIntersections(Context& ctx)
      {
        invoke_if_possible(LocalOperator::finishIntersections(),*_lop,ctx.cellContext());
      }

      template<typename Context>
      void volumePostIntersections(Context& ctx)
      {
        invoke_if_possible(LocalOperator::volumeJacobianPostIntersections(),*_lop,ctx.cellContext());
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
        constrain_residual(trialConstraints(),result());
      }

      template<typename Context>
      void result(Context& ctx)
      {}

    };


    template<typename TrialVector, typename TestVector, typename LOP>
    ApplyJacobianEngine(
        const TrialVector&,
        const TrialVector&,
        TestVector&,
        LOP&
      )
      -> ApplyJacobianEngine<
        TrialVector,
        TestVector,
        LOP,
        EmptyTransformation,
        EmptyTransformation,
        false,
        Galerkin::automatic
        >;

    template<typename TrialVector, typename TestVector, typename LOP, Galerkin galerkin>
    ApplyJacobianEngine(
        const TrialVector&,
        const TrialVector&,
        TestVector&,
        LOP&,
        std::integral_constant<Galerkin,galerkin>
      )
      -> ApplyJacobianEngine<
        TrialVector,
        TestVector,
        LOP,
        EmptyTransformation,
        EmptyTransformation,
        false,
        galerkin
        >;


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_APPLYJACOBIANENGINE_HH
