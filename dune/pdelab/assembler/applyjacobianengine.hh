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

namespace Dune::PDELab::Experimental {

  template<
    typename TrialVector_,
    typename TestVector_,
    typename LOP,
    typename TrialConstraints_ = EmptyTransformation,
    typename TestConstraints_ = EmptyTransformation,
    typename EngineParameters = DefaultApplyJacobianEngineParameters<false,Galerkin::automatic>
    >
  class ApplyJacobianEngine
    : public InstationaryEngineBase<typename TrialVector_::value_type,EngineParameters::instationary>
    , public FunctionSpaceProvider<typename TrialVector_::GridFunctionSpace,
                                   typename TestVector_::GridFunctionSpace,
                                   TrialConstraints_,
                                   TestConstraints_,
                                   not LocalOperator::disableFunctionSpaceFlavors<LOP>()
                                   >
  {

    static constexpr bool enable_flavors = not LocalOperator::disableFunctionSpaceFlavors<LOP>();

    using IEB = InstationaryEngineBase<typename TrialVector_::value_type,EngineParameters::instationary>;

    using FSP = FunctionSpaceProvider<
      typename TrialVector_::GridFunctionSpace,
      typename TestVector_::GridFunctionSpace,
      TrialConstraints_,
      TestConstraints_,
      enable_flavors
      >;

    using Types = typename FSP::Types;

  public:

    using IEB::instationary;

    using FSP::unconstrained;
    using FSP::trialConstraints;
    using FSP::testConstraints;

    using size_type        = std::size_t;

    using TrialVector      = TrialVector_;
    using TestVector       = TestVector_;
    using TrialConstraints = typename Types::TrialConstraints;
    using TestConstraints  = typename Types::TestConstraints;

    using TestSpace        = typename Types::TestSpace;

    using TrialSpace       = typename Types::TrialSpace;

    using EntitySet        = typename TestSpace::Traits::EntitySet;

    static constexpr
    std::bool_constant<EngineParameters::template galerkin<TrialSpace,TestSpace>>
    isGalerkin()
    {
      return {};
    }

  private:

    LOP* _lop;

    const TrialVector* _coefficient;
    const TrialVector* _linearization_point;
    TestVector* _result;

    EmptyTransformation _empty_constraints;

    bool _symmetric_dirichlet_constraints;

  public:

    template<typename CellFlavor>
    struct Data
      : public Context::RootContext
    {

      using Flavor = CellFlavor;
      using Engine = ApplyJacobianEngine;
      using EntitySet = typename Engine::EntitySet;

      static constexpr bool assembleVariablePart()
      {
        return EngineParameters::assembleVariablePart;
      }

      static constexpr bool assembleConstantPart()
      {
        return EngineParameters::assembleConstantPart;
      }

      static constexpr bool assembleOffDiagonalSkeletonPart()
      {
        return EngineParameters::assembleOffDiagonalSkeletonPart;
      }

      static constexpr bool assembleDiagonalSkeletonPart()
      {
        return EngineParameters::assembleDiagonalSkeletonPart;
      }

      static constexpr auto isGalerkin()
      {
        return ApplyJacobianEngine::isGalerkin();
      }

      static constexpr std::bool_constant<EngineParameters::fastDG> fastDG()
      {
        return {};
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
      const TrialVector& coefficient,
      const TrialVector& linearization_point,
      TestVector& result,
      LOP& lop,
      const TrialConstraints& trial_constraints,
      const TestConstraints& test_constraints,
      EngineParameters = {}
      )
      : FSP(&trial_constraints,&test_constraints)
      , _lop(&lop)
      , _coefficient(&coefficient)
      , _linearization_point(&linearization_point)
      , _result(&result)
      , _symmetric_dirichlet_constraints(false)
    {}

    template<typename LOP_>
    ApplyJacobianEngine(
      const TrialVector& coefficient,
      const TrialVector& linearization_point,
      TestVector& result,
      LOP_& lop,
      std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,EngineParameters> = {}
      )
      : _lop(&lop)
      , _coefficient(&coefficient)
      , _linearization_point(&linearization_point)
      , _result(&result)
      , _symmetric_dirichlet_constraints(false)
    {}


    ApplyJacobianEngine(
      LOP& lop,
      const TrialConstraints& trial_constraints,
      const TestConstraints& test_constraints,
      EngineParameters = {}
      )
      : FSP(&trial_constraints,&test_constraints)
      , _lop(&lop)
      , _coefficient(nullptr)
      , _linearization_point(nullptr)
      , _result(nullptr)
      , _symmetric_dirichlet_constraints(false)
    {}

    template<typename LOP_>
    ApplyJacobianEngine(
      LOP_& lop,
      std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,EngineParameters> = {}
      )
      : _lop(&lop)
      , _coefficient(nullptr)
      , _linearization_point(nullptr)
      , _result(nullptr)
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
      return isNonlinear(localOperator());
    }

    const TrialVector& coefficient() const
    {
      return *_coefficient;
    }

    void setCoefficient(const TrialVector& coefficient)
    {
      _coefficient = &coefficient;
    }

    const TrialVector& linearizationPoint() const
    {
#if DUNE_PDELAB_ENABLE_CHECK_ASSEMBLY
      if (not isNonlinear(localOperator()))
        DUNE_THROW(AssemblyError, "Not allowed to call linearizationPoint() for linear operators");
#endif
      return *_linearization_point;
    }

    void setLinearizationPoint(const TrialVector& linearization_point)
    {
#if DUNE_PDELAB_ENABLE_CHECK_ASSEMBLY
      if (not isNonlinear(localOperator()))
        DUNE_THROW(AssemblyError, "Not allowed to call setLinearizationPoint() for linear operators");
#endif
      _linearization_point = &linearization_point;
    }

    const TestSpace& testSpace() const
    {
      return _result->gridFunctionSpace();
    }

    const TrialSpace& trialSpace() const
    {
      return _coefficient->gridFunctionSpace();
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
                      vectorData<TestVector,Flavor::Test,LocalViewDataMode::accumulate>(
                        cellLinearizationPointData(
                          models<Concept::PossiblyNonLinear,LOP>(),
                          vectorData<TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                            cellCoefficientData(
                              vectorData<TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                                trialSpaceData(
                                  testSpaceData(
                                    cellGridData(
                                      timeData(
                                        Data<CellFlavor::Outside<enable_flavors>>(*this)
                                        )))))))))))),
                insideCell(
                  extractCellContext(
                    *_lop,
                    cellTimeResultData(
                      std::bool_constant<instationary()>{},
                      cellResultData(
                        vectorData<TestVector,Flavor::Test,LocalViewDataMode::accumulate>(
                          cellLinearizationPointData(
                            models<Concept::PossiblyNonLinear,LOP>(),
                            vectorData<TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                              cellCoefficientData(
                                vectorData<TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                                  trialSpaceData(
                                    testSpaceData(
                                      cellGridData(
                                        timeData(
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
      DefaultApplyJacobianEngineParameters<false,Galerkin::automatic>
      >;

  template<typename TrialVector, typename TestVector, typename LOP, typename EngineParameters>
  ApplyJacobianEngine(
      const TrialVector&,
      const TrialVector&,
      TestVector&,
      LOP&,
      EngineParameters
    )
    -> ApplyJacobianEngine<
      TrialVector,
      TestVector,
      LOP,
      EmptyTransformation,
      EmptyTransformation,
      EngineParameters
      >;

} // namespace Dune::PDELab::Experimental

#endif // DUNE_PDELAB_ASSEMBLER_APPLYJACOBIANENGINE_HH
