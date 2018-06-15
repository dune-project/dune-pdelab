// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_RESIDUALENGINE_HH
#define DUNE_PDELAB_ASSEMBLER_RESIDUALENGINE_HH

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/assembler/context.hh>
#include <dune/pdelab/assembler/celldata.hh>
#include <dune/pdelab/assembler/functionspacedata.hh>
#include <dune/pdelab/assembler/vectordata.hh>
#include <dune/pdelab/localoperator/guardedcalls.hh>

namespace Dune {
  namespace PDELab {

    template<bool is_galerkin = true>
    struct ResidualEngineAssemblyFlags
    {
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

      static constexpr std::bool_constant<is_galerkin> isGalerkin()
      {
        return {};
      }
    };

    enum class Galerkin { disable, enable, automatic };

    static constexpr auto enableGalerkin = std::integral_constant<Galerkin,Galerkin::enable>{};
    static constexpr auto disableGalerkin = std::integral_constant<Galerkin,Galerkin::disable>{};
    static constexpr auto automaticGalerkin = std::integral_constant<Galerkin,Galerkin::automatic>{};

    template<typename TrialVector, typename TestVector, typename LOP, Galerkin galerkin = Galerkin::automatic>
    class ResidualEngine
    {

    public:

      using TestSpace       = typename TestVector::GridFunctionSpace;
      using TestLocalSpace  = LocalFunctionSpace<TestSpace>;
      using TestLFS         = LocalFunctionSpace<TestSpace>;
      using TestLFSCache    = LFSIndexCache<TestLFS>;

      using TrialSpace      = typename TrialVector::GridFunctionSpace;
      using TrialLocalSpace = LocalFunctionSpace<TrialSpace>;
      using TrialLFS        = LocalFunctionSpace<TrialSpace>;
      using TrialLFSCache   = LFSIndexCache<TrialLFS>;

      static constexpr std::bool_constant<galerkin == Galerkin::automatic ?
                                          std::is_same<TrialSpace,TestSpace>::value : bool(galerkin)> isGalerkin()
      {
        return {};
      }

    private:

      LOP* _lop;

      TrialVector* _trial_vector;
      TestVector* _test_vector;

    public:

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

      ResidualEngine(TrialVector& trial_vector, TestVector& test_vector, LOP& lop,
                     std::integral_constant<Galerkin,galerkin> = std::integral_constant<Galerkin,galerkin>{})
        : _lop(&lop)
        , _trial_vector(&trial_vector)
        , _test_vector(&test_vector)
      {}

      TestVector& residual()
      {
        return *_test_vector;
      }

      TrialVector& argument()
      {
        return *_trial_vector;
      }

      const TestSpace& testSpace() const
      {
        return _test_vector->gridFunctionSpace();
      }

      const TrialSpace& trialSpace() const
      {
        return _trial_vector->gridFunctionSpace();
      }

      template<typename Assembler>
      auto context(const Assembler& assembler)
      {
        return makeContext(
          insideCell(
            assembler.entitySet(),
            cellGridData(assembler.entitySet()),
            trialSpaceData<TrialLFS,TrialLFSCache>(isGalerkin(),testSpaceData<TestLFS,TestLFSCache,CellType::Inside>(testSpace()),trialSpace()),
            cellArgumentData(CachedVectorData<UncachedVectorView<TestVector,TestLFSCache>,LocalViewDataMode::read>()),
            cellResidualData(CachedVectorData<UncachedVectorView<TestVector,TestLFSCache>,LocalViewDataMode::accumulate>()),
            extractCellContext<CellDataHolder>(*this,*_lop)
            ),
          outsideCell(
            assembler.entitySet(),
            cellGridData(assembler.entitySet()),
            trialSpaceData<TrialLFS,TrialLFSCache>(isGalerkin(),testSpaceData<TestLFS,TestLFSCache,CellType::Outside>(testSpace()),trialSpace()),
            cellArgumentData(CachedVectorData<UncachedVectorView<TestVector,TestLFSCache>,LocalViewDataMode::read>()),
            cellResidualData(CachedVectorData<UncachedVectorView<TestVector,TestLFSCache>,LocalViewDataMode::accumulate>()),
            extractCellContext<CellDataHolder>(*this,*_lop)
            ),
          cellDomainData(assembler.entitySet()),
          intersectionDomainData(assembler.entitySet()),
          ResidualEngineAssemblyFlags<isGalerkin()>(),
          extractContext<GlobalDataHolder>(*this,*_lop)
          );
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
      }

      template<typename Context>
      void result(Context& ctx)
      {}

    };



  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_RESIDUALENGINE_HH
