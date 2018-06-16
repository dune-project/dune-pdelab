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
      using TestSpaceCache  = LFSIndexCache<TestLFS>;

      using TrialSpace      = typename TrialVector::GridFunctionSpace;
      using TrialLocalSpace = LocalFunctionSpace<TrialSpace>;
      using TrialLFS        = LocalFunctionSpace<TrialSpace>;
      using TrialSpaceCache = LFSIndexCache<TrialLFS>;

      using EntitySet       = typename TestSpace::Traits::EntitySet;

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

      TestSpaceCache makeTestSpaceCache() const
      {
        return TestSpaceCache();
      }

      const TrialSpace& trialSpace() const
      {
        return _trial_vector->gridFunctionSpace();
      }

      TrialSpaceCache makeTrialSpaceCache() const
      {
        return TrialSpaceCache();
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
                    cellResidualData(
                      cachedVectorData<UncachedVectorView<TestVector,TestSpaceCache>,LocalViewDataMode::accumulate>(
                        cellArgumentData(
                          cachedVectorData<UncachedVectorView<TrialVector,TrialSpaceCache>,LocalViewDataMode::read>(
                            trialSpaceData(
                              testSpaceData(
                                cellGridData(
                                  Data<CellFlavor::Outside>(*this)
                                  )))))))),
                  insideCell(
                    extractCellContext(
                      *_lop,
                      cellResidualData(
                        cachedVectorData<UncachedVectorView<TestVector,TestSpaceCache>,LocalViewDataMode::accumulate>(
                          cellArgumentData(
                            cachedVectorData<UncachedVectorView<TrialVector,TrialSpaceCache>,LocalViewDataMode::read>(
                              trialSpaceData(
                                testSpaceData(
                                  cellGridData(
                                    Data<CellFlavor::Inside>(*this)
                                    )))))))))))));
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
