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
#include <dune/pdelab/gridfunctionspace/flavor.hh>

namespace Dune {
  namespace PDELab {

    enum class Galerkin { disable, enable, automatic };

    static constexpr auto enableGalerkin = std::integral_constant<Galerkin,Galerkin::enable>{};
    static constexpr auto disableGalerkin = std::integral_constant<Galerkin,Galerkin::disable>{};
    static constexpr auto automaticGalerkin = std::integral_constant<Galerkin,Galerkin::automatic>{};

    template<typename TestSpace_, typename TrialSpace_, bool enable_flavors>
    struct LocalFunctionSpaceTypes;

    template<typename TestSpace_, typename TrialSpace_>
    struct LocalFunctionSpaceTypes<TestSpace_,TrialSpace_,true>
    {

      using TestSpace       = TestSpace_;

      template<typename Flavor = Flavor::Generic>
      using TestLocalSpace  = LocalFunctionSpace<TestSpace,Flavor>;

      template<typename Flavor = Flavor::Generic>
      using TestSpaceCache  = LFSIndexCache<TestLocalSpace<Flavor>>;


      using TrialSpace      = TrialSpace_;

      template<typename Flavor = Flavor::Generic>
      using TrialLocalSpace = LocalFunctionSpace<TrialSpace,Flavor>;

      template<typename Flavor = Flavor::Generic>
      using TrialSpaceCache = LFSIndexCache<TrialLocalSpace<Flavor>>;

    };

    template<typename TestSpace_, typename TrialSpace_>
    struct LocalFunctionSpaceTypes<TestSpace_,TrialSpace_,false>
    {

      using TestSpace       = TestSpace_;

      template<typename = Flavor::Generic>
      using TestLocalSpace  = LocalFunctionSpace<TestSpace,Flavor::Generic>;

      template<typename = Flavor::Generic>
      using TestSpaceCache  = LFSIndexCache<TestLocalSpace<Flavor::Generic>>;


      using TrialSpace      = TrialSpace_;

      template<typename = Flavor::Generic>
      using TrialLocalSpace = LocalFunctionSpace<TrialSpace,Flavor::Generic>;

      template<typename = Flavor::Generic>
      using TrialSpaceCache = LFSIndexCache<TrialLocalSpace<Flavor::Generic>>;

    };


    template<typename TrialVector, typename TestVector, typename LOP, Galerkin galerkin = Galerkin::automatic>
    class ResidualEngine
    {

      static constexpr bool enable_flavors = not LocalOperator::disableFunctionSpaceFlavors<LOP>();

      using Types = LocalFunctionSpaceTypes<
        typename TestVector::GridFunctionSpace,
        typename TrialVector::GridFunctionSpace,
        enable_flavors
        >;

    public:

      using TestSpace       = typename Types::TestSpace;

      template<typename Flavor = Flavor::Generic>
      using TestLocalSpace  = typename Types::template TestLocalSpace<Flavor>;

      template<typename Flavor = Flavor::Generic>
      using TestSpaceCache  = typename Types::template TestSpaceCache<Flavor>;


      using TrialSpace      = typename Types::TrialSpace;

      template<typename Flavor = Flavor::Generic>
      using TrialLocalSpace = typename Types::template TrialLocalSpace<Flavor>;

      template<typename Flavor = Flavor::Generic>
      using TrialSpaceCache = typename Types::template TrialSpaceCache<Flavor>;


      using EntitySet       = typename TestSpace::Traits::EntitySet;

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

      template<typename Flavor_>
      std::enable_if_t<
        Std::to_true_type<Flavor_>::value and enable_flavors,
        TestSpaceCache<Flavor_>
        >
      makeTestSpaceCache(Flavor_) const
      {
        return TestSpaceCache<Flavor_>();
      }

      template<typename Flavor_>
      std::enable_if_t<
        Std::to_true_type<Flavor_>::value and not enable_flavors,
        TestSpaceCache<Flavor::Generic>
        >
      makeTestSpaceCache(Flavor_) const
      {
        return TestSpaceCache<Flavor::Generic>();
      }

      const TrialSpace& trialSpace() const
      {
        return _trial_vector->gridFunctionSpace();
      }

      template<typename Flavor_>
      std::enable_if_t<
        Std::to_true_type<Flavor_>::value and enable_flavors,
        TrialSpaceCache<Flavor_>
        >
      makeTrialSpaceCache(Flavor_) const
      {
        return TrialSpaceCache<Flavor_>();
      }

      template<typename Flavor_>
      std::enable_if_t<
        Std::to_true_type<Flavor_>::value and not enable_flavors,
        TrialSpaceCache<Flavor::Generic>
        >
      makeTrialSpaceCache(Flavor_) const
      {
        return TrialSpaceCache<Flavor::Generic>();
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
                    cellResidualData(
                      cachedVectorData<UncachedVectorView,TestVector,Flavor::Test,LocalViewDataMode::accumulate>(
                        cellArgumentData(
                          cachedVectorData<UncachedVectorView,TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                            trialSpaceData(
                              testSpaceData(
                                cellGridData(
                                  Data<CellFlavor::Outside<enable_flavors>>(*this)
                                  )))))))),
                  insideCell(
                    extractCellContext(
                      *_lop,
                      cellResidualData(
                        cachedVectorData<UncachedVectorView,TestVector,Flavor::Test,LocalViewDataMode::accumulate>(
                          cellArgumentData(
                            cachedVectorData<UncachedVectorView,TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                              trialSpaceData(
                                testSpaceData(
                                  cellGridData(
                                    Data<CellFlavor::Inside<enable_flavors>>(*this)
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
