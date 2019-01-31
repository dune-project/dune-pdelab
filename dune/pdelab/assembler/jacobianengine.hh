// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_JACOBIANENGINE_HH
#define DUNE_PDELAB_ASSEMBLER_JACOBIANENGINE_HH

#include <dune/pdelab/common/checks.hh>
#include <dune/pdelab/common/exceptions.hh>
#include <dune/pdelab/common/typetraits.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/common/uncachedvectorview.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/assembler/context.hh>
#include <dune/pdelab/assembler/celldata.hh>
#include <dune/pdelab/assembler/functionspacedata.hh>
#include <dune/pdelab/assembler/matrixdata.hh>
#include <dune/pdelab/assembler/vectordata.hh>
#include <dune/pdelab/assembler/utility.hh>
#include <dune/pdelab/localoperator/guardedcalls.hh>
#include <dune/pdelab/gridfunctionspace/flavor.hh>
#include <dune/pdelab/assembler/enginebase.hh>

namespace Dune {
  namespace PDELab {

    template<
      typename TrialVector_,
      typename Jacobian_,
      typename LOP,
      typename TrialConstraints_ = EmptyTransformation,
      typename TestConstraints_ = EmptyTransformation,
      typename EngineParameters = DefaultJacobianEngineParameters<false,Galerkin::automatic>
      >
    class JacobianEngine
      : public InstationaryEngineBase<typename TrialVector_::value_type,EngineParameters::instationary>
    {

      static constexpr bool enable_flavors = not LocalOperator::disableFunctionSpaceFlavors<LOP>();

      using Types = LocalFunctionSpaceTypes<
        typename Jacobian_::TestSpace,
        typename TrialVector_::GridFunctionSpace,
        TrialConstraints_,
        TestConstraints_,
        enable_flavors
        >;

      using IEB = InstationaryEngineBase<typename TrialVector_::value_type,EngineParameters::instationary>;

    public:

      using IEB::instationary;

      using size_type        = std::size_t;

      using TrialVector      = TrialVector_;
      using Jacobian         = Jacobian_;
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
      std::bool_constant<EngineParameters::template galerkin<TrialSpace,TestSpace>>
      isGalerkin()
      {
        return {};
      }

    private:

      LOP* _lop;

      const TrialVector* _linearization_point;
      Jacobian* _jacobian;

      EmptyTransformation _empty_constraints;

      const TrialConstraints* _trial_constraints;
      const TestConstraints* _test_constraints;

      bool _symmetric_dirichlet_constraints;

    public:

      template<typename CellFlavor>
      struct Data
        : public Context::RootContext
      {

        using Flavor = CellFlavor;
        using Engine = JacobianEngine;
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
          return JacobianEngine::isGalerkin();
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

      JacobianEngine(
        const TrialVector& linearization_point,
        Jacobian& jacobian,
        LOP& lop,
        const TrialConstraints& trial_constraints,
        const TestConstraints& test_constraints,
        EngineParameters = {}
        )
        : _lop(&lop)
        , _linearization_point(&linearization_point)
        , _jacobian(&jacobian)
        , _trial_constraints(&trial_constraints)
        , _test_constraints(&test_constraints)
        , _symmetric_dirichlet_constraints(false)
      {}

      template<typename LOP_>
      JacobianEngine(
        const TrialVector& linearization_point,
        Jacobian& jacobian,
        LOP_& lop,
        std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,EngineParameters> = {}
        )
        : _lop(&lop)
        , _linearization_point(&linearization_point)
        , _jacobian(&jacobian)
        , _trial_constraints(&_empty_constraints)
        , _test_constraints(&_empty_constraints)
        , _symmetric_dirichlet_constraints(false)
      {}

      template<typename LOP_>
      JacobianEngine(
        Jacobian& jacobian,
        LOP_& lop,
        const TrialConstraints& trial_constraints,
        const TestConstraints& test_constraints,
        std::enable_if_t<not models<Concept::PossiblyNonLinear,LOP>() and std::is_same_v<LOP_,LOP>,EngineParameters> = {}
        )
        : _lop(&lop)
        , _linearization_point(nullptr)
        , _jacobian(&jacobian)
        , _trial_constraints(&trial_constraints)
        , _test_constraints(&test_constraints)
        , _symmetric_dirichlet_constraints(false)
      {}

      template<typename LOP_>
      JacobianEngine(
        Jacobian& jacobian,
        LOP_& lop,
        std::enable_if_t<not models<Concept::PossiblyNonLinear,LOP>() and unconstrained() and std::is_same_v<LOP_,LOP>,EngineParameters> = {}
        )
        : _lop(&lop)
        , _linearization_point(nullptr)
        , _jacobian(&jacobian)
        , _trial_constraints(&_empty_constraints)
        , _test_constraints(&_empty_constraints)
        , _symmetric_dirichlet_constraints(false)
      {}

      JacobianEngine(
        LOP& lop,
        const TrialConstraints& trial_constraints,
        const TestConstraints& test_constraints,
        EngineParameters = {}
        )
        : _lop(&lop)
        , _linearization_point(nullptr)
        , _jacobian(nullptr)
        , _trial_constraints(&trial_constraints)
        , _test_constraints(&test_constraints)
        , _symmetric_dirichlet_constraints(false)
      {}

      template<typename LOP_>
      JacobianEngine(
        LOP_& lop,
        std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,EngineParameters> = {}
        )
        : _lop(&lop)
        , _linearization_point(nullptr)
        , _jacobian(nullptr)
        , _trial_constraints(&_empty_constraints)
        , _test_constraints(&_empty_constraints)
        , _symmetric_dirichlet_constraints(false)
      {}

      void setLinearizationPoint(const TrialVector& linearization_point)
      {
        _linearization_point = &linearization_point;
      }

      void setJacobian(Jacobian& jacobian)
      {
        _jacobian = &jacobian;
      }

      Jacobian& jacobian()
      {
        return *_jacobian;
      }

      constexpr bool bindLinearizationPoint() const
      {
        return isNonLinear(localOperator());
      }

      const TrialVector& linearizationPoint() const
      {
#if DUNE_PDELAB_ENABLE_CHECK_ASSEMBLY
        if (not isNonLinear(localOperator()))
          DUNE_THROW(AssemblyError, "Not allowed to call linearizationPoint() for linear operators");
#endif
        return *_linearization_point;
      }

      const TestSpace& testSpace() const
      {
        return _jacobian->testSpace();
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
        return _linearization_point->gridFunctionSpace();
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
            skeletonTimeJacobianData(
              std::bool_constant<instationary()>{},
              skeletonJacobianData(
                cachedSkeletonMatrixData<UncachedMatrixView,Jacobian,LocalViewDataMode::accumulate>(
                  intersectionDomainData(
                    cellDomainData(
                      outsideCell(
                        extractCellContext(
                          *_lop,
                          cellTimeJacobianData(
                            std::bool_constant<instationary()>{},
                            cellJacobianData(
                              cachedMatrixData<UncachedMatrixView,Jacobian,LocalViewDataMode::accumulate>(
                                cellLinearizationPointData(
                                  models<Concept::PossiblyNonLinear,LOP>(),
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
                            cellTimeJacobianData(
                              std::bool_constant<instationary()>{},
                              cellJacobianData(
                                cachedMatrixData<UncachedMatrixView,Jacobian,LocalViewDataMode::accumulate>(
                                  cellLinearizationPointData(
                                    models<Concept::PossiblyNonLinear,LOP>(),
                                    cachedVectorData<ConstUncachedVectorView,TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                                      trialSpaceData(
                                        testSpaceData(
                                          cellGridData(
                                            timeData(
                                              Data<CellFlavor::Inside<enable_flavors>>(*this)
                                              ))))))))))))))))));
      }


      template<typename LocalContainer, typename View>
      void scatterJacobianUnconstrained(LocalContainer& local_container, View& view)
      {
        // write out all entries without considering constraints.
        for (size_type i = 0, rows = local_container.nrows() ; i < rows ; ++i)
          for (size_type j = 0, cols = local_container.ncols() ; j < cols ; ++j)
            {
              const auto& entry = local_container.getEntry(i,j);
              // skip 0 entries because they might not be present in the pattern
              if (entry == 0.0)
                continue;
              view.add(i,j,entry);
            }
      }


      template<typename LocalContainer, typename View>
      std::enable_if_t<
        Std::to_true_v<LocalContainer> and unconstrained()
        >
      scatterJacobian(LocalContainer& local_container, View& view)
      {
        scatterJacobianUnconstrained(local_container,view);
      }

      template<typename LocalContainer, typename View>
      std::enable_if_t<
        Std::to_true_v<LocalContainer> and not unconstrained()
        >
      scatterJacobian(LocalContainer& local_container, View& view)
      {

        const auto& test_indices  = view.rowIndexCache();
        const auto& trial_indices = view.colIndexCache();

        if (_symmetric_dirichlet_constraints)
          {
            for (size_type j = 0, cols = trial_indices.size() ; j < cols ; ++j)
              {
                if (trial_indices.isConstrained(j) and trial_indices.isDirichletConstraint(j))
                  {
                    // clear out the current column
                    for (size_t i = 0, rows = test_indices.size() ; i < rows ; ++i)
                      {
                        // we do not need to update the residual, since the solution
                        // (i.e. the correction) for the Dirichlet DOF is 0 by definition
                        local_container.getEntry(i,j) = 0.0;
                      }
                  }
              }
          }

        if (not (_trial_constraints->containsNonDirichletConstraints() or
                 _test_constraints->containsNonDirichletConstraints()))
          {
            // Dirichlet constraints are applied in batch after assembly is complete
            scatterJacobianUnconstrained(local_container,view);
            return;
          }

        assert(test_indices.constraintsCachingEnabled());
        assert(trial_indices.constraintsCachingEnabled());

        for (size_t i = 0, rows = test_indices.size() ; i < rows ; ++i)
          {
            bool constrained_v = test_indices.isConstrained(i);

            for (size_t j = 0, cols = trial_indices.size() ; j < cols ; ++j)
              {

                const auto& entry = local_container.getEntry(i,j);

                if (entry == 0.0)
                  continue;

                bool constrained_u = trial_indices.isConstrained(j);

                if (constrained_v)
                  {
                    if (test_indices.isDirichletConstraint(i))
                      continue;

                    for (auto vcit = test_indices.constraintsBegin(i),
                           vcend = test_indices.constraintsEnd(i) ;
                         vcit != vcend ;
                         ++vcit)
                      {
                        if (constrained_u)
                          {
                            if (trial_indices.isDirichletConstraint(j))
                              {
                                auto value = entry * vcit->weight();
                                if (value != 0.0)
                                  view.add(vcit->containerIndex(),j,value);
                              }
                            else
                              {
                                for (auto ucit = trial_indices.constraintsBegin(j),
                                       ucend = trial_indices.constraintsEnd(j) ;
                                     ucit != ucend ;
                                     ++ucit)
                                  {
                                    auto value = entry * vcit->weight() * ucit->weight();
                                    if (value != 0.0)
                                      view.add(vcit->containerIndex(),ucit->containerIndex(),value);
                                  }
                              }
                          }
                        else
                          {
                            auto value = entry * vcit->weight();
                            if (value != 0.0)
                              view.add(vcit->containerIndex(),j,value);
                          }
                      }
                  }
                else
                  {
                    if (constrained_u)
                      {
                        if (trial_indices.isDirichletConstraint(j))
                          {
                            view.add(i,j,entry);
                          }
                        else
                          {
                            for (auto ucit = trial_indices.constraintsBegin(j),
                                   ucend = trial_indices.constraintsEnd(j) ;
                                 ucit != ucend ;
                                 ++ucit)
                              {
                                auto value = entry * ucit->weight();
                                if (value != 0.0)
                                  view.add(i,ucit->containerIndex(),value);
                              }
                          }
                      }
                    else
                      view.add(i,j,entry);
                  }
              }
          }
      }

      void handleDirichletConstraints()
      {
        jacobian().flush();
        set_trivial_rows(testSpace(),jacobian(),testConstraints());
        jacobian().finalize();
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
        invoke_if_possible(LocalOperator::volumeJacobian(),*_lop,ctx.cellContext());
      }

      template<typename Context>
      void startIntersections(Context& ctx)
      {
        invoke_if_possible(LocalOperator::startIntersections(),*_lop,ctx.cellContext());
      }

      template<typename Context>
      void skeleton(Context& ctx)
      {
        invoke_if_possible(LocalOperator::skeletonJacobian(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void periodic(Context& ctx)
      {
        invoke_if_possible(LocalOperator::skeletonJacobian(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void boundary(Context& ctx)
      {
        invoke_if_possible(LocalOperator::boundaryJacobian(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void processor(Context& ctx)
      {
        invoke_if_possible(LocalOperator::boundaryJacobian(),*_lop,ctx.intersectionContext());
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
        handleDirichletConstraints();
      }

      template<typename Context>
      void result(Context& ctx)
      {}

    };


    template<typename TrialVector, typename Jacobian, typename LOP>
    JacobianEngine(
        const TrialVector&,
        Jacobian&,
        LOP&
      )
      -> JacobianEngine<
        TrialVector,
        Jacobian,
        LOP,
        EmptyTransformation,
        EmptyTransformation,
        DefaultJacobianEngineParameters<false,Galerkin::automatic>
        >;

    template<typename TrialVector, typename Jacobian, typename LOP, typename EngineParameters>
    JacobianEngine(
        const TrialVector&,
        Jacobian&,
        LOP&,
        EngineParameters
      )
      -> JacobianEngine<
        TrialVector,
        Jacobian,
        LOP,
        EmptyTransformation,
        EmptyTransformation,
        EngineParameters
        >;

    template<typename Jacobian, typename LOP>
    JacobianEngine(
        Jacobian&,
        LOP&
      )
      -> JacobianEngine<
        Backend::Vector<
          typename Jacobian::TrialGridFunctionSpace,
          typename Jacobian::value_type
          >,
        Jacobian,
        LOP,
        EmptyTransformation,
        EmptyTransformation,
        DefaultJacobianEngineParameters<false,Galerkin::automatic>
        >;


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_JACOBIANENGINE_HH
