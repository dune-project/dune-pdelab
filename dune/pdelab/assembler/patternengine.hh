// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_PATTERNENGINE_HH
#define DUNE_PDELAB_ASSEMBLER_PATTERNENGINE_HH

#include <memory>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/assembler/context.hh>
#include <dune/pdelab/assembler/celldata.hh>
#include <dune/pdelab/assembler/functionspacedata.hh>
#include <dune/pdelab/assembler/vectordata.hh>
#include <dune/pdelab/assembler/patterndata.hh>
#include <dune/pdelab/localoperator/guardedcalls.hh>
#include <dune/pdelab/gridfunctionspace/flavor.hh>
#include <dune/pdelab/assembler/enginebase.hh>

namespace Dune {
  namespace PDELab {

    template<
      typename TrialSpace_,
      typename TestSpace_,
      typename LOP,
      typename MBE,
      typename JF,
      typename TrialConstraints_ = EmptyTransformation,
      typename TestConstraints_ = EmptyTransformation,
      Galerkin galerkin = Galerkin::automatic
      >
    class PatternEngine
    {

      static constexpr bool enable_flavors = not LocalOperator::disableFunctionSpaceFlavors<LOP>();

      using Types = LocalFunctionSpaceTypes<
        TrialSpace_,
        TestSpace_,
        TrialConstraints_,
        TestConstraints_,
        enable_flavors
        >;

    public:

      using size_type        = std::size_t;

      using TrialVector      = Dune::PDELab::Backend::Vector<TrialSpace_,JF>;
      using TrialConstraints = TrialConstraints_;

      using TestVector       = Dune::PDELab::Backend::Vector<TestSpace_,JF>;
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

      using MatrixBackend    = MBE;
      using Matrix           = Dune::PDELab::Backend::Matrix<MBE,TrialVector,TestVector,JF>;
      using Pattern          = typename Matrix::Pattern;

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
      MatrixBackend _matrix_backend;

      const TrialSpace* _trial_space;
      const TestSpace* _test_space;
      std::shared_ptr<Pattern> _pattern;

      EmptyTransformation _empty_constraints;

      const TrialConstraints* _trial_constraints;
      const TestConstraints* _test_constraints;


    public:

      template<typename CellFlavor>
      struct Data
        : public Context::RootContext
      {

        using Flavor = CellFlavor;
        using Engine = PatternEngine;
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
          return false;
        }

        static constexpr auto isGalerkin()
        {
          return PatternEngine::isGalerkin();
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
        return std::is_invocable_v<decltype(LocalOperator::boundaryPattern()),LOP,PlaceHolder&>;
      }

      static constexpr bool visitSkeletonIntersections()
      {
        return std::is_invocable_v<decltype(LocalOperator::skeletonPattern()),LOP,PlaceHolder&>;
      }

      static constexpr bool visitPeriodicIntersections()
      {
        return visitSkeletonIntersections();
      }

      static constexpr bool visitProcessorIntersections()
      {
        return visitSkeletonIntersections();
      }

      PatternEngine(
        const TrialSpace& trial_space,
        const TestSpace& test_space,
        LOP& lop,
        const MatrixBackend& matrix_backend,
        const TrialConstraints& trial_constraints,
        const TestConstraints& test_constraints,
        JF,
        std::integral_constant<Galerkin,galerkin> = std::integral_constant<Galerkin,galerkin>{}
        )
        : _lop(&lop)
        , _matrix_backend(matrix_backend)
        , _trial_space(&trial_space)
        , _test_space(&test_space)
        , _trial_constraints(&trial_constraints)
        , _test_constraints(&test_constraints)
      {}

      template<typename LOP_>
      PatternEngine(
        const TrialSpace& trial_space,
        const TestSpace& test_space,
        LOP_& lop,
        const MatrixBackend& matrix_backend,
        JF,
        std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,std::integral_constant<Galerkin,galerkin>> = std::integral_constant<Galerkin,galerkin>{}
        )
        : _lop(&lop)
        , _matrix_backend(matrix_backend)
        , _trial_space(&trial_space)
        , _test_space(&test_space)
        , _trial_constraints(&_empty_constraints)
        , _test_constraints(&_empty_constraints)
      {}

      const TestSpace& testSpace() const
      {
        return *_test_space;
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
        return *_trial_space;
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

      template<typename Assembler>
      auto context(const Assembler& assembler)
      {
        return
          skeletonPatternData(
            intersectionDomainData(
              cellDomainData(
                outsideCell(
                  cellPatternData(
                    trialSpaceData(
                      testSpaceData(
                        cellGridData(
                          Data<CellFlavor::Outside<enable_flavors>>(*this)
                          )))),
                  insideCell(
                    cellPatternData(
                      trialSpaceData(
                        testSpaceData(
                          cellGridData(
                            Data<CellFlavor::Inside<enable_flavors>>(*this)
                            )))))))));
      }


      template<typename LocalPattern, typename TestIndices, typename TrialIndices>
      std::enable_if_t<Std::to_true_v<LocalPattern> and unconstrained()>
      scatterPattern(LocalPattern& local_pattern, const TestIndices& test_indices, const TrialIndices& trial_indices)
      {
        for (size_type i = 0, rows = test_indices.size() ; i < rows ; ++i)
          for (size_type j = 0, cols = trial_indices.size() ; j < cols ; ++j)
            _pattern->add_link(test_indices.containerIndex(i),trial_indices.containerIndex(j));
      }


      template<typename LocalPattern, typename TestIndices, typename TrialIndices>
      std::enable_if_t<Std::to_true_v<LocalPattern> and not unconstrained()>
      scatterPattern(LocalPattern& local_pattern, const TestIndices& test_indices, const TrialIndices& trial_indices)
      {
        for (size_type i = 0, rows = test_indices.size() ; i < rows ; ++i)
          {
            auto constrained_v = test_indices.isConstrained(i);
            for (size_type j = 0, cols = trial_indices.size() ; j < cols ; ++j)
              {
                auto constrained_u = trial_indices.isConstrained(j);
                add_diagonal_entry(*_pattern,test_indices.containerIndex(i),trial_indices.containerIndex(j));

                if (not constrained_v)
                  {
                    if (constrained_u and not trial_indices.isDirichletConstraint(j))
                      {
                        for (auto it = trial_indices.constraintsBegin(j),
                               end = trial_indices.constraintsEnd(j) ;
                             it != end ;
                             ++it)
                          _pattern->add_link(test_indices.containerIndex(i),it->containerIndex());
                      }
                    else
                      {
                        _pattern->add_link(test_indices.containerIndex(i),trial_indices.containerIndex(j));
                      }
                  }
                else if (not test_indices.isDirichletConstraint(i))
                  {
                    for (auto vit = test_indices.constraintsBegin(i),
                           end = test_indices.constraintsEnd(i) ;
                         vit != end ;
                         ++vit)
                      {
                        if (not constrained_u or trial_indices.isDirichletConstraint(j))
                          {
                            _pattern->add_link(vit->containerIndex(),trial_indices.containerIndex(j));
                          }
                        else
                          {
                            for (auto uit = trial_indices.constraintsBegin(j),
                                   end = trial_indices.constraintsEnd(j) ;
                                 uit != end ;
                                 ++uit)
                              _pattern->add_link(vit->containerIndex(),uit->containerIndex());
                          }
                      }
                  }
                else
                {
                  _pattern->add_link(test_indices.containerIndex(i),trial_indices.containerIndex(j));
                }
              }
          }
      }



      template<typename Context>
      void start(Context& ctx)
      {
        _pattern = _matrix_backend.template makePattern<Matrix>(testSpace(),trialSpace());
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
        invoke_if_possible(LocalOperator::volumePattern(),*_lop,ctx.cellContext());
      }

      template<typename Context>
      void startIntersections(Context& ctx)
      {
        invoke_if_possible(LocalOperator::startIntersections(),*_lop,ctx.cellContext());
      }

      template<typename Context>
      void skeleton(Context& ctx)
      {
        invoke_if_possible(LocalOperator::skeletonPattern(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void periodic(Context& ctx)
      {
        invoke_if_possible(LocalOperator::skeletonPattern(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void boundary(Context& ctx)
      {
        invoke_if_possible(LocalOperator::boundaryPattern(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void processor(Context& ctx)
      {
        invoke_if_possible(LocalOperator::boundaryPattern(),*_lop,ctx.intersectionContext());
      }

      template<typename Context>
      void finishIntersections(Context& ctx)
      {
        invoke_if_possible(LocalOperator::finishIntersections(),*_lop,ctx.cellContext());
      }
      template<typename Context>
      void volumePostIntersections(Context& ctx)
      {
        invoke_if_possible(LocalOperator::volumePatternPostIntersections(),*_lop,ctx.cellContext());
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
      std::shared_ptr<Pattern> result(Context& ctx)
      {
        return std::move(_pattern);
      }

    };

    template<typename TrialSpace, typename TestSpace, typename LOP, typename MBE, typename JF>
    PatternEngine(
        const TrialSpace&,
        const TestSpace&,
        LOP&,
        const MBE&,
        JF
      )
      -> PatternEngine<
        TrialSpace,
        TestSpace,
        LOP,
        MBE,
        JF,
        EmptyTransformation,
        EmptyTransformation,
        Galerkin::automatic
        >;

    template<
      typename TrialSpace,
      typename TestSpace,
      typename LOP,
      typename MBE,
      typename JF,
      typename TrialConstraints,
      typename TestConstraints
      >
    PatternEngine(
        const TrialSpace&,
        const TestSpace&,
        LOP&,
        const MBE&,
        const TrialConstraints&,
        const TestConstraints&,
        JF
      )
      -> PatternEngine<
        TrialSpace,
        TestSpace,
        LOP,
        MBE,
        JF,
        TrialConstraints,
        TestConstraints,
        Galerkin::automatic
        >;

    template<typename TrialSpace, typename TestSpace, typename LOP, typename MBE, typename JF, Galerkin galerkin>
    PatternEngine(
        const TrialSpace&,
        const TestSpace&,
        LOP&,
        const MBE&,
        JF,
        std::integral_constant<Galerkin,galerkin>
      )
      -> PatternEngine<
        TrialSpace,
        TestSpace,
        LOP,
        MBE,
        JF,
        EmptyTransformation,
        EmptyTransformation,
        galerkin
        >;

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_PATTERNENGINE_HH
