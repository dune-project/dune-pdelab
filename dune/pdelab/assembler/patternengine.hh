// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_PATTERNENGINE_HH
#define DUNE_PDELAB_ASSEMBLER_PATTERNENGINE_HH

#include <memory>

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
      typename TrialVector_,
      typename TestVector_,
      typename LOP,
      typename MBE,
      typename TrialConstraints_ = EmptyTransformation,
      typename TestConstraints_ = EmptyTransformation,
      typename JF = typename TrialVector_::value_type,
      Galerkin galerkin = Galerkin::automatic
      >
    class PatternEngine
    {

      static constexpr bool enable_flavors = not LocalOperator::disableFunctionSpaceFlavors<LOP>();

      using Types = LocalFunctionSpaceTypes<
        typename TestVector_::GridFunctionSpace,
        typename TrialVector_::GridFunctionSpace,
        enable_flavors
        >;

    public:

      using size_type        = std::size_t;

      using TrialVector      = TrialVector_;
      using TrialConstraints = TrialConstraints_;

      using TestVector       = TestVector_;
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
        const TrialVector& trial_vector,
        const TestVector& test_vector,
        LOP& lop,
        const MatrixBackend& matrix_backend,
        const TrialConstraints& trial_constraints,
        const TestConstraints& test_constraints,
        std::integral_constant<Galerkin,galerkin> = std::integral_constant<Galerkin,galerkin>{}
        )
        : _lop(&lop)
        , _matrix_backend(matrix_backend)
        , _trial_space(&trial_vector.gridFunctionSpace())
        , _test_space(&test_vector.gridFunctionSpace())
        , _trial_constraints(&trial_constraints)
        , _test_constraints(&test_constraints)
      {}

      template<typename LOP_>
      PatternEngine(
        const TrialVector& trial_vector,
        const TestVector& test_vector,
        LOP_& lop,
        const MatrixBackend& matrix_backend,
        std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,std::integral_constant<Galerkin,galerkin>> = std::integral_constant<Galerkin,galerkin>{}
        )
        : _lop(&lop)
        , _matrix_backend(matrix_backend)
        , _trial_space(&trial_vector.gridFunctionSpace())
        , _test_space(&test_vector.gridFunctionSpace())
        , _trial_constraints(&_empty_constraints)
        , _test_constraints(&_empty_constraints)
      {}

      const TestSpace& testSpace() const
      {
        return *_test_space;
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
        return *_trial_space;
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
                        for (auto it = trial_indices.constraintsBegin(),
                               end = trial_indices.constraintsEnd() ;
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
                    for (auto vit = test_indices.constraintsBegin(),
                           end = test_indices.constraintsEnd() ;
                         vit != end ;
                         ++vit)
                      {
                        if (not constrained_u or trial_indices.isDirichletConstraint(j))
                          {
                            _pattern->add_link(vit->containerIndex(),trial_indices.containerIndex(j));
                          }
                        else
                          {
                            for (auto uit = trial_indices.constraintsBegin(),
                                   end = trial_indices.constraintsEnd() ;
                                 uit != end ;
                                 ++uit)
                              _pattern->add_link(vit->containerIndex(),uit->containerIndex());
                          }
                      }
                  }
              }
          }
      }



      template<typename Context>
      void start(Context& ctx)
      {
        _pattern = std::make_shared<Pattern>(_matrix_backend.template makePattern<Matrix>(testSpace(),trialSpace()));
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

    template<typename TrialVector, typename TestVector, typename LOP, typename MBE>
    PatternEngine(
        const TrialVector&,
        const TestVector&,
        LOP&,
        const MBE&
      )
      -> PatternEngine<
        TrialVector,
        TestVector,
        LOP,
        MBE,
        EmptyTransformation,
        EmptyTransformation,
        typename TrialVector::value_type,
        Galerkin::automatic
        >;

    template<typename TrialVector, typename TestVector, typename LOP, typename MBE, Galerkin galerkin>
    PatternEngine(
        const TrialVector&,
        const TestVector&,
        LOP&,
        const MBE&,
        std::integral_constant<Galerkin,galerkin>
      )
      -> PatternEngine<
        TrialVector,
        TestVector,
        LOP,
        MBE,
        EmptyTransformation,
        EmptyTransformation,
        typename TrialVector::value_type,
        galerkin
        >;

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_PATTERNENGINE_HH
