// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_PATTERNENGINE_HH
#define DUNE_PDELAB_ASSEMBLER_PATTERNENGINE_HH

#include <memory>
#include <vector>

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
      typename EngineParameters = DefaultPatternEngineParameters<Galerkin::automatic>
      >
    class PatternEngine
      : public FunctionSpaceProvider<TrialSpace_,
                                     TestSpace_,
                                     TrialConstraints_,
                                     TestConstraints_,
                                     not LocalOperator::disableFunctionSpaceFlavors<LOP>()
                                     >
    {

      static constexpr bool enable_flavors = not LocalOperator::disableFunctionSpaceFlavors<LOP>();

      using FSP = FunctionSpaceProvider<
        TrialSpace_,
        TestSpace_,
        TrialConstraints_,
        TestConstraints_,
        enable_flavors
        >;

      using Types = typename FSP::Types;

    public:

      using FSP::unconstrained;
      using FSP::trialConstraints;
      using FSP::testConstraints;

      using size_type        = std::size_t;

      using TrialVector      = Dune::PDELab::Backend::Vector<TrialSpace_,JF>;
      using TrialConstraints = typename Types::TrialConstraints;

      using TestVector       = Dune::PDELab::Backend::Vector<TestSpace_,JF>;
      using TestConstraints  = typename Types::TestConstraints;

      using TestSpace        = typename Types::TestSpace;

      using TrialSpace       = typename Types::TrialSpace;

      using MatrixBackend    = MBE;
      using Matrix           = Dune::PDELab::Backend::Matrix<MBE,TrialVector,TestVector,JF>;
      using Pattern          = typename Matrix::Pattern;

      using EntitySet        = typename TestSpace::Traits::EntitySet;

      static constexpr
      std::bool_constant<EngineParameters::template galerkin<TrialSpace,TestSpace>>
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

    public:

      template<typename CellFlavor>
      struct Data
        : public Context::RootContext
      {

        using Flavor = CellFlavor;
        using Engine = PatternEngine;
        using EntitySet = typename Engine::EntitySet;

        static constexpr auto isGalerkin()
        {
          return PatternEngine::isGalerkin();
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
        EngineParameters = {}
        )
        : FSP(&trial_constraints,&test_constraints)
        , _lop(&lop)
        , _matrix_backend(matrix_backend)
        , _trial_space(&trial_space)
        , _test_space(&test_space)
      {}

      template<typename LOP_>
      PatternEngine(
        const TrialSpace& trial_space,
        const TestSpace& test_space,
        LOP_& lop,
        const MatrixBackend& matrix_backend,
        JF,
        std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,EngineParameters> = {}
        )
        : _lop(&lop)
        , _matrix_backend(matrix_backend)
        , _trial_space(&trial_space)
        , _test_space(&test_space)
      {}

      const TestSpace& testSpace() const
      {
        return *_test_space;
      }

      const TrialSpace& trialSpace() const
      {
        return *_trial_space;
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


      template<typename Context, typename LocalPattern, typename TestIndices, typename TrialIndices>
      std::enable_if_t<Std::to_true_type_v<LocalPattern> and unconstrained()>
      scatterPattern(
        const Context& ctx,
        LocalPattern& local_pattern,
        const TestIndices& test_indices,
        const TrialIndices& trial_indices
        )
      {
        // in fastDG mode, there can only be one entry with index (0,0)
        if (ctx.fastDG() and local_pattern.size() > 0)
        {
          assert (local_pattern.size() == 1);
          auto [i,j] = *local_pattern.begin();
          assert (i == 0);
          assert (j == 0);
        }

        for (auto [i,j] : local_pattern)
          _pattern->add_link(test_indices.containerIndex(i),trial_indices.containerIndex(j));
      }

      template<typename Context, typename LocalPattern, typename TestIndices, typename TrialIndices>
      std::enable_if_t<Std::to_true_type_v<LocalPattern> and not unconstrained()>
      scatterPattern(
        const Context& ctx,
        LocalPattern& local_pattern,
        const TestIndices& test_indices,
        const TrialIndices& trial_indices
        )
      {
        // in fastDG mode, there can only be one entry with index (0,0)
        if (ctx.fastDG() and local_pattern.size() > 0)
        {
          assert (local_pattern.size() == 1);
          auto [i,j] = *local_pattern.begin();
          assert (i == 0);
          assert (j == 0);
        }

        for (auto [i,j] : local_pattern)
        {
            auto constrained_v = test_indices.isConstrained(i);
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
        DefaultPatternEngineParameters<Galerkin::automatic>
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
        DefaultPatternEngineParameters<Galerkin::automatic>
        >;

    template<typename TrialSpace, typename TestSpace, typename LOP, typename MBE, typename JF, typename EngineParameters>
    PatternEngine(
        const TrialSpace&,
        const TestSpace&,
        LOP&,
        const MBE&,
        JF,
        EngineParameters
      )
      -> PatternEngine<
        TrialSpace,
        TestSpace,
        LOP,
        MBE,
        JF,
        EmptyTransformation,
        EmptyTransformation,
        EngineParameters
        >;

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_PATTERNENGINE_HH
