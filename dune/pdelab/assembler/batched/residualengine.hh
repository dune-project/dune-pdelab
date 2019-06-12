// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_BATCHED_RESIDUALENGINE_HH
#define DUNE_PDELAB_ASSEMBLER_BATCHED_RESIDUALENGINE_HH

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
      typename EngineParameters = DefaultResidualEngineParameters<false,Galerkin::automatic>
      >
    class BatchedResidualEngine
      : public ResidualEngine<TrialVector_,TestVector_,LOP,TrialConstraints_,TestConstraints_,EngineParameters>
    {


      template<typename LOP_>
      ResidualEngine(
        const TrialVector& trial_vector,
        TestVector& test_vector,
        LOP_& lop,
        std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,EngineParameters> = {}
        )
        : _lop(&lop)
        , _argument(&trial_vector)
        , _residual(&test_vector)
      {}

      ResidualEngine(
        const TrialVector& trial_vector,
        TestVector& test_vector,
        LOP& lop,
        const TrialConstraints& trial_constraints,
        const TestConstraints& test_constraints,
        EngineParameters = {}
        )
        : FSP(&trial_constraints,&test_constraints)
        , _lop(&lop)
        , _argument(&trial_vector)
        , _residual(&test_vector)
      {}

      template<typename LOP_>
      ResidualEngine(
        LOP_& lop,
        std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,EngineParameters> = {}
        )
        : _lop(&lop)
        , _argument(nullptr)
        , _residual(nullptr)
      {}

      ResidualEngine(
        LOP& lop,
        const TrialConstraints& trial_constraints,
        const TestConstraints& test_constraints,
        EngineParameters = {}
        )
        : FSP(&trial_constraints,&test_constraints)
        , _lop(&lop)
        , _argument(nullptr)
        , _residual(nullptr)
      {}

      template<typename Context>
      void volume(Context& ctx)
      {
        invoke_if_possible(LocalOperator::volumeIntegral(),*_lop,ctx);
      }

      template<typename Context>
      void startIntersections(Context& ctx)
      {
        invoke_if_possible(LocalOperator::startIntersections(),*_lop,ctx);
      }

      template<typename Context>
      void skeleton(Context& ctx)
      {
        invoke_if_possible(LocalOperator::skeletonIntegral(),*_lop,ctx);
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
        DefaultResidualEngineParameters<false,Galerkin::automatic>
        >;


    template<typename Argument, typename Residual, typename LOP, typename EngineParameters>
    ResidualEngine(
        const Argument&,
        Residual&,
        LOP&,
        EngineParameters
      )
      -> ResidualEngine<
        Argument,
        Residual,
        LOP,
        EmptyTransformation,
        EmptyTransformation,
        EngineParameters
        >;


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_BATCHED_RESIDUALENGINE_HH
