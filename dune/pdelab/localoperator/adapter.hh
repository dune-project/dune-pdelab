// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_ADAPTER_HH
#define DUNE_PDELAB_LOCALOPERATOR_ADAPTER_HH

#include <type_traits>
#include <dune/pdelab/common/typetraits.hh>

#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/localoperator/callswitch.hh>

namespace Dune {
  namespace PDELab {

    template<typename LOP>
    class LocalOperatorAdapter
    {

      template<bool enabled>
      using Call = LocalAssemblerCallSwitch<LOP,enabled>;

    public:

      static constexpr bool disableFunctionSpaceFlavors()
      {
        return true;
      }

      static constexpr bool intersectionsTwoSided()
      {
        return LOP::doSkeletonTwoSided;
      }

      static constexpr bool isNonLinear()
      {
        return true;
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and (LOP::doAlphaVolume or LOP::doLambdaVolume)>
      volumeIntegral(Context& ctx) const
      {
        auto residual = ctx.residual();
        auto eg = ElementGeometry<typename Context::Entity>(ctx.entity());
        Call<LOP::doAlphaVolume and Context::assembleVariablePart()>::alpha_volume(
          _lop,eg,
          ctx.trial().functionSpace(),ctx.coefficient(),
          ctx.test().functionSpace(),residual
          );
        Call<LOP::doLambdaVolume and Context::assembleConstantPart()>::lambda_volume(
          _lop,eg,
          ctx.test().functionSpace(),residual
          );
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and (LOP::doAlphaVolumePostSkeleton or LOP::doLambdaVolumePostSkeleton)>
      volumeIntegralPostIntersections(Context& ctx) const
      {
        auto residual = ctx.residual();
        auto eg = ElementGeometry<typename Context::Entity>(ctx.entity());
        Call<LOP::doAlphaVolumePostSkeleton and Context::assembleVariablePart()>::alpha_volume_post_skeleton(
          _lop,eg,
          ctx.trial().functionSpace(),ctx.coefficient(),
          ctx.test().functionSpace(),residual
          );
        Call<LOP::doLambdaVolumePostSkeleton and Context::assembleConstantPart()>::lambda_volume_post_skeleton(
          _lop,eg,
          ctx.test().functionSpace(),residual
          );
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and (LOP::doAlphaBoundary or LOP::doLambdaBoundary)>
      boundaryIntegral(Context& ctx) const
      {
        auto residual = ctx.inside().residual();
        auto ig = IntersectionGeometry<typename Context::Domain::Intersection>(ctx.domain().intersection(),ctx.domain().index());
        Call<LOP::doAlphaBoundary and Context::assembleVariablePart()>::alpha_boundary(
          _lop,ig,
          ctx.inside().trial().functionSpace(),ctx.inside().coefficient(),
          ctx.inside().test().functionSpace(),residual
          );
        Call<LOP::doLambdaBoundary and Context::assembleConstantPart()>::lambda_boundary(
          _lop,ig,
          ctx.inside().test().functionSpace(),residual
          );
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and (LOP::doAlphaSkeleton or LOP::doLambdaSkeleton)>
      skeletonIntegral(Context& ctx) const
      {
        auto inside_residual = ctx.inside().residual();
        auto outside_residual = ctx.outside().residual();
        auto ig = IntersectionGeometry<typename Context::Domain::Intersection>(ctx.domain().intersection(),ctx.domain().index());
        Call<LOP::doAlphaSkeleton and Context::assembleVariablePart()>::alpha_skeleton(
          _lop,ig,
          ctx.inside().trial().functionSpace(),ctx.inside().coefficient(),ctx.inside().test().functionSpace(),
          ctx.outside().trial().functionSpace(),ctx.outside().coefficient(),ctx.outside().test().functionSpace(),
          inside_residual,outside_residual
          );
        Call<LOP::doLambdaSkeleton and Context::assembleConstantPart()>::lambda_skeleton(
          _lop,ig,
          ctx.inside().test().functionSpace(),ctx.outside().test().functionSpace(),
          inside_residual,outside_residual
          );
      }


      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and LOP::doAlphaVolume>
      volumeApplyJacobian(Context& ctx) const
      {
        auto result = ctx.result();
        auto eg = ElementGeometry<typename Context::Entity>(ctx.entity());
        Call<LOP::doAlphaVolume and Context::assembleVariablePart()>::jacobian_apply_volume(
          _lop,eg,
          ctx.trial().functionSpace(),ctx.coefficient(),
          ctx.test().functionSpace(),result
          );
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and LOP::doAlphaVolumePostSkeleton>
      volumeApplyJacobianPostIntersections(Context& ctx) const
      {
        auto result = ctx.result();
        auto eg = ElementGeometry<typename Context::Entity>(ctx.entity());
        Call<LOP::doAlphaVolumePostSkeleton and Context::assembleVariablePart()>::jacobian_apply_volume_post_skeleton(
          _lop,eg,
          ctx.trial().functionSpace(),ctx.coefficient(),
          ctx.test().functionSpace(),result
          );
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and LOP::doAlphaBoundary>
      boundaryApplyJacobian(Context& ctx) const
      {
        auto result = ctx.inside().result();
        auto ig = IntersectionGeometry<typename Context::Domain::Intersection>(ctx.domain().intersection(),ctx.domain().index());
        Call<LOP::doAlphaBoundary and Context::assembleVariablePart()>::jacobian_apply_boundary(
          _lop,ig,
          ctx.inside().trial().functionSpace(),ctx.inside().coefficient(),
          ctx.inside().test().functionSpace(),result
          );
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and LOP::doAlphaSkeleton>
      skeletonApplyJacobian(Context& ctx) const
      {
        auto inside_result = ctx.inside().result();
        auto outside_result = ctx.outside().result();
        auto ig = IntersectionGeometry<typename Context::Domain::Intersection>(ctx.domain().intersection(),ctx.domain().index());
        Call<LOP::doAlphaSkeleton and Context::assembleVariablePart()>::jacobian_apply_skeleton(
          _lop,ig,
          ctx.inside().trial().functionSpace(),ctx.inside().coefficient(),ctx.inside().test().functionSpace(),
          ctx.outside().trial().functionSpace(),ctx.outside().coefficient(),ctx.outside().test().functionSpace(),
          inside_result,outside_result
          );
      }


      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and LOP::doAlphaVolume>
      volumeJacobian(Context& ctx) const
      {
        auto jacobian = ctx.inside().jacobian();
        auto eg = ElementGeometry<typename Context::Entity>(ctx.entity());
        Call<LOP::doAlphaVolume>::jacobian_volume(
          _lop,eg,
          ctx.trial().functionSpace(),ctx.linearizationPoint(),
          ctx.test().functionSpace(),jacobian
          );
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and LOP::doAlphaVolumePostSkeleton>
      volumeJacobianPostIntersections(Context& ctx) const
      {
        auto jacobian = ctx.inside().jacobian();
        auto eg = ElementGeometry<typename Context::Entity>(ctx.entity());
        Call<LOP::doAlphaVolume>::jacobian_volume_post_skeleton(
          _lop,eg,
          ctx.trial().functionSpace(),ctx.linearizationPoint(),
          ctx.test().functionSpace(),jacobian
          );
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and LOP::doAlphaBoundary>
      boundaryJacobian(Context& ctx) const
      {
        auto jacobian = ctx.inside().jacobian();
        auto ig = IntersectionGeometry<typename Context::Domain::Intersection>(ctx.domain().intersection(),ctx.domain().index());
        Call<LOP::doAlphaBoundary>::jacobian_boundary(
          _lop,ig,
          ctx.inside().trial().functionSpace(),ctx.inside().linearizationPoint(),
          ctx.inside().test().functionSpace(),jacobian
          );
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and LOP::doAlphaSkeleton>
      skeletonJacobian(Context& ctx) const
      {
        auto jac_ss = ctx.inside().jacobian();
        auto jac_nn = ctx.outside().jacobian();
        auto jac_sn = ctx.jacobian(ctx.inside(),ctx.outside());
        auto jac_ns = ctx.jacobian(ctx.outside(),ctx.inside());
        auto ig = IntersectionGeometry<typename Context::Domain::Intersection>(ctx.domain().intersection(),ctx.domain().index());
        Call<LOP::doAlphaSkeleton>::jacobian_skeleton(
          _lop,ig,
          ctx.inside().trial().functionSpace(),ctx.inside().linearizationPoint(),ctx.inside().test().functionSpace(),
          ctx.outside().trial().functionSpace(),ctx.outside().linearizationPoint(),ctx.outside().test().functionSpace(),
          jac_ss,jac_sn,jac_ns,jac_nn
          );
      }



      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and LOP::doPatternVolume>
      volumePattern(Context& ctx) const
      {
        auto& pattern = ctx.inside().pattern();
        _lop.pattern_volume(
          ctx.trial().functionSpace(),
          ctx.test().functionSpace(),
          pattern
          );
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and LOP::doPatternVolumePostSkeleton>
      volumePatternPostIntersections(Context& ctx) const
      {
        auto& pattern = ctx.inside().pattern();
        _lop.pattern_volume_post_skeleton(
          ctx.trial().functionSpace(),
          ctx.test().functionSpace(),
          pattern
          );
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and LOP::doPatternBoundary>
      boundaryPattern(Context& ctx) const
      {
        auto& pattern = ctx.inside().pattern();
        _lop.pattern_boundary(
          ctx.inside().trial().functionSpace(),
          ctx.inside().test().functionSpace(),
          pattern
          );
      }

      template<typename Context>
      std::enable_if_t<Std::to_true_type_v<Context> and LOP::doPatternSkeleton>
      skeletonPattern(Context& ctx) const
      {
        auto& pattern_sn = ctx.pattern(ctx.inside(),ctx.outside());
        auto& pattern_ns = ctx.pattern(ctx.outside(),ctx.inside());
        _lop.pattern_skeleton(
          ctx.inside().trial().functionSpace(),ctx.inside().test().functionSpace(),
          ctx.outside().trial().functionSpace(),ctx.outside().test().functionSpace(),
          pattern_sn,
          pattern_ns
          );
      }

      LocalOperatorAdapter(const LOP& lop)
        : _lop(lop)
      {}

    private:

      const LOP& _lop;

    };

  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_ADAPTER_HH
