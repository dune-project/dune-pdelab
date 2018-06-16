// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_ADAPTER_HH
#define DUNE_PDELAB_LOCALOPERATOR_ADAPTER_HH

#include <type_traits>

#include <dune/pdelab/common/geometrywrapper.hh>

namespace Dune {
  namespace PDELab {

    template<typename LOP>
    class LocalOperatorAdapter
    {

      template<bool enabled>
      using Call = LocalAssemblerCallSwitch<LOP,enabled>;

    public:

      template<typename Context>
      std::enable_if_t<LOP::doAlphaVolume or LOP::doLambdaVolume>
      volumeIntegral(Context& ctx) const
      {
        auto residual = ctx.residual();
        auto eg = ElementGeometry<typename Context::Entity>(ctx.entity());
        Call<LOP::doAlphaVolume and not ctx.skipVariablePart()>::alpha_volume(
          _lop,eg,
          ctx.trial().functionSpace(),ctx.argument(),
          ctx.test().functionSpace(),residual
          );
        Call<LOP::doLambdaVolume and not ctx.skipConstantPart()>::lambda_volume(
          _lop,eg,
          ctx.test().functionSpace(),residual
          );
      }

      template<typename Context>
      std::enable_if_t<LOP::doAlphaVolumePostSkeleton or LOP::doLambdaVolumePostSkeleton>
      volumeIntegralPostIntersections(Context& ctx) const
      {
        auto residual = ctx.residual();
        auto eg = ElementGeometry<typename Context::Entity>(ctx.entity());
        Call<LOP::doAlphaVolumePostSkeleton and not ctx.skipVariablePart()>::alpha_volume_post_skeleton(
          _lop,eg,
          ctx.trial().functionSpace(),ctx.argument(),
          ctx.test().functionSpace(),residual
          );
        Call<LOP::doLambdaVolumePostSkeleton and not ctx.skipConstantPart()>::lambda_volume_post_skeleton(
          _lop,eg,
          ctx.test().functionSpace(),residual
          );
      }

      template<typename Context>
      std::enable_if_t<LOP::doAlphaBoundary or LOP::doLambdaBoundary>
      boundaryIntegral(Context& ctx) const
      {
        auto residual = ctx.inside().residual();
        auto ig = IntersectionGeometry<typename Context::Domain::Intersection>(ctx.domain().intersection());
        Call<LOP::doAlphaBoundary and not ctx.skipVariablePart()>::alpha_boundary(
          _lop,ig,
          ctx.inside().trial().functionSpace(),ctx.inside().argument(),
          ctx.inside().test().functionSpace(),residual
          );
        Call<LOP::doLambdaBoundary and not ctx.skipConstantPart()>::lambda_boundary(
          _lop,ig,
          ctx.inside().test().functionSpace(),residual
          );
      }

      template<typename Context>
      std::enable_if_t<LOP::doAlphaSkeleton or LOP::doLambdaSkeleton>
      skeletonIntegral(Context& ctx) const
      {
        auto inside_residual = ctx.inside().residual();
        auto outside_residual = ctx.outside().residual();
        auto ig = IntersectionGeometry<typename Context::Domain::Intersection>(ctx.domain().intersection());
        Call<LOP::doAlphaSkeleton and not ctx.skipVariablePart()>::alpha_skeleton(
          _lop,ig,
          ctx.inside().trial().functionSpace(),ctx.inside().argument(),ctx.inside().test().functionSpace(),
          ctx.outside().trial().functionSpace(),ctx.outside().argument(),ctx.outside().test().functionSpace(),
          inside_residual,outside_residual
          );
        Call<LOP::doLambdaSkeleton and not ctx.skipConstantPart()>::lambda_skeleton(
          _lop,ig,
          ctx.inside().test().functionSpace(),ctx.outside().test().functionSpace(),
          inside_residual,outside_residual
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
