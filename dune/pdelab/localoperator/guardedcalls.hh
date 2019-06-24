// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_GUARDEDCALLS_HH
#define DUNE_PDELAB_LOCALOPERATOR_GUARDEDCALLS_HH

#include <cassert>
#include <type_traits>

#include <dune/typetree/childextraction.hh>

#include <dune/pdelab/assembler/context.hh>
#include <dune/pdelab/assembler/utility.hh>
#include <dune/pdelab/assembler/finiteelementwrapper.hh>

namespace Dune::PDELab::Experimental {

  namespace LocalOperator {

    inline constexpr auto start()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.start(ctx)) {
               return lop.start(ctx);
             };
    }

    inline constexpr auto finish()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.finish(ctx)) {
               return lop.finish(ctx);
             };
    }

    inline constexpr auto skipCell()
    {
      return [&](const auto& lop, auto& ctx, const auto& element, auto index) -> decltype(lop.skipCell(ctx,element,index)) {
               return lop.skipCell(ctx,element,index);
             };
    }

    inline constexpr auto startCell()
    {
      return [&](const auto& lop, auto& ctx, const auto& element, auto index) -> decltype(lop.startCell(ctx,element,index)) {
               return lop.startCell(ctx,element,index);
             };
    }

    inline constexpr auto finishCell()
    {
      return [&](const auto& lop, auto& ctx, const auto& element, auto index) -> decltype(lop.finishCell(ctx,element,index)) {
               return lop.finishCell(ctx,element,index);
             };
    }

    inline constexpr auto startIntersections()
    {
      return [&](const auto& lop, auto& ctx, const auto& element, auto index) -> decltype(lop.startIntersections(ctx,element,index)) {
               return lop.startIntersections(ctx,element,index);
             };
    }

    inline constexpr auto finishIntersections()
    {
      return [&](const auto& lop, auto& ctx, const auto& element, auto index) -> decltype(lop.finishIntersections(ctx,element,index)) {
               return lop.finishIntersections(ctx,element,index);
             };
    }

    inline constexpr auto volumeIntegral()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.volumeIntegral(ctx)) {
               return lop.volumeIntegral(ctx);
             };
    }

    inline constexpr auto skeletonIntegral()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.skeletonIntegral(ctx)) {
               return lop.skeletonIntegral(ctx);
             };
    }

    inline constexpr auto boundaryIntegral()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.boundaryIntegral(ctx)) {
               return lop.boundaryIntegral(ctx);
             };
    }

    inline constexpr auto volumeIntegralPostIntersections()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.volumeIntegralPostIntersections(ctx)) {
               return lop.volumeIntegralPostIntersections(ctx);
             };
    }

    inline constexpr auto volumeJacobian()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.volumeJacobian(ctx)) {
               return lop.volumeJacobian(ctx);
             };
    }

    inline constexpr auto skeletonJacobian()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.skeletonJacobian(ctx)) {
               return lop.skeletonJacobian(ctx);
             };
    }

    inline constexpr auto boundaryJacobian()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.boundaryJacobian(ctx)) {
               return lop.boundaryJacobian(ctx);
             };
    }

    inline constexpr auto volumeJacobianPostIntersections()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.volumeJacobianPostIntersections(ctx)) {
               return lop.volumeJacobianPostIntersections(ctx);
             };
    }


    inline constexpr auto volumeApplyJacobian()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.volumeApplyJacobian(ctx)) {
               return lop.volumeApplyJacobian(ctx);
             };
    }

    inline constexpr auto skeletonApplyJacobian()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.skeletonApplyJacobian(ctx)) {
               return lop.skeletonApplyJacobian(ctx);
             };
    }

    inline constexpr auto boundaryApplyJacobian()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.boundaryApplyJacobian(ctx)) {
               return lop.boundaryApplyJacobian(ctx);
             };
    }

    inline constexpr auto volumeApplyJacobianPostIntersections()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.volumeApplyJacobianPostIntersections(ctx)) {
               return lop.volumeApplyJacobianPostIntersections(ctx);
             };
    }


    inline constexpr auto volumePattern()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.volumePattern(ctx)) {
               return lop.volumePattern(ctx);
             };
    }

    inline constexpr auto skeletonPattern()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.skeletonPattern(ctx)) {
               return lop.skeletonPattern(ctx);
             };
    }

    inline constexpr auto boundaryPattern()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.boundaryPattern(ctx)) {
               return lop.boundaryPattern(ctx);
             };
    }

    inline constexpr auto volumePatternPostIntersections()
    {
      return [&](const auto& lop, auto& ctx) -> decltype(lop.volumePatternPostIntersections(ctx)) {
               return lop.volumePatternPostIntersections(ctx);
             };
    }


    template<typename LOP>
    inline constexpr auto intersectionsTwoSided()
    {
      auto call = [&](auto h) -> decltype(deduce_to_first_t<LOP,decltype(h)>::intersectionsTwoSided()) {
                    return deduce_to_first_t<LOP,decltype(h)>::intersectionsTwoSided();
                  };
      return constexpr_invoke_or(call,false,0);
    }

    template<typename LOP>
    inline constexpr auto intersectionsTwoSided(const LOP&)
    {
      return intersectionsTwoSided<LOP>();
    }

    template<typename LOP>
    inline constexpr auto disableFunctionSpaceFlavors()
    {
      constexpr auto call = [&](auto h) -> decltype(deduce_to_first_t<LOP,decltype(h)>::disableFunctionSpaceFlavors()) {
                              return deduce_to_first_t<LOP,decltype(h)>::disableFunctionSpaceFlavors();
                            };
      return constexpr_invoke_or(call,false,0);
    }

    template<typename LOP>
    inline constexpr auto disableFunctionSpaceFlavors(const LOP&)
    {
      return disableFunctionSpaceFlavors<LOP>();
    }

  } // namespace LocalOperator

} // namespace Dune::PDELab


#endif // DUNE_PDELAB_LOCALOPERATOR_GUARDEDCALLS_HH
