// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_GUARDEDCALLS_HH
#define DUNE_PDELAB_LOCALOPERATOR_GUARDEDCALLS_HH

#include <cassert>
#include <type_traits>

#include <dune/pdelab/assembler/context.hh>
#include <dune/pdelab/assembler/utility.hh>

namespace Dune {
  namespace PDELab {

    namespace LocalOperator {

      namespace Impl {

        template<typename LFS, typename = void>
        struct extract_basis_traits
        {
          using type = typename LFS::Traits::FiniteElement::Traits::LocalBasisType::Traits;
        };

        template<typename LFS>
        struct extract_basis_traits<LFS,std::void_t<TypeTree::Child<LFS,0>>>
          : public extract_basis_traits<TypeTree::Child<LFS,0>>
        {};

        template<typename T>
        struct extract_local_test_space
        {
          using type = typename T::TestLocalSpace;
        };

        template<typename Ctx>
        struct extract_local_test_space<Context::Context<Ctx>>
        {
          using type = typename Ctx::Test::FunctionSpace;
        };

        template<typename Ctx>
        struct extract_local_test_space<Context::CellContext<Ctx>>
        {
          using type = typename Ctx::Test::FunctionSpace;
        };

        template<typename Ctx>
        struct extract_local_test_space<Context::IntersectionContext<Ctx>>
        {
          using type = typename Ctx::Test::FunctionSpace;
        };

      }

      template<typename T>
      using TestSpace = typename Impl::extract_local_test_space<T>::type;

      template<typename T, std::size_t... I>
      using Range = typename Impl::extract_basis_traits<TypeTree::Child<TestSpace<T>,I...>>::type::RangeType;

      template<typename T, std::size_t... I>
      using RangeField = typename Impl::extract_basis_traits<TypeTree::Child<TestSpace<T>,I...>>::type::RangeFieldType;

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

      template<typename LOP>
      inline constexpr auto intersectionsTwoSided()
      {
        auto call = [&](const auto* lop) -> decltype(lop->intersectionsTwoSided()) {
          return std::remove_pointer_t<std::decay_t<decltype(lop)>>::intersectionsTwoSided();
        };
        return invoke_or(call,false,nullptr);
      }

      template<typename LOP>
      inline constexpr auto intersectionsTwoSided(const LOP&)
      {
        return intersectionsTwoSided<LOP>();
      }

    }

  } // namespace PDELab
} // namespace Dune


#endif // DUNE_PDELAB_LOCALOPERATOR_GUARDEDCALLS_HH
