// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_UTILITY_HH
#define DUNE_PDELAB_ASSEMBLER_UTILITY_HH

#include <functional>
#include <type_traits>
#include <tuple>

#include <dune/common/concept.hh>

#include <dune/pdelab/constraints/common/constraintstransformation.hh>

namespace Dune {
  namespace PDELab {

    namespace Concept {

      struct IntersectionEmbedding
      {
        template<class Embedding>
        auto require(Embedding&&) -> decltype(
          Dune::Concept::requireType<typename Embedding::Intersection>()
        );
      };

      struct PossiblyNonLinear
      {
        template<class LocalOperator>
        auto require(const LocalOperator& lop) -> decltype(
          lop.isNonlinear()
        );
      };

      struct ReturnsMatrixPattern
      {
        template<typename GridOperator>
        auto require(const GridOperator& go) -> decltype(
          go.matrixPattern()
        );
      };

    }


    namespace Impl {

      template<typename F, typename... Args>
      decltype(auto) invoke_if_possible(std::true_type,F&& f, Args&&... args)
      {
        return std::invoke(std::forward<F>(f),std::forward<Args>(args)...);
      }

      template<typename F, typename... Args>
      void invoke_if_possible(std::false_type,F&& f, Args&&... args)
      {}

      template<typename F, typename R, typename... Args>
      decltype(auto) invoke_or(std::true_type,F&& f, R&& r, Args&&... args)
      {
        return std::invoke(std::forward<F>(f),std::forward<Args>(args)...);
      }

      template<typename F, typename R, typename... Args>
      decltype(auto) invoke_or(std::false_type,F&& f, R&& r, Args&&... args)
      {
        return std::forward<R>(r);
      }

      template<typename F, typename... Args>
      constexpr decltype(auto) constexpr_invoke_if_possible(std::true_type,F&& f, Args&&... args)
      {
        return f(std::forward<Args>(args)...);
      }

      template<typename F, typename... Args>
      constexpr void constexpr_invoke_if_possible(std::false_type,F&& f, Args&&... args)
      {}

      template<typename F, typename R, typename... Args>
      constexpr decltype(auto) constexpr_invoke_or(std::true_type,F&& f, R&& r, Args&&... args)
      {
        return f(std::forward<Args>(args)...);
      }

      template<typename F, typename R, typename... Args>
      constexpr decltype(auto) constexpr_invoke_or(std::false_type,F&& f, R&& r, Args&&... args)
      {
        return std::forward<R>(r);
      }

    }

    template<typename F, typename... Args>
    constexpr decltype(auto) invoke_if_possible(F&& f, Args&&... args)
    {
      return Impl::invoke_if_possible(std::is_invocable<F,decltype(std::forward<Args>(args))...>{},std::forward<F>(f),std::forward<Args>(args)...);
    }

    template<typename F, typename... Args>
    constexpr int invoke_if_possible_discard_return(F&& f, Args&&... args)
    {
      invoke_if_possible(std::forward<F>(f),std::forward<Args>(args)...);
      return 0;
    }

    template<typename F, typename R, typename... Args>
    constexpr decltype(auto) invoke_or(F&& f, R&& r, Args&&... args)
    {
      return Impl::invoke_or(std::is_invocable<F,decltype(std::forward<Args>(args))...>{},std::forward<F>(f),std::forward<R>(r),std::forward<Args>(args)...);
    }

    template<typename F, typename... Args>
    constexpr decltype(auto) constexpr_invoke_if_possible(F&& f, Args&&... args)
    {
      return Impl::constexpr_invoke_if_possible(std::is_invocable<F,decltype(std::forward<Args>(args))...>{},std::forward<F>(f),std::forward<Args>(args)...);
    }

    template<typename F, typename... Args>
    constexpr int constexpr_invoke_if_possible_discard_return(F&& f, Args&&... args)
    {
      constexpr_invoke_if_possible(std::forward<F>(f),std::forward<Args>(args)...);
      return 0;
    }

    template<typename F, typename R, typename... Args>
    constexpr decltype(auto) constexpr_invoke_or(F&& f, R&& r, Args&&... args)
    {
      return Impl::constexpr_invoke_or(std::is_invocable<F,decltype(std::forward<Args>(args))...>{},std::forward<F>(f),std::forward<R>(r),std::forward<Args>(args)...);
    }

    struct applyToVariadicArguments
    {
      template<typename... T>
      applyToVariadicArguments(T...)
      {}
    };

    template<typename T, typename U>
    struct deduce_to_first
    {
      using type = T;
    };

    template<typename T, typename U>
    using deduce_to_first_t = typename deduce_to_first<T,U>::type;

    namespace Impl {

      template<typename T, T i, T end, T increment, T... sequence>
      struct make_ascending_integer_sequence
        : public make_ascending_integer_sequence<T,i + increment,end,increment,sequence...,i>
      {};

      template<typename T, T end, T increment, T... sequence>
      struct make_ascending_integer_sequence<T,end,end,increment,sequence...>
      {
        using type = std::integer_sequence<T,sequence...>;
      };

      template<typename T, T i, T end, T increment, T... sequence>
      struct make_descending_integer_sequence
        : public make_descending_integer_sequence<T,i + increment,end,increment,sequence...,i>
      {};

      template<typename T, T end, T increment, T... sequence>
      struct make_descending_integer_sequence<T,end,end,increment,sequence...>
      {
        using type = std::integer_sequence<T,sequence...,end>;
      };

      template<typename T, T start, T end, std::make_signed_t<T> increment = 1>
      struct make_general_integer_sequence
        : public std::conditional_t<
                   (increment > 0),
                   make_ascending_integer_sequence<T,start,end,increment>,
                   make_descending_integer_sequence<T,start + increment,end,increment>
                   >
      {};

    }

    template<typename T, T start, T end, std::make_signed_t<T> increment = 1>
    using make_general_integer_sequence = typename Impl::make_general_integer_sequence<T,start,end,increment>::type;

    template<std::size_t start, std::size_t end, std::make_signed_t<std::size_t> increment = 1>
    using make_general_index_sequence = make_general_integer_sequence<std::size_t,start,end,increment>;

    template<typename... T>
    using reverse_index_sequence_for = make_general_index_sequence<sizeof...(T),0,-1>;

    template<typename F, typename Tuple, std::size_t... indices>
    void applyToVariadicArgumentsWithOrder(F&& f, Tuple&& tuple, std::index_sequence<indices...>)
    {
      applyToVariadicArguments{f(std::get<indices>(std::forward<Tuple>(tuple)))...};
    }

    struct PlaceHolder
    {};

    namespace Impl {

      template<typename Target, typename Context>
      auto extractCellContext(PriorityTag<2>, Target& target, Context&& ctx)
        -> decltype(typename Target::template CellContext<Context>(std::move(ctx)))
      {
        return typename Target::template CellContext<Context>(std::move(ctx));
      }

      template<typename Target, typename Context>
      auto extractCellContext(PriorityTag<1>, Target& target, Context&& ctx)
        -> Context&&
      {
        return std::move(ctx);
      }

    }

    template<typename Target, typename Context>
    auto extractCellContext(Target& target, Context&& ctx)
    {
      return Impl::extractCellContext(PriorityTag<2>{},target, std::move(ctx));
    }


    namespace Impl {

      template<typename Target, typename Context>
      auto extractContext(PriorityTag<2>, Target& target, Context&& ctx)
        -> decltype(typename Target::template Context<Context>(std::move(ctx)))
      {
        return typename Target::template Context<Context>(std::move(ctx));
      }

      template<typename Target, typename Context>
      auto extractContext(PriorityTag<1>, Target& target, Context&& ctx)
        -> Context&&
      {
        return std::move(ctx);
      }

    }

    template<typename Target, typename Context>
    auto extractContext(Target& target, Context&& ctx)
    {
      return Impl::extractContext(PriorityTag<2>{},target, std::move(ctx));
    }



    template<typename Pattern, typename RI, typename CI>
    typename std::enable_if<
      std::is_same<RI,CI>::value
      >::type
    add_diagonal_entry(Pattern& pattern, const RI& ri, const CI& ci)
    {
      if (ri == ci)
        pattern.add_link(ri,ci);
    }

    template<typename Pattern, typename RI, typename CI>
    typename std::enable_if<
      !std::is_same<RI,CI>::value
      >::type
    add_diagonal_entry(Pattern& pattern, const RI& ri, const CI& ci)
    {}

    template<typename GFSV, typename GC, typename C>
    void set_trivial_rows(const GFSV& gfsv, GC& globalcontainer, const C& c)
    {
      typedef typename C::const_iterator global_row_iterator;
      for (global_row_iterator cit = c.begin(); cit != c.end(); ++cit)
        globalcontainer.clear_row(cit->first,1);
    }

    template<typename GFSV, typename GC>
    void set_trivial_rows(const GFSV& gfsv, GC& globalcontainer, const EmptyTransformation& c)
    {
    }

    template<typename LOP>
    constexpr std::enable_if_t<models<Concept::PossiblyNonLinear,LOP>(),bool>
    isNonlinear(const LOP& lop)
    {
      return lop.isNonlinear();
    }

    template<typename LOP>
    constexpr std::enable_if_t<not models<Concept::PossiblyNonLinear,LOP>(),bool>
    isNonlinear(const LOP& lop)
    {
      return false;
    }

    enum class LocalViewDataMode { read, write, accumulate, readWrite, readAccumulate };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_UTILITY_HH
