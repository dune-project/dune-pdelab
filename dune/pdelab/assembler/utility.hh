// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_UTILITY_HH
#define DUNE_PDELAB_ASSEMBLER_UTILITY_HH

#include <functional>
#include <type_traits>
#include <tuple>

namespace Dune {
  namespace PDELab {

    namespace Impl {

      template<typename F, typename... Args>
      constexpr decltype(auto) invoke_if_possible(std::true_type,F&& f, Args&&... args)
      {
        return std::invoke(std::forward<F>(f),std::forward<Args>(args)...);
      }

      template<typename F, typename... Args>
      constexpr void invoke_if_possible(std::false_type,F&& f, Args&&... args)
      {}

      template<typename F, typename R, typename... Args>
      constexpr decltype(auto) invoke_or(std::true_type,F&& f, R&& r, Args&&... args)
      {
        return std::invoke(std::forward<F>(f),std::forward<Args>(args)...);
      }

      template<typename F, typename R, typename... Args>
      constexpr decltype(auto) invoke_or(std::false_type,F&& f, R&& r, Args&&... args)
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

    struct applyToVariadicArguments
    {
      template<typename... T>
      applyToVariadicArguments(T...)
      {}
    };

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


    template<typename Data>
    struct DataHolder
    {

    protected:

      DataHolder(Data&& data)
        : _data(std::move(data))
      {}

      Data& data()
      {
        return _data;
      }

    private:

      Data _data;

    };

    template<typename Data_>
    struct CellDataHolder
      : public DataHolder<Data_>
    {

      using Base = DataHolder<Data_>;

    public:

      using Data = Data_;

      CellDataHolder(Data&& data)
        : Base(std::forward<Data>(data))
      {}

      Data& data()
      {
        return Base::data();
      }

    };

    struct CellDataPlaceHolder
    {};

    template<typename Data>
    struct GlobalDataHolder
      : public DataHolder<Data>
    {

      using Base = DataHolder<Data>;

    public:

      using Global = Data;

      GlobalDataHolder(Data&& data)
        : Base(std::forward<Data>(data))
      {}

      Global& global()
      {
        return Base::data();
      }

    };

    struct GlobalDataPlaceHolder
    {};

    struct PlaceHolder
    {};

    namespace Impl {

      template<template<typename> typename Holder, typename Engine, typename Target>
      auto extractCellContext(PriorityTag<5>, const Engine& engine, Target& target)
        -> decltype(Holder<typename Target::template CellContext<Engine>>(typename Target::template CellContext<Engine>(engine,target)))
      {
        return Holder<typename Target::template CellContext<Engine>>(typename Target::template CellContext<Engine>(engine,target));
      }

      template<template<typename> typename Holder, typename Engine, typename Target>
      auto extractCellContext(PriorityTag<4>, const Engine& engine, Target& target)
        -> decltype(Holder<typename Target::template CellContext<Engine>>(typename Target::template CellContext<Engine>(engine)))
      {
        return Holder<typename Target::template CellContext<Engine>>(typename Target::template CellContext<Engine>(engine));
      }

      template<template<typename> typename Holder, typename Engine, typename Target>
      auto extractCellContext(PriorityTag<3>, const Engine& engine, Target& target)
        -> decltype(Holder<typename Target::template CellContext<Engine>>(typename Target::template CellContext<Engine>(target)))
      {
        return Holder<typename Target::template CellContext<Engine>>(typename Target::template CellContext<Engine>(target));
      }

      template<template<typename> typename Holder, typename Engine, typename Target>
      auto extractCellContext(PriorityTag<2>, const Engine& engine, Target& target)
        -> decltype(Holder<typename Target::template CellContext<Engine>>(typename Target::template CellContext<Engine>()))
      {
        return Holder<typename Target::template CellContext<Engine>>(typename Target::template CellContext<Engine>());
      }

      template<template<typename> typename Holder, typename Engine, typename Target>
      auto extractCellContext(PriorityTag<1>, const Engine& engine, Target& target)
        -> CellDataPlaceHolder
      {
        return {};
      }

    }

    template<template<typename> class Holder, typename Engine, typename Target>
    auto extractCellContext(const Engine& engine, Target& target)
    {
      return Impl::extractCellContext<Holder>(PriorityTag<5>{},engine,target);
    }


    namespace Impl {

      template<template<typename> typename Holder, typename Engine, typename Target>
      auto extractContext(PriorityTag<5>, const Engine& engine, Target& target)
        -> decltype(Holder<typename Target::template Context<Engine>>(typename Target::template Context<Engine>(engine,target)))
      {
        return Holder<typename Target::template Context<Engine>>(typename Target::template Context<Engine>(engine,target));
      }

      template<template<typename> typename Holder, typename Engine, typename Target>
      auto extractContext(PriorityTag<4>, const Engine& engine, Target& target)
        -> decltype(Holder<typename Target::template Context<Engine>>(typename Target::template Context<Engine>(engine)))
      {
        return Holder<typename Target::template Context<Engine>>(typename Target::template Context<Engine>(engine));
      }

      template<template<typename> typename Holder, typename Engine, typename Target>
      auto extractContext(PriorityTag<3>, const Engine& engine, Target& target)
        -> decltype(Holder<typename Target::template Context<Engine>>(typename Target::template Context<Engine>(target)))
      {
        return Holder<typename Target::template Context<Engine>>(typename Target::template Context<Engine>(target));
      }

      template<template<typename> typename Holder, typename Engine, typename Target>
      auto extractContext(PriorityTag<2>, const Engine& engine, Target& target)
        -> decltype(Holder<typename Target::template Context<Engine>>(typename Target::template Context<Engine>()))
      {
        return Holder<typename Target::template Context<Engine>>(typename Target::template Context<Engine>());
      }

      template<template<typename> typename Holder, typename Engine, typename Target>
      auto extractContext(PriorityTag<1>, const Engine& engine, Target& target)
        -> GlobalDataPlaceHolder
      {
        return {};
      }

    }

    template<template<typename> class Holder, typename Engine, typename Target>
    auto extractContext(const Engine& engine, Target& target)
    {
      return Impl::extractContext<Holder>(PriorityTag<5>{},engine,target);
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_UTILITY_HH
