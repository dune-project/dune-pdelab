// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_CONTEXT_HH
#define DUNE_PDELAB_ASSEMBLER_CONTEXT_HH

#include <type_traits>

#include <dune/pdelab/assembler/utility.hh>

namespace Dune {
  namespace PDELab {

    template<typename... Components>
    class ContextBase
      : public Components...
    {

    public:

      ContextBase(Components&&... components)
        : Components(std::forward<Components>(components))...
      {}

      template<typename Engine>
      void setup(Engine& engine)
      {
        auto setup = [&](auto&& c)
          -> decltype(c.setup(*this,engine))
          { return c.setup(*this,engine); };
        applyToVariadicArguments{invoke_if_possible_discard_return(setup,static_cast<Components&>(*this))...};
      }

      template<typename... Args>
      void bind(Args&&... args)
      {
        auto bind = [&](auto&& c)
          -> decltype(c.bind(*this,std::forward<Args>(args)...))
          { return c.bind(*this,std::forward<Args>(args)...); };
        applyToVariadicArguments{invoke_if_possible_discard_return(bind,static_cast<Components&>(*this))...};
      }

      template<typename... Args>
      void unbind(Args&&... args)
      {
        auto unbind = [&](auto&& c)
          -> decltype(c.unbind(*this,std::forward<Args>(args)...))
          { return c.unbind(*this,std::forward<Args>(args)...); };
        auto apply = [&](auto&& c) { return invoke_if_possible_discard_return(unbind,std::forward<decltype(c)>(c)); };
        applyToVariadicArgumentsWithOrder(apply,std::forward_as_tuple(static_cast<Components&>(*this)...),reverse_index_sequence_for<Components...>{});
      }

    };

    template<typename... Components>
    class CellContext
      : public ContextBase<Components...>
    {

      using Base = ContextBase<Components...>;

    public:

      using Domain = typename Base::CellDomain;

      Domain& domain()
      {
        return Base::cellDomain();
      }

      const Domain& domain() const
      {
        return Base::cellDomain();
      }

      CellContext(Components&&... components)
        : Base(std::forward<Components>(components)...)
      {}

      CellContext& cellContext()
      {
        return *this;
      }

      const CellContext& cellContext() const
      {
        return *this;
      }

    };

    template<typename... Components>
    class IntersectionContext
      : public CellContext<Components...>
    {

      using Base = CellContext<Components...>;

    public:

      using Domain = typename Base::IntersectionDomain;

      Domain& domain()
      {
        return Base::intersectionDomain();
      }

      const Domain& domain() const
      {
        return Base::intersectionDomain();
      }

      IntersectionContext(Components&&... components)
        : Base(std::forward<Components>(components)...)
      {}

      IntersectionContext& intersectionContext()
      {
        return *this;
      }

      const IntersectionContext& intersectionContext() const
      {
        return *this;
      }

    };

    template<typename... Components>
    class Context
      : public IntersectionContext<Components...>
    {

      using Base = IntersectionContext<Components...>;
      using Domain = int;

    public:

      Domain& domain() = delete;
      const Domain& domain() const = delete;

      Context(Components&&... components)
        : Base(std::forward<Components>(components)...)
      {}

    };

    template<typename... Components>
    auto makeContext(Components&&... components)
    {
      return Context<std::decay_t<Components>...>{std::forward<Components>(components)...};
    }


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_CELLDATA_HH
