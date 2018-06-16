// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_CONTEXT_HH
#define DUNE_PDELAB_ASSEMBLER_CONTEXT_HH

#include <type_traits>

#include <dune/pdelab/assembler/utility.hh>

namespace Dune {
  namespace PDELab {

    namespace Context {

      struct RootContext
      {

        constexpr RootContext* setup()
        {
          return this;
        }

        template<typename... Args>
        constexpr RootContext* bind(Args&&... args)
        {
          return this;
        }

        template<typename... Args>
        constexpr RootContext* unbind(Args&&... args)
        {
          return this;
        }

      };


      template<typename Context>
      std::enable_if_t<std::is_same<Context,RootContext>::value> setup(Context& ctx)
      {}

      template<typename Context>
      std::enable_if_t<not std::is_same<Context,RootContext>::value> setup(Context& ctx)
      {
        setup(static_cast<decltype(*ctx.setup())>(ctx));
        ctx.setup();
      }

      template<typename Context, typename... Args>
      std::enable_if_t<std::is_same<Context,RootContext>::value> bind(Context& ctx, Args&&... args)
      {}

      template<typename Context, typename... Args>
      std::enable_if_t<not std::is_same<Context,RootContext>::value> bind(Context& ctx, Args&&... args)
      {
        bind(static_cast<decltype(*ctx.bind(std::forward<Args>(args)...))>(ctx),std::forward<Args>(args)...);
        ctx.bind(std::forward<Args>(args)...);
      }

      template<typename Context, typename... Args>
      std::enable_if_t<std::is_same<Context,RootContext>::value> unbind(Context& ctx, Args&&... args)
      {}

      template<typename Context, typename... Args>
      std::enable_if_t<not std::is_same<Context,RootContext>::value> unbind(Context& ctx, Args&&... args)
      {
        unbind(*ctx.unbind(std::forward<Args>(args)...),std::forward<Args>(args)...);
      }

      template<typename Context>
      class CellContext
        : public Context
      {

      public:

        using Domain = typename Context::CellDomain;

        Domain domain() const
        {
          return Context::cellDomain();
        }

        CellContext(Context&& ctx)
          : Context(std::move(ctx))
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

      template<typename Context>
      class IntersectionContext
        : public Context
      {

      public:

        using Domain = typename Context::IntersectionDomain;

        Domain domain() const
        {
          return Context::intersectionDomain();
        }

        IntersectionContext(Context&& ctx)
          : Context(std::move(ctx))
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

      template<typename Base>
      class Context
        : public Base
      {

        using Domain = int;

      public:

        Domain& domain() = delete;
        const Domain& domain() const = delete;

        Context(Base&& base)
          : Base(std::move(base))
        {}

        void setup()
        {
          Dune::PDELab::Context::setup(*static_cast<Base*>(this));
        }

        template<typename... Args>
        std::enable_if_t<sizeof...(Args) != 6> bind(Args&&... args)
        {
          Dune::PDELab::Context::bind(*static_cast<Base*>(this),std::forward<Args>(args)...);
        }

        template<typename IntersectionType, typename Intersection, typename Index, typename Entity>
        void bind(IntersectionType type, const Intersection& is, Index intersection_index, const Entity& entity, Index entity_index, Index unique_index)
        {
          Dune::PDELab::Context::bind(*static_cast<Base*>(this),type,is,intersection_index,entity,entity_index,unique_index);
        }

        template<typename... Args>
        std::enable_if_t<sizeof...(Args) != 6> unbind(Args&&... args)
        {
          Dune::PDELab::Context::unbind(*static_cast<Base*>(this),std::forward<Args>(args)...);
        }

        template<typename IntersectionType, typename Intersection, typename Index, typename Entity>
        void unbind(IntersectionType type, const Intersection& is, Index intersection_index, const Entity& entity, Index entity_index, Index unique_index)
        {
          Dune::PDELab::Context::unbind(*static_cast<Base*>(this),type,is,intersection_index,entity,entity_index,unique_index);
        }

      };

      template<typename Context_>
      auto makeContext(Context_&& ctx)
      {
        return Context<IntersectionContext<CellContext<Context_>>>{{{std::move(ctx)}}};
      }

    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_CELLDATA_HH
