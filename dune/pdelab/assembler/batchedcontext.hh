// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_BATCHEDCONTEXT_HH
#define DUNE_PDELAB_ASSEMBLER_BATCHEDCONTEXT_HH

#include <type_traits>

#include <dune/pdelab/assembler/utility.hh>

namespace Dune {
  namespace PDELab {

    namespace Context {

      template<typename Ctx, typename BatchConfig>
      class BatchedContext
        : public Base
      {

        std::array<Ctx,BatchConfig::batchSize()> _contexts;

        static_assert(BatchConfig::batchSize() > 0, "Cannot have empty batches");

      public:

        using Context = Ctx;
        using iterator = typename std::array<Ctx,BatchConfig::batchSize()>::iterator;

        static constexpr bool batched()
        {
          return true;
        }

        static constexpr int batchSize()
        {
          return BatchConfig::batchSize();
        }

        Context& operator[](int i)
        {
          return _contexts[i];
        }

        iterator begin()
        {
          return begin(_contexts);
        }

        iterator end()
        {
          return end(_contexts);
        }

        template<typename Factory>
        BatchedContext(Factory factory)
        {
          std::generate(begin(_contexts),_end(contexts),factory);
        }

        void setup()
        {
          for (auto& ctx : _contexts)
            ctx.setup();
        }

        template<typename... Args>
        bind(Args&&... args)
        {
          for (int i = 0 ; i < batchSize() ; ++i)
            _contexts[i].bind(std::forward<Args>(args)[i]...);
        }

        template<typename... Args>
        unbind(Args&&... args)
        {
          for (int i = 0 ; i < batchSize() ; ++i)
            _contexts[i].unbind(std::forward<Args>(args)[i]...);
        }

      };

      template<typename Factory, typename BatchConfig>
      auto makeBatchedContext(Factory factory, const BatchConfig& bc)
      {
        return BatchedContext<decltype(factory()),BatchConfig>(std::move(factory));
      }

    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_BATCHEDCONTEXT_HH
