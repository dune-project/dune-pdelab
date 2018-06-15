// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_VECTORDATA_HH
#define DUNE_PDELAB_ASSEMBLER_VECTORDATA_HH

#include <cassert>
#include <type_traits>
#include <optional>

#include <dune/geometry/identitygeometry.hh>

#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/backend/common/uncachedvectorview.hh>
#include <dune/pdelab/assembler/utility.hh>

namespace Dune {
  namespace PDELab {

    enum class LocalViewDataMode { read, write, accumulate, readWrite, readAccumulate };

    template<typename LV, LocalViewDataMode _mode>
    class CachedVectorData
    {

    protected:

      struct Traits
      {
        using LocalView        = LV;
        using Vector           = typename LocalView::Container;
        using Container        = LocalVector<typename LV::value_type>;
        using AccumulationView = typename Container::WeightedAccumulationView;
        using Weight           = typename AccumulationView::weight_type;
        using value_type       = typename LV::value_type;
      };

      void setWeight(const typename Traits::Weight& weight)
      {
        _weight = weight;
      }

      template<typename Context, typename CellContext, typename Engine>
      void setup(Context& ctx, CellContext& cell_ctx, Engine& engine, typename Traits::Vector& vector)
      {
        _local_view.attach(vector);
      }

      template<typename Context, typename CellContext, typename LFSCache>
      void bind(Context& ctx, CellContext& cell_ctx, const LFSCache& lfs_cache)
      {
        _local_view.bind(lfs_cache);

        if constexpr (_mode == LocalViewDataMode::read or _mode == LocalViewDataMode::readWrite or _mode == LocalViewDataMode::readAccumulate)
          {
            _container.resize(lfs_cache.size());
            _local_view.read(_container);
          }
        else if constexpr (_mode == LocalViewDataMode::write or _mode == LocalViewDataMode::accumulate)
          {
            _container.assign(lfs_cache.size(),_initial);
          }
      }

      typename Traits::AccumulationView accumulationView()
      {
        return _container.weightedAccumulationView(_weight);
      }

      const typename Traits::Container& readOnlyView()
      {
        return _container;
      }

    public:

      CachedVectorData(typename Traits::value_type initial = typename Traits::value_type(0))
        : _weight(1.0)
        , _initial(initial)
      {}

      template<typename Context, typename CellContext, typename Entity, typename Index>
      void unbind(Context& ctx, CellContext& cell_ctx, const Entity&, Index, Index)
      {
        if constexpr (_mode == LocalViewDataMode::write or _mode == LocalViewDataMode::readWrite)
          {
            _local_view.write(_container);
            _local_view.commit();
          }
        else if constexpr (_mode == LocalViewDataMode::accumulate or _mode == LocalViewDataMode::readAccumulate)
          {
            _local_view.add(_container);
            _local_view.commit();
          }

        _local_view.unbind();
      }

    private:

      typename Traits::LocalView _local_view;
      typename Traits::Container _container;
      typename Traits::Weight _weight;
      typename Traits::value_type _initial;

    };


    template<typename Implementation>
    struct CellResidualData
      : public Implementation
    {

      using Residual = typename Implementation::Traits;

      typename Residual::AccumulationView residual()
      {
        return Implementation::accumulationView();
      }

      template<typename Context, typename CellContext, typename Engine>
      void setup(Context& ctx, CellContext& cell_ctx, Engine& engine)
      {
        Implementation::setup(ctx,cell_ctx,engine,engine.residual());
      }

      template<typename Context, typename CellContext, typename Entity, typename Index>
      void bind(Context& ctx, CellContext& cell_ctx, const Entity&, Index, Index)
      {
        Implementation::bind(ctx,cell_ctx,cell_ctx.test().cache());
      }

      CellResidualData(Implementation&& implementation)
        : Implementation(std::move(implementation))
      {}

    };

    template<typename Implementation>
    auto cellResidualData(Implementation&& implementation)
    {
      return CellResidualData<Implementation>(std::forward<Implementation>(implementation));
    }

    template<typename Implementation>
    struct CellArgumentData
      : public Implementation
    {

      using Argument = typename Implementation::Traits;

      const typename Argument::Container& argument()
      {
        return Implementation::readOnlyView();
      }

      template<typename Context, typename CellContext, typename Engine>
      void setup(Context& ctx, CellContext& cell_ctx, Engine& engine)
      {
        Implementation::setup(ctx,cell_ctx,engine,engine.argument());
      }

      template<typename Context, typename CellContext, typename Entity, typename Index>
      void bind(Context& ctx, CellContext& cell_ctx, const Entity&, Index, Index)
      {
        Implementation::bind(ctx,cell_ctx,cell_ctx.trial().cache());
      }

      CellArgumentData(Implementation&& implementation)
        : Implementation(std::move(implementation))
      {}

    };

    template<typename Implementation>
    auto cellArgumentData(Implementation&& implementation)
    {
      return CellArgumentData<Implementation>(std::forward<Implementation>(implementation));
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_VECTORDATA_HH
