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

    template<typename Context, typename LV, LocalViewDataMode _mode>
    class CachedVectorData
      : public Context
    {

    protected:

      using Context_ = Context;

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

      void setup(typename Traits::Vector& vector)
      {
        _local_view.attach(vector);
      }

      using Context::bind;
      using Context::unbind;

      template<typename LFSCache>
      void bind(LFSCache& lfs_cache)
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

      template<typename LFSCache>
      void unbind(LFSCache& lfs_cache)
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

      typename Traits::AccumulationView accumulationView()
      {
        return _container.weightedAccumulationView(_weight);
      }

      const typename Traits::Container& readOnlyView()
      {
        return _container;
      }

    public:

      CachedVectorData(Context&& ctx, typename Traits::value_type initial = typename Traits::value_type(0))
        : Context(std::move(ctx))
        , _weight(1.0)
        , _initial(initial)
      {}

    private:

      typename Traits::LocalView _local_view;
      typename Traits::Container _container;
      typename Traits::Weight _weight;
      typename Traits::value_type _initial;

    };


    template<
      template<typename,typename> typename LV,
      typename Vector,
      typename Flavor,
      LocalViewDataMode _mode,
      typename Field,
      typename Context
      >
    auto cachedVectorData(const Field& initial, Context&& ctx)
    {
      return CachedVectorData<Context,LV<Vector,std::decay_t<decltype(ctx.cache(Flavor{}))>>,_mode>{std::move(ctx),initial};
    }

    template<
      template<typename,typename> typename LV,
      typename Vector,
      typename Flavor,
      LocalViewDataMode _mode,
      typename Context
      >
    auto cachedVectorData(Context&& ctx)
    {
      return CachedVectorData<Context,LV<Vector,std::decay_t<decltype(ctx.cache(Flavor{}))>>,_mode>{std::move(ctx)};
    }

    template<typename Implementation>
    struct CellResidualData
      : public Implementation
    {

      using Context_ = typename Implementation::Context_;

      using Residual = typename Implementation::Traits;

      typename Residual::AccumulationView residual()
      {
        return Implementation::accumulationView();
      }

      Context_* setup()
      {
        Implementation::setup(Context_::engine().residual());
        return this;
      }

      using Context_::bind;
      using Context_::unbind;

      Context_* bind(const typename Context_::Entity&, typename Context_::Index, typename Context_::Index)
      {
        Implementation::bind(Context_::test().cache());
        return this;
      }

      Context_* unbind(const typename Context_::Entity&, typename Context_::Index, typename Context_::Index)
      {
        Implementation::unbind(Context_::test().cache());
        return this;
      }

      CellResidualData(Implementation&& implementation)
        : Implementation(std::move(implementation))
      {}

    };

    template<typename Implementation>
    auto cellResidualData(Implementation&& implementation)
    {
      return CellResidualData<Implementation>(std::move(implementation));
    }

    template<typename Implementation>
    struct CellArgumentData
      : public Implementation
    {

      using Context_ = typename Implementation::Context_;

      using Argument = typename Implementation::Traits;

      const typename Argument::Container& argument()
      {
        return Implementation::readOnlyView();
      }

      Context_* setup()
      {
        Implementation::setup(Context_::engine().argument());
        return this;
      }

      using Context_::bind;
      using Context_::unbind;

      Context_* bind(const typename Context_::Entity&, typename Context_::Index, typename Context_::Index)
      {
        Implementation::bind(Context_::trial().cache());
        return this;
      }

      Context_* unbind(const typename Context_::Entity&, typename Context_::Index, typename Context_::Index)
      {
        Implementation::unbind(Context_::trial().cache());
        return this;
      }

      CellArgumentData(Implementation&& implementation)
        : Implementation(std::move(implementation))
      {}

    };

    template<typename Implementation>
    auto cellArgumentData(Implementation&& implementation)
    {
      return CellArgumentData<Implementation>(std::move(implementation));
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_VECTORDATA_HH
