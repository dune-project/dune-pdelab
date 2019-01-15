// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_VECTORDATA_HH
#define DUNE_PDELAB_ASSEMBLER_VECTORDATA_HH

#include <cassert>
#include <type_traits>
#include <optional>

#include <dune/geometry/identitygeometry.hh>

#include <dune/pdelab/common/checks.hh>
#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/backend/common/uncachedvectorview.hh>
#include <dune/pdelab/assembler/utility.hh>
#include <dune/pdelab/assembler/localviewproxy.hh>

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

    private:

      using VectorStorage = std::conditional_t<_mode == LocalViewDataMode::read, const typename Traits::Vector, typename Traits::Vector>;

    protected:

      void setWeight(const typename Traits::Weight& weight)
      {
        _weight = weight;
      }

      void setup(VectorStorage& vector)
      {
        Context::setup();
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

      typename Traits::Container& readWriteView()
      {
        return _container;
      }

      typename Traits::Weight weight() const
      {
        return _weight;
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
      return CachedVectorData<Context,LV<std::conditional_t<_mode == LocalViewDataMode::read,const Vector, Vector>,std::decay_t<decltype(ctx.cache(Flavor{}))>>,_mode>{std::move(ctx),initial};
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
      return CachedVectorData<Context,LV<std::conditional_t<_mode == LocalViewDataMode::read,const Vector, Vector>,std::decay_t<decltype(ctx.cache(Flavor{}))>>,_mode>{std::move(ctx)};
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

      template<typename LFS>
      LocalVectorProxy<typename Residual::Container,LFS> residual(const LFS& lfs)
      {
        return {Implementation::readWriteView(),lfs,Implementation::weight()};
      }

      void setup()
      {
        Implementation::setup(Context_::engine().residual());
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
    struct CellResultData
      : public Implementation
    {

      using Context_ = typename Implementation::Context_;

      using Result = typename Implementation::Traits;
      using Residual = Result;

      typename Result::AccumulationView result()
      {
        return Implementation::accumulationView();
      }

      typename Result::AccumulationView residual()
      {
        return Implementation::accumulationView();
      }

      template<typename LFS>
      LocalVectorProxy<typename Result::Container,LFS> result(const LFS& lfs)
      {
        return {Implementation::readWriteView(),lfs,Implementation::weight()};
      }

      template<typename LFS>
      LocalVectorProxy<typename Result::Container,LFS> residual(const LFS& lfs)
      {
        return {Implementation::readWriteView(),lfs,Implementation::weight()};
      }

      void setup()
      {
        Implementation::setup(Context_::engine().result());
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

      CellResultData(Implementation&& implementation)
        : Implementation(std::move(implementation))
      {}

    };

    template<typename Implementation>
    auto cellResultData(Implementation&& implementation)
    {
      return CellResultData<Implementation>(std::move(implementation));
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

      void setup()
      {
        Implementation::setup(Context_::engine().argument());
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



    template<bool enabled, typename Implementation>
    struct CellLinearizationPointData
      : public Implementation
    {

      using Context_ = typename Implementation::Context_;

      using LinearizationPoint = typename Implementation::Traits;

      template<typename T = int>
      const typename LinearizationPoint::Container& linearizationPoint(T dummy = 0)
      {
        static_assert(Std::to_true_type<T>::value and enabled, "Calling linearizationPoint() is not allowed for linear problems!");
#if DUNE_PDELAB_ENABLE_CHECK_ASSEMBLY
        if (not Context_::engine().bindLinearizationPoint())
          DUNE_THROW(AssemblyError, "Not allowed to call linearizationPoint() for linear operators");
#endif
        return Implementation::readOnlyView();
      }

      void setup()
      {
        if (enabled and Context_::engine().bindLinearizationPoint())
          Implementation::setup(Context_::engine().linearizationPoint());
        else
          Context_::setup();
      }

      using Context_::bind;
      using Context_::unbind;

      Context_* bind(const typename Context_::Entity&, typename Context_::Index, typename Context_::Index)
      {
        if (enabled and Context_::engine().bindLinearizationPoint())
          Implementation::bind(Context_::trial().cache());
        return this;
      }

      Context_* unbind(const typename Context_::Entity&, typename Context_::Index, typename Context_::Index)
      {
        if (enabled and Context_::engine().bindLinearizationPoint())
          Implementation::unbind(Context_::trial().cache());
        return this;
      }

      CellLinearizationPointData(Implementation&& implementation)
        : Implementation(std::move(implementation))
      {}

    };

    template<bool enabled, typename Implementation>
    auto cellLinearizationPointData(std::bool_constant<enabled>, Implementation&& implementation)
    {
      return CellLinearizationPointData<enabled,Implementation>(std::move(implementation));
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_VECTORDATA_HH
