// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_VECTORDATA_HH
#define DUNE_PDELAB_ASSEMBLER_VECTORDATA_HH

#include <cassert>
#include <type_traits>
#include <optional>

#include <dune/geometry/identitygeometry.hh>

#include <dune/pdelab/common/checks.hh>
#include <dune/pdelab/common/exceptions.hh>
#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/common/typetraits.hh>
#include <dune/pdelab/backend/common/aliasedvectorview.hh>
#include <dune/pdelab/backend/common/uncachedvectorview.hh>
#include <dune/pdelab/assembler/utility.hh>
#include <dune/pdelab/assembler/localviewproxy.hh>

namespace Dune {
  namespace PDELab {

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
        using View             = Container;
        using AccumulationView = typename Container::WeightedAccumulationView;
        using Weight           = typename AccumulationView::weight_type;
        using value_type       = typename LV::value_type;

        static constexpr LocalViewDataMode dataMode = _mode;
      };

    private:

      using VectorStorage = std::conditional_t<_mode == LocalViewDataMode::read, const typename Traits::Vector, typename Traits::Vector>;

    protected:

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

      typename Traits::AccumulationView accumulationView(typename Traits::Weight weight)
      {
        return _container.weightedAccumulationView(weight);
      }

      const typename Traits::Container& readOnlyView()
      {
        return _container;
      }

      typename Traits::Container& readWriteView()
      {
        return _container;
      }

    public:

      CachedVectorData(Context&& ctx, typename Traits::value_type initial = typename Traits::value_type(0))
        : Context(std::move(ctx))
        , _initial(initial)
      {}

    private:

      typename Traits::LocalView _local_view;
      typename Traits::Container _container;
      typename Traits::value_type _initial;

    };


    template<typename Context, typename LV, LocalViewDataMode _mode>
    class AliasedVectorData
      : public Context
    {

    protected:

      using Context_ = Context;

      struct Traits
      {
        using LocalView        = LV;
        using Vector           = typename LocalView::Container;
        using View             = LocalView;
        using AccumulationView = LocalView;
        using Weight           = typename AccumulationView::weight_type;
        using value_type       = typename LV::value_type;

        static constexpr LocalViewDataMode dataMode = _mode;
      };

    private:

      using VectorStorage = std::conditional_t<_mode == LocalViewDataMode::read, const typename Traits::Vector, typename Traits::Vector>;

    protected:

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
      }

      template<typename LFSCache>
      void unbind(LFSCache& lfs_cache)
      {
        if constexpr (
          _mode == LocalViewDataMode::write or
          _mode == LocalViewDataMode::readWrite or
          _mode == LocalViewDataMode::accumulate or
          _mode == LocalViewDataMode::readAccumulate
          )
        {
          _local_view.commit();
        }
        _local_view.unbind();
      }

      typename Traits::AccumulationView accumulationView(typename Traits::Weight weight)
      {
        DUNE_THROW(NotImplemented,"accumulationView() not implemented");
      }

      const typename Traits::LocalView& readOnlyView()
      {
        return _local_view;
      }

      typename Traits::LocalView& readWriteView()
      {
        return _local_view;
      }

    public:

      AliasedVectorData(Context&& ctx, typename Traits::value_type initial = typename Traits::value_type(0))
        : Context(std::move(ctx))
      {}

    private:

      typename Traits::LocalView _local_view;

    };

    template<
      typename Vector,
      typename Flavor,
      LocalViewDataMode _mode,
      typename Field,
      typename Context
      >
    auto vectorData(const Field& initial, Context&& ctx)
    {
      if constexpr (Context::fastDG())
      {
        using VectorView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstAliasedVectorView<const Vector,std::decay_t<decltype(ctx.cache(Flavor{}))>>,
          AliasedVectorView<Vector,std::decay_t<decltype(ctx.cache(Flavor{}))>>
          >;
        return AliasedVectorData<Context,VectorView,_mode>{std::move(ctx),initial};
      }
      else
      {
        using VectorView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstUncachedVectorView<const Vector,std::decay_t<decltype(ctx.cache(Flavor{}))>>,
          UncachedVectorView<Vector,std::decay_t<decltype(ctx.cache(Flavor{}))>>
          >;
        return CachedVectorData<Context,VectorView,_mode>{std::move(ctx),initial};
      }
    }

    template<
      typename Vector,
      typename Flavor,
      LocalViewDataMode _mode,
      typename Context
      >
    auto vectorData(Context&& ctx)
    {
      if constexpr (Context::fastDG())
      {
        using VectorView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstAliasedVectorView<const Vector,std::decay_t<decltype(ctx.cache(Flavor{}))>>,
          AliasedVectorView<Vector,std::decay_t<decltype(ctx.cache(Flavor{}))>>
          >;
        return AliasedVectorData<Context,VectorView,_mode>{std::move(ctx)};
      }
      else
      {
        using VectorView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstUncachedVectorView<const Vector,std::decay_t<decltype(ctx.cache(Flavor{}))>>,
          UncachedVectorView<Vector,std::decay_t<decltype(ctx.cache(Flavor{}))>>
          >;
        return CachedVectorData<Context,VectorView,_mode>{std::move(ctx)};
      }
    }

    template<typename Implementation>
    struct CellResidualData
      : public Implementation
    {

      using Context_ = typename Implementation::Context_;

      using Residual = typename Implementation::Traits;

      typename Residual::AccumulationView residual()
      {
        return Implementation::accumulationView(Context_::engine().weight());
      }

      template<typename LFS>
      LocalVectorProxy<typename Residual::View,LFS,Residual::dataMode> residual(const LFS& lfs)
      {
        return {Implementation::readWriteView(),lfs,Context_::engine().weight()};
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
    struct CellTimeResidualData
      : public Implementation
    {

      using Context_ = typename Implementation::Context_;

      using Residual = typename Implementation::Traits;

      typename Residual::AccumulationView timeResidual()
      {
        return Implementation::accumulationView(Context_::engine().timeWeight());
      }

      template<typename LFS>
      LocalVectorProxy<typename Residual::View,LFS,Residual::dataMode> residual(const LFS& lfs)
      {
        return {Implementation::readWriteView(),lfs,Context_::engine().timeWeight()};
      }

      void setup()
      {
        Implementation::setup(Context_::engine().timeResidual());
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

      CellTimeResidualData(Implementation&& implementation)
        : Implementation(std::move(implementation))
      {}

    };

    template<
      typename Vector,
      typename Flavor,
      typename Context
      >
    auto cellTimeResidualData(std::true_type, Context&& context)
    {
      using Implementation = std::decay_t<
        decltype(vectorData<Vector,Flavor,LocalViewDataMode::accumulate>(std::move(context)))
        >;
      return CellTimeResidualData<Implementation>(
        vectorData<Vector,Flavor,LocalViewDataMode::accumulate>(std::move(context))
        );
    }

    template<
      typename Vector,
      typename Flavor,
      typename Context
      >
    auto cellTimeResidualData(std::false_type, Context&& context)
    {
      return std::move(context);
    }


    template<typename Implementation>
    struct CellResultData
      : public Implementation
    {

    protected:

      using ResultImplementation = Implementation;

    public:

      using Context_ = typename Implementation::Context_;

      using Result = typename Implementation::Traits;
      using Residual = Result;

      typename Result::AccumulationView result()
      {
        return Implementation::accumulationView(Context_::engine().weight());
      }

      typename Result::AccumulationView residual()
      {
        return Implementation::accumulationView(Context_::engine().weight());
      }

      template<typename LFS>
      LocalVectorProxy<typename Result::View,LFS,Result::dataMode> result(const LFS& lfs)
      {
        return {Implementation::readWriteView(),lfs,Context_::engine().weight()};
      }

      template<typename LFS>
      LocalVectorProxy<typename Result::View,LFS,Result::dataMode> residual(const LFS& lfs)
      {
        return {Implementation::readWriteView(),lfs,Context_::engine().weight()};
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


    template<typename Context>
    struct CellTimeResultData
      : public Context
    {

      using Context_ = Context;
      using Result   = typename Context_::Result;

      typename Result::AccumulationView timeResult()
      {
        return Context_::ResultImplementation::accumulationView(Context_::engine().timeWeight());
      }

      typename Result::AccumulationView timeResidual()
      {
        return timeResult();
      }

      template<typename LFS>
      LocalVectorProxy<typename Result::View,LFS,Result::dataMode> timeResult(const LFS& lfs)
      {
        return {Context_::ResultImplementation::readWriteView(),lfs,Context_::engine().timeWeight()};
      }

      template<typename LFS>
      LocalVectorProxy<typename Result::View,LFS,Result::dataMode> timeResidual(const LFS& lfs)
      {
        return {Context_::ResultImplementation::readWriteView(),lfs,Context_::engine().timeWeight()};
      }

      CellTimeResultData(Context&& context)
        : Context(std::move(context))
      {}

    };

    template<typename Context>
    auto cellTimeResultData(std::true_type, Context&& context)
    {
      return CellTimeResultData<Context>(std::move(context));
    }

    template<typename Context>
    auto cellTimeResultData(std::false_type, Context&& context)
    {
      return std::move(context);
    }

    template<typename Implementation>
    struct CellCoefficientData
      : public Implementation
    {

      using Context_ = typename Implementation::Context_;

      using Coefficient = typename Implementation::Traits;

      const typename Coefficient::View& coefficient()
      {
        return Implementation::readOnlyView();
      }

      template<typename LFS>
      LocalVectorProxy<const typename Coefficient::View,LFS,Coefficient::dataMode> coefficient(const LFS& lfs)
      {
        return {Implementation::readOnlyView(),lfs};
      }

      void setup()
      {
        Implementation::setup(Context_::engine().coefficient());
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

      CellCoefficientData(Implementation&& implementation)
        : Implementation(std::move(implementation))
      {}

    };

    template<typename Implementation>
    auto cellCoefficientData(Implementation&& implementation)
    {
      return CellCoefficientData<Implementation>(std::move(implementation));
    }



    template<bool enabled, typename Implementation>
    struct CellLinearizationPointData
      : public Implementation
    {

      using Context_ = typename Implementation::Context_;

      using LinearizationPoint = typename Implementation::Traits;

      template<typename T = int>
      const typename LinearizationPoint::View& linearizationPoint(T dummy = 0)
      {
        static_assert(Std::to_true_type_v<T> and enabled, "Calling linearizationPoint() is not allowed for linear problems!");
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
