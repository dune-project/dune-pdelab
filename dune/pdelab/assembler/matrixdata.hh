// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_MATRIXDATA_HH
#define DUNE_PDELAB_ASSEMBLER_MATRIXDATA_HH

#include <cassert>
#include <type_traits>
#include <optional>

#include <dune/geometry/identitygeometry.hh>

#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/backend/common/aliasedmatrixview.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/assembler/utility.hh>
#include <dune/pdelab/assembler/vectordata.hh>
#include <dune/pdelab/gridfunctionspace/flavor.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>

namespace Dune {
  namespace PDELab {

    template<typename LM>
    struct CachedMatrixDataTraits
    {
      using LocalView        = LM;
      using Matrix           = typename LocalView::Container;
      using Container        = LocalMatrix<typename LM::value_type>;
      using View             = Container;
      using AccumulationView = typename Container::WeightedAccumulationView;
      using Weight           = typename AccumulationView::weight_type;
      using value_type       = typename LM::value_type;
    };

    template<typename LM>
    struct AliasedMatrixDataTraits
    {
      using LocalView        = LM;
      using Matrix           = typename LocalView::Container;
      using View             = LocalView;
      using AccumulationView = LocalView;
      using Weight           = typename AccumulationView::weight_type;
      using value_type       = typename LM::value_type;
    };

    template<typename Context, typename LM, LocalViewDataMode _mode>
    class CachedMatrixData
      : public Context
    {

    protected:

      using Context_ = Context;

      using Traits = CachedMatrixDataTraits<LM>;

      void setup(typename Traits::Matrix& matrix)
      {
        Context::setup();
        _local_view.attach(matrix);
      }

      using Context::bind;
      using Context::unbind;

      template<typename TestCache, typename TrialCache>
      void bind(TestCache& test_cache, TrialCache& trial_cache)
      {
        _local_view.bind(test_cache,trial_cache);

        if constexpr (_mode == LocalViewDataMode::read or _mode == LocalViewDataMode::readWrite or _mode == LocalViewDataMode::readAccumulate)
          {
            _container.resize(test_cache.size(),trial_cache.size());
            _local_view.read(_container);
          }
        else if constexpr (_mode == LocalViewDataMode::write or _mode == LocalViewDataMode::accumulate)
          {
            _container.assign(test_cache.size(),trial_cache.size(),_initial);
          }
      }

      template<typename TestCache, typename TrialCache>
      void unbind(TestCache& test_cache, TrialCache& trial_cache)
      {
        if constexpr (_mode == LocalViewDataMode::write or _mode == LocalViewDataMode::readWrite)
          {
            Context::engine().scatterJacobian(_container,_local_view);
            _local_view.commit();
          }
        else if constexpr (_mode == LocalViewDataMode::accumulate or _mode == LocalViewDataMode::readAccumulate)
          {
            Context::engine().scatterJacobian(_container,_local_view);
            _local_view.commit();
          }

        _local_view.unbind();
      }

      typename Traits::AccumulationView accumulationView(typename Traits::Weight weight)
      {
        return _container.weightedAccumulationView(weight);
      }

      const typename Traits::View& readOnlyView()
      {
        return _container;
      }

      typename Traits::View& readWriteView()
      {
        return _container;
      }

    public:

      CachedMatrixData(Context&& ctx, typename Traits::value_type initial = typename Traits::value_type(0))
        : Context(std::move(ctx))
        , _initial(initial)
      {}

    private:

      typename Traits::LocalView _local_view;
      typename Traits::Container _container;
      typename Traits::value_type _initial;

    };


    template<typename Context, typename LM, LocalViewDataMode _mode>
    class AliasedMatrixData
      : public Context
    {

    protected:

      using Context_ = Context;

      using Traits = AliasedMatrixDataTraits<LM>;

      void setup(typename Traits::Matrix& matrix)
      {
        Context::setup();
        _local_view.attach(matrix);
      }

      using Context::bind;
      using Context::unbind;

      template<typename TestCache, typename TrialCache>
      void bind(TestCache& test_cache, TrialCache& trial_cache)
      {
        _local_view.bind(test_cache,trial_cache);
      }

      template<typename TestCache, typename TrialCache>
      void unbind(TestCache& test_cache, TrialCache& trial_cache)
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

      AliasedMatrixData(Context&& ctx, typename Traits::value_type initial = typename Traits::value_type(0))
        : Context(std::move(ctx))
      {}

    private:

      typename Traits::LocalView _local_view;

    };



    template<
      typename Matrix,
      LocalViewDataMode _mode,
      typename Field,
      typename Context>
    auto matrixData(const Field& initial, Context&& ctx)
    {
      if (Context::fastDG())
      {
        using MatrixView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstAliasedMatrixView<
            const Matrix,
            std::decay_t<decltype(ctx.cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.cache(Flavor::Trial{}))>
            >,
          AliasedMatrixView<
            Matrix,
            std::decay_t<decltype(ctx.cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.cache(Flavor::Trial{}))>
            >
          >;
        return AliasedMatrixData<Context,MatrixView,_mode>{std::move(ctx),initial};
      }
      else
      {
        using MatrixView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstUncachedMatrixView<
            const Matrix,
            std::decay_t<decltype(ctx.cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.cache(Flavor::Trial{}))>
            >,
          UncachedMatrixView<
            Matrix,
            std::decay_t<decltype(ctx.cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.cache(Flavor::Trial{}))>
            >
          >;
        return CachedMatrixData<Context,MatrixView,_mode>{std::move(ctx),initial};
      }
    }

    template<
      typename Matrix,
      LocalViewDataMode _mode,
      typename Context>
    auto matrixData(Context&& ctx)
    {
      if constexpr (Context::fastDG())
      {
        using MatrixView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstAliasedMatrixView<
            const Matrix,
            std::decay_t<decltype(ctx.cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.cache(Flavor::Trial{}))>
            >,
          AliasedMatrixView<
            Matrix,
            std::decay_t<decltype(ctx.cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.cache(Flavor::Trial{}))>
            >
          >;
        return AliasedMatrixData<Context,MatrixView,_mode>{std::move(ctx)};
      }
      else
      {
        using MatrixView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstUncachedMatrixView<
            const Matrix,
            std::decay_t<decltype(ctx.cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.cache(Flavor::Trial{}))>
            >,
          UncachedMatrixView<
            Matrix,
            std::decay_t<decltype(ctx.cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.cache(Flavor::Trial{}))>
            >
          >;
        return CachedMatrixData<Context,MatrixView,_mode>{std::move(ctx)};
      }
    }

    template<typename Implementation>
    struct CellJacobianData
      : public Implementation
    {

    protected:

      using JacobianImplementation = Implementation;

    public:

      using Context_ = typename Implementation::Context_;

      using Jacobian = typename Implementation::Traits;

      typename Jacobian::AccumulationView jacobian()
      {
        return Implementation::accumulationView(Context_::engine().weight());
      }

      template<typename TestSpace, typename TrialSpace>
      std::enable_if_t<
        std::is_same_v<typename TestSpace::Flavor,typename Context_::Flavor::Test> and
        (std::is_same_v<typename TrialSpace::Flavor,typename Context_::Flavor::Trial> or
         (Context_::isGalerkin() and std::is_same_v<typename TrialSpace::Flavor,typename Context_::Flavor::Test>)),
        LocalMatrixProxy<typename Jacobian::View,TestSpace,TrialSpace>
        >
      jacobian(const TestSpace& test_space, const TrialSpace& trial_space)
      {
        return {Implementation::readWriteView(),test_space,trial_space,Context_::engine().weight()};
      }

      void setup()
      {
        Implementation::setup(Context_::engine().jacobian());
      }

      using Context_::bind;
      using Context_::unbind;

      Context_* bind(const typename Context_::Entity&, typename Context_::Index, typename Context_::Index)
      {
        Implementation::bind(Context_::test().cache(),Context_::trial().cache());
        return this;
      }

      Context_* unbind(const typename Context_::Entity&, typename Context_::Index, typename Context_::Index)
      {
        Implementation::unbind(Context_::test().cache(),Context_::trial().cache());
        return this;
      }

      CellJacobianData(Implementation&& implementation)
        : Implementation(std::move(implementation))
      {
        static_assert(std::is_same_v<typename Context_::Test,typename Context_::Trial>,"foo");

      }

    };

    template<typename Implementation>
    auto cellJacobianData(Implementation&& implementation)
    {
      return CellJacobianData<Implementation>(std::move(implementation));
    }

    template<typename Context>
    struct CellTimeJacobianData
      : public Context
    {

      using Context_ = Context;

      using Jacobian = typename Context::Jacobian;

      typename Jacobian::AccumulationView timeJacobian()
      {
        return Context_::JacobianImplementation::accumulationView(Context_::engine().timeWeight());
      }

      template<typename TestSpace, typename TrialSpace>
      std::enable_if_t<
        std::is_same_v<typename TestSpace::Flavor,typename Context_::Flavor::Test> and
        (std::is_same_v<typename TrialSpace::Flavor,typename Context_::Flavor::Trial> or
         (Context_::isGalerkin() and std::is_same_v<typename TrialSpace::Flavor,typename Context_::Flavor::Test>)),
        LocalMatrixProxy<typename Jacobian::View,TestSpace,TrialSpace>
        >
      timeJacobian(const TestSpace& test_space, const TrialSpace& trial_space)
      {
        return {Jacobian::readWriteView(),test_space,trial_space,Context_::engine().timeWeight()};
      }

      CellTimeJacobianData(Context&& context)
        : Context(std::move(context))
      {}

    };

    template<typename Context>
    auto cellTimeJacobianData(std::true_type,Context&& context)
    {
      return CellTimeJacobianData<Context>(std::move(context));
    }

    template<typename Context>
    auto cellTimeJacobianData(std::false_type,Context&& context)
    {
      return std::move(context);
    }

    template<typename Context, typename LM_IO, typename LM_OI, LocalViewDataMode _mode>
    class CachedSkeletonMatrixData
      : public Context
    {

    protected:

      using Context_ = Context;

      using InsideOutsideTraits = CachedMatrixDataTraits<LM_IO>;
      using OutsideInsideTraits = CachedMatrixDataTraits<LM_OI>;

      void setup(typename InsideOutsideTraits::Matrix& matrix)
      {
        Context::setup();
        _local_view_io.attach(matrix);
        _local_view_oi.attach(matrix);
      }

      using Context::bind;
      using Context::unbind;

      Context* bind()
      {
        auto&& inside_test_cache = Context::inside().test().cache();
        auto&& inside_trial_cache = Context::inside().trial().cache();
        auto&& outside_test_cache = Context::outside().test().cache();
        auto&& outside_trial_cache = Context::outside().trial().cache();

        _local_view_io.bind(inside_test_cache,outside_trial_cache);
        _local_view_oi.bind(outside_test_cache,inside_trial_cache);

        if constexpr (_mode == LocalViewDataMode::read or _mode == LocalViewDataMode::readWrite or _mode == LocalViewDataMode::readAccumulate)
          {
            _container_io.resize(inside_test_cache.size(),outside_trial_cache.size());
            _container_oi.resize(outside_test_cache.size(),inside_trial_cache.size());

            _local_view_io.read(_container_io);
            _local_view_oi.read(_container_oi);
          }
        else if constexpr (_mode == LocalViewDataMode::write or _mode == LocalViewDataMode::accumulate)
          {
            _container_io.assign(inside_test_cache.size(),outside_trial_cache.size(),_initial);
            _container_oi.assign(outside_test_cache.size(),inside_trial_cache.size(),_initial);
          }
        return this;
      }

      Context* unbind()
      {
        if constexpr (_mode == LocalViewDataMode::write or _mode == LocalViewDataMode::readWrite)
          {
            Context::engine().scatterJacobian(_container_io,_local_view_io);
            Context::engine().scatterJacobian(_container_oi,_local_view_oi);

            _local_view_io.commit();
            _local_view_oi.commit();
          }
        else if constexpr (_mode == LocalViewDataMode::accumulate or _mode == LocalViewDataMode::readAccumulate)
          {
            Context::engine().scatterJacobian(_container_io,_local_view_io);
            Context::engine().scatterJacobian(_container_oi,_local_view_oi);

            _local_view_io.commit();
            _local_view_oi.commit();
          }

        _local_view_io.unbind();
        _local_view_oi.unbind();
        return this;
      }

      typename InsideOutsideTraits::AccumulationView insideOutsideAccumulationView(typename InsideOutsideTraits::Weight weight)
      {
        return _container_io.weightedAccumulationView(weight);
      }

      typename OutsideInsideTraits::AccumulationView outsideInsideAccumulationView(typename InsideOutsideTraits::Weight weight)
      {
        return _container_oi.weightedAccumulationView(weight);
      }

      const typename InsideOutsideTraits::View& insideOutsideReadOnlyView()
      {
        return _container_io;
      }

      typename InsideOutsideTraits::View& insideOutsideReadWriteView()
      {
        return _container_io;
      }

      const typename OutsideInsideTraits::View& outsideInsideReadOnlyView()
      {
        return _container_oi;
      }

      typename OutsideInsideTraits::View& outsideInsideReadWriteView()
      {
        return _container_oi;
      }

    public:

      CachedSkeletonMatrixData(Context&& ctx, typename InsideOutsideTraits::value_type initial = typename InsideOutsideTraits::value_type(0))
        : Context(std::move(ctx))
        , _initial(initial)
      {}

    private:

      typename InsideOutsideTraits::LocalView _local_view_io;
      typename OutsideInsideTraits::LocalView _local_view_oi;
      typename InsideOutsideTraits::Container _container_io;
      typename OutsideInsideTraits::Container _container_oi;
      typename InsideOutsideTraits::value_type _initial;

    };


    template<typename Context, typename LM_IO, typename LM_OI, LocalViewDataMode _mode>
    class AliasedSkeletonMatrixData
      : public Context
    {

    protected:

      using Context_ = Context;

      using InsideOutsideTraits = AliasedMatrixDataTraits<LM_IO>;
      using OutsideInsideTraits = AliasedMatrixDataTraits<LM_OI>;

      void setup(typename InsideOutsideTraits::Matrix& matrix)
      {
        Context::setup();
        _local_view_io.attach(matrix);
        _local_view_oi.attach(matrix);
      }

      using Context::bind;
      using Context::unbind;

      Context* bind()
      {
        auto&& inside_test_cache = Context::inside().test().cache();
        auto&& inside_trial_cache = Context::inside().trial().cache();
        auto&& outside_test_cache = Context::outside().test().cache();
        auto&& outside_trial_cache = Context::outside().trial().cache();

        _local_view_io.bind(inside_test_cache,outside_trial_cache);
        _local_view_oi.bind(outside_test_cache,inside_trial_cache);

        return this;
      }

      Context* unbind()
      {
        if constexpr (
          _mode == LocalViewDataMode::write or
          _mode == LocalViewDataMode::readWrite or
          _mode == LocalViewDataMode::accumulate or
          _mode == LocalViewDataMode::readAccumulate
          )
        {
          _local_view_io.commit();
          _local_view_oi.commit();
        }

        _local_view_io.unbind();
        _local_view_oi.unbind();

        return this;
      }

      typename InsideOutsideTraits::AccumulationView insideOutsideAccumulationView(typename InsideOutsideTraits::Weight weight)
      {
        DUNE_THROW(NotImplemented,"accumulationView() not implemented");
      }

      typename OutsideInsideTraits::AccumulationView outsideInsideAccumulationView(typename InsideOutsideTraits::Weight weight)
      {
        DUNE_THROW(NotImplemented,"accumulationView() not implemented");
      }

      const typename InsideOutsideTraits::LocalView& insideOutsideReadOnlyView()
      {
        return _local_view_io;
      }

      typename InsideOutsideTraits::LocalView& insideOutsideReadWriteView()
      {
        return _local_view_io;
      }

      const typename OutsideInsideTraits::LocalView& outsideInsideReadOnlyView()
      {
        return _local_view_oi;
      }

      typename OutsideInsideTraits::LocalView& outsideInsideReadWriteView()
      {
        return _local_view_oi;
      }

    public:

      AliasedSkeletonMatrixData(Context&& ctx, typename InsideOutsideTraits::value_type initial = typename InsideOutsideTraits::value_type(0))
        : Context(std::move(ctx))
      {}

    private:

      typename InsideOutsideTraits::LocalView _local_view_io;
      typename OutsideInsideTraits::LocalView _local_view_oi;

    };


    template<
      typename Matrix,
      LocalViewDataMode _mode,
      typename Field,
      typename Context>
    auto skeletonMatrixData(const Field& initial, Context&& ctx)
    {
      if constexpr (Context::fastDG())
      {
        using InsideOutsideMatrixView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstAliasedMatrixView<
            const Matrix,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Trial{}))>
            >,
          AliasedMatrixView<
            Matrix,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Trial{}))>
            >
          >;
        using OutsideInsideMatrixView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstAliasedMatrixView<
            const Matrix,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Trial{}))>
            >,
          AliasedMatrixView<
            Matrix,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Trial{}))>
            >
          >;
        return AliasedSkeletonMatrixData<
          Context,
          InsideOutsideMatrixView,
          OutsideInsideMatrixView,
          _mode
          >{std::move(ctx),initial};
      }
      else
      {
        using InsideOutsideMatrixView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstUncachedMatrixView<
            const Matrix,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Trial{}))>
            >,
          UncachedMatrixView<
            Matrix,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Trial{}))>
            >
          >;
        using OutsideInsideMatrixView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstUncachedMatrixView<
            const Matrix,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Trial{}))>
            >,
          UncachedMatrixView<
            Matrix,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Trial{}))>
            >
          >;
        return CachedSkeletonMatrixData<
          Context,
          InsideOutsideMatrixView,
          OutsideInsideMatrixView,
          _mode
          >{std::move(ctx),initial};
      }
    }


    template<
      typename Matrix,
      LocalViewDataMode _mode,
      typename Context>
    auto skeletonMatrixData(Context&& ctx)
    {
      if constexpr (Context::fastDG())
      {
        using InsideOutsideMatrixView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstAliasedMatrixView<
            const Matrix,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Trial{}))>
            >,
          AliasedMatrixView<
            Matrix,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Trial{}))>
            >
          >;
        using OutsideInsideMatrixView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstAliasedMatrixView<
            const Matrix,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Trial{}))>
            >,
          AliasedMatrixView<
            Matrix,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Trial{}))>
            >
          >;
        return AliasedSkeletonMatrixData<
          Context,
          InsideOutsideMatrixView,
          OutsideInsideMatrixView,
          _mode
          >{std::move(ctx)};
      }
      else
      {
        using InsideOutsideMatrixView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstUncachedMatrixView<
            const Matrix,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Trial{}))>
            >,
          UncachedMatrixView<
            Matrix,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Trial{}))>
            >
          >;
        using OutsideInsideMatrixView = std::conditional_t<
          _mode == LocalViewDataMode::read,
          ConstUncachedMatrixView<
            const Matrix,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Trial{}))>
            >,
          UncachedMatrixView<
            Matrix,
            std::decay_t<decltype(ctx.outside().cache(Flavor::Test{}))>,
            std::decay_t<decltype(ctx.inside().cache(Flavor::Trial{}))>
            >
          >;
        return CachedSkeletonMatrixData<
          Context,
          InsideOutsideMatrixView,
          OutsideInsideMatrixView,
          _mode
          >{std::move(ctx)};
      }
    }


    template<typename Implementation>
    struct SkeletonJacobianData
      : public Implementation
    {

      using Context_ = typename Implementation::Context_;

      using InsideOutsideJacobian = typename Implementation::InsideOutsideTraits;
      using OutsideInsideJacobian = typename Implementation::OutsideInsideTraits;

      SkeletonJacobianData(Implementation&& implementation)
        : Implementation(std::move(implementation))
      {}


      template<typename Inside>
      std::enable_if_t<
        std::is_same_v<Inside,typename Implementation::Inside>,
        typename Context_::Inside::Jacobian::AccumulationView
        >
      jacobian(const Inside&, const Inside&)
      {
        return Context_::inside().accumulationView(Context_::engine().weight());
      }

      template<typename Outside>
      std::enable_if_t<
        std::is_same_v<Outside,typename Implementation::Outside>,
        typename Context_::Outside::Jacobian::AccumulationView
        >
      jacobian(const Outside&, const Outside&)
      {
        return Context_::outside().accumulationView(Context_::engine().weight());
      }

      template<typename Inside, typename Outside>
      std::enable_if_t<
        std::is_same_v<Inside,typename Implementation::Inside> and std::is_same_v<Outside,typename Implementation::Outside>,
        typename InsideOutsideJacobian::AccumulationView
        >
      jacobian(const Inside&, const Outside&)
      {
        return Implementation::insideOutsideAccumulationView(Context_::engine().weight());
      }

      template<typename Outside, typename Inside>
      std::enable_if_t<
        std::is_same_v<Outside,typename Implementation::Outside> and std::is_same_v<Inside,typename Implementation::Inside>,
        typename OutsideInsideJacobian::AccumulationView
        >
      jacobian(const Outside&, const Inside&)
      {
        return Implementation::outsideInsideAccumulationView(Context_::engine().weight());
      }

      template<typename TestSpace, typename TrialSpace>
      std::enable_if_t<
        std::is_same_v<typename TestSpace::Flavor,Flavor::InsideTest> and
        (std::is_same_v<typename TrialSpace::Flavor,Flavor::InsideTrial> or
         (Context_::isGalerkin() and std::is_same_v<typename TrialSpace::Flavor,Flavor::InsideTest>)),
        LocalMatrixProxy<typename Context_::Inside::Jacobian::View,TestSpace,TrialSpace>
        >
      jacobian(const TestSpace& test_space, const TrialSpace& trial_space)
      {
        return Context_::inside().jacobian(test_space,trial_space);
      }

      template<typename TestSpace, typename TrialSpace>
      std::enable_if_t<
        std::is_same_v<typename TestSpace::Flavor,Flavor::OutsideTest> and
        (std::is_same_v<typename TrialSpace::Flavor,Flavor::OutsideTrial> or
         (Context_::isGalerkin() and std::is_same_v<typename TrialSpace::Flavor,Flavor::OutsideTest>)),
        LocalMatrixProxy<typename Context_::Outside::Jacobian::View,TestSpace,TrialSpace>
        >
      jacobian(const TestSpace& test_space, const TrialSpace& trial_space)
      {
        return Context_::outside().jacobian(test_space,trial_space);
      }

      template<typename TestSpace, typename TrialSpace>
      std::enable_if_t<
        std::is_same_v<typename TestSpace::Flavor,Flavor::InsideTest> and
        (std::is_same_v<typename TrialSpace::Flavor,Flavor::OutsideTrial> or
         (Context_::isGalerkin() and std::is_same_v<typename TrialSpace::Flavor,Flavor::OutsideTest>)),
        LocalMatrixProxy<typename InsideOutsideJacobian::View,TestSpace,TrialSpace>
        >
      jacobian(const TestSpace& test_space, const TrialSpace& trial_space)
      {
        return {Implementation::insideOutsideReadWriteView(),test_space,trial_space,Context_::engine().weight()};
      }

      template<typename TestSpace, typename TrialSpace>
      std::enable_if_t<
        std::is_same_v<typename TestSpace::Flavor,Flavor::OutsideTest> and
        (std::is_same_v<typename TrialSpace::Flavor,Flavor::InsideTrial> or
         (Context_::isGalerkin() and std::is_same_v<typename TrialSpace::Flavor,Flavor::InsideTest>)),
        LocalMatrixProxy<typename OutsideInsideJacobian::View,TestSpace,TrialSpace>
        >
      jacobian(const TestSpace& test_space, const TrialSpace& trial_space)
      {
        return {Implementation::outsideInsideReadWriteView(),test_space,trial_space,Context_::engine().weight()};
      }

      void setup()
      {
        Implementation::setup(Context_::engine().jacobian());
      }

      using Context_::bind;
      using Context_::unbind;


      template<typename IntersectionType_>
      Context_* bind(
        IntersectionType_ type,
        const typename Context_::IntersectionDomain::Intersection&,
        typename Context_::IntersectionDomain::Index,
        const typename Context_::IntersectionDomain::Entity& entity,
        typename Context_::IntersectionDomain::Index entity_index,
        typename Context_::IntersectionDomain::Index unique_index
        )
      {
        if (type == IntersectionType::skeleton or type == IntersectionType::periodic)
          Implementation::bind();
        return this;
      }

      template<typename IntersectionType_>
      Context_* unbind(
        IntersectionType_ type,
        const typename Context_::IntersectionDomain::Intersection&,
        typename Context_::IntersectionDomain::Index,
        const typename Context_::IntersectionDomain::Entity& entity,
        typename Context_::IntersectionDomain::Index entity_index,
        typename Context_::IntersectionDomain::Index unique_index
        )
      {
        if (type == IntersectionType::skeleton or type == IntersectionType::periodic)
          Implementation::unbind();
        return this;
      }

    };

    template<typename Implementation>
    auto skeletonJacobianData(Implementation&& implementation)
    {
      return SkeletonJacobianData<Implementation>(std::move(implementation));
    }


    template<typename Context>
    struct SkeletonTimeJacobianData
      : public Context
    {

      using Context_ = Context;

      using InsideOutsideJacobian = typename Context_::InsideOutsideTraits;
      using OutsideInsideJacobian = typename Context_::OutsideInsideTraits;

      SkeletonTimeJacobianData(Context&& context)
        : Context(std::move(context))
      {}


      template<typename Inside>
      std::enable_if_t<
        std::is_same_v<Inside,typename Context_::Inside>,
        typename Context_::Inside::Jacobian::AccumulationView
        >
      timeJacobian(const Inside&, const Inside&)
      {
        return Context_::inside().accumulationView(Context_::engine().timeWeight());
      }

      template<typename Outside>
      std::enable_if_t<
        std::is_same_v<Outside,typename Context_::Outside>,
        typename Context_::Outside::Jacobian::AccumulationView
        >
      timeJacobian(const Outside&, const Outside&)
      {
        return Context_::outside().accumulationView(Context_::engine().timeWeight());
      }

      template<typename Inside, typename Outside>
      std::enable_if_t<
        std::is_same_v<Inside,typename Context_::Inside> and std::is_same_v<Outside,typename Context_::Outside>,
        typename InsideOutsideJacobian::AccumulationView
        >
      timeJacobian(const Inside&, const Outside&)
      {
        return Context_::insideOutsideAccumulationView(Context_::engine().timeWeight());
      }

      template<typename Outside, typename Inside>
      std::enable_if_t<
        std::is_same_v<Outside,typename Context_::Outside> and std::is_same_v<Inside,typename Context_::Inside>,
        typename OutsideInsideJacobian::AccumulationView
        >
      timeJacobian(const Outside&, const Inside&)
      {
        return Context_::outsideInsideAccumulationView(Context_::engine().timeWeight());
      }

      template<typename TestSpace, typename TrialSpace>
      std::enable_if_t<
        std::is_same_v<typename TestSpace::Flavor,Flavor::InsideTest> and
        (std::is_same_v<typename TrialSpace::Flavor,Flavor::InsideTrial> or
         (Context_::isGalerkin() and std::is_same_v<typename TrialSpace::Flavor,Flavor::InsideTest>)),
        LocalMatrixProxy<typename Context_::Inside::Jacobian::View,TestSpace,TrialSpace>
        >
      timeJacobian(const TestSpace& test_space, const TrialSpace& trial_space)
      {
        return Context_::inside().timeJacobian(test_space,trial_space);
      }

      template<typename TestSpace, typename TrialSpace>
      std::enable_if_t<
        std::is_same_v<typename TestSpace::Flavor,Flavor::OutsideTest> and
        (std::is_same_v<typename TrialSpace::Flavor,Flavor::OutsideTrial> or
         (Context_::isGalerkin() and std::is_same_v<typename TrialSpace::Flavor,Flavor::OutsideTest>)),
        LocalMatrixProxy<typename Context_::Outside::Jacobian::View,TestSpace,TrialSpace>
        >
      timeJacobian(const TestSpace& test_space, const TrialSpace& trial_space)
      {
        return Context_::outside().timeJacobian(test_space,trial_space);
      }

      template<typename TestSpace, typename TrialSpace>
      std::enable_if_t<
        std::is_same_v<typename TestSpace::Flavor,Flavor::InsideTest> and
        (std::is_same_v<typename TrialSpace::Flavor,Flavor::OutsideTrial> or
         (Context_::isGalerkin() and std::is_same_v<typename TrialSpace::Flavor,Flavor::OutsideTest>)),
        LocalMatrixProxy<typename InsideOutsideJacobian::View,TestSpace,TrialSpace>
        >
      timeJacobian(const TestSpace& test_space, const TrialSpace& trial_space)
      {
        return {Context_::insideOutsideReadWriteView(),test_space,trial_space,Context_::engine().timeWeight()};
      }

      template<typename TestSpace, typename TrialSpace>
      std::enable_if_t<
        std::is_same_v<typename TestSpace::Flavor,Flavor::OutsideTest> and
        (std::is_same_v<typename TrialSpace::Flavor,Flavor::InsideTrial> or
         (Context_::isGalerkin() and std::is_same_v<typename TrialSpace::Flavor,Flavor::InsideTest>)),
        LocalMatrixProxy<typename OutsideInsideJacobian::View,TestSpace,TrialSpace>
        >
      timeJacobian(const TestSpace& test_space, const TrialSpace& trial_space)
      {
        return {Context_::outsideInsideReadWriteView(),test_space,trial_space,Context_::engine().timeWeight()};
      }

    };

    template<typename Context>
    auto skeletonTimeJacobianData(std::true_type, Context&& context)
    {
      return SkeletonTimeJacobianData<Context>(std::move(context));
    }

    template<typename Context>
    auto skeletonTimeJacobianData(std::false_type, Context&& context)
    {
      return std::move(context);
    }


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_MATRIXDATA_HH
