// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_FUNCTIONSPACEDATA_HH
#define DUNE_PDELAB_ASSEMBLER_FUNCTIONSPACEDATA_HH

#include <cassert>
#include <type_traits>

#include <dune/typetree/childextraction.hh>
#include <dune/typetree/transformation.hh>
#include <dune/typetree/transformationutilities.hh>

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/assembler/utility.hh>
#include <dune/pdelab/assembler/celldata.hh>
#include <dune/pdelab/assembler/finiteelementwrapper.hh>
#include <dune/pdelab/assembler/quadraturerule.hh>

namespace Dune {
  namespace PDELab {

    template<typename Context>
    struct lfs_to_finite_elements
    {};

    template<typename PowerLocalFunctionSpace, typename Context>
    TypeTree::SimplePowerNodeTransformation<
      PowerLocalFunctionSpace,
      lfs_to_finite_elements<Context>,
      TypeTree::PowerNode
      >
    registerNodeTransformation(PowerLocalFunctionSpace* plfs, lfs_to_finite_elements<Context>* t, PowerLocalFunctionSpaceTag* tag);

    template<typename CompositeLocalFunctionSpace, typename Context>
    TypeTree::SimpleCompositeNodeTransformation<
      CompositeLocalFunctionSpace,
      lfs_to_finite_elements<Context>,
      TypeTree::CompositeNode
      >
    registerNodeTransformation(CompositeLocalFunctionSpace* plfs, lfs_to_finite_elements<Context>* t, CompositeLocalFunctionSpaceTag* tag);

    template<typename LocalFunctionSpace, typename Context>
    TypeTree::SimpleLeafNodeTransformation<
      LocalFunctionSpace,
      lfs_to_finite_elements<Context>,
      FiniteElementWrapper<typename LocalFunctionSpace::Traits::FiniteElement,Context>
      >
    registerNodeTransformation(LocalFunctionSpace* plfs, lfs_to_finite_elements<Context>* t, LeafLocalFunctionSpaceTag* tag);


    template<typename Context, typename LFS, typename LFSCache>
    struct FunctionSpaceData
      : public Context
    {

      using FunctionSpace      = LFS;
      using FunctionSpaceCache = LFSCache;
      using Cache              = LFSCache;

    private:

      using fe_transformation = TypeTree::TransformTree<LFS,lfs_to_finite_elements<Context>>;

    public:

      using FiniteElements    = typename fe_transformation::transformed_type;

      template<std::size_t... I>
      using SubSpace = TypeTree::Child<LFS,I...>;

      template<std::size_t... I>
      using FiniteElement = TypeTree::Child<FiniteElements,I...>;

      template<std::size_t... I>
      using Basis = typename FiniteElement<I...>::Basis;

      using Context::bind;

      Context* bind(const typename Context::Entity& element, typename Context::Index, typename Context::Index)
      {
        _function_space.bind(element);
        TypeTree::applyToTreePair(_function_space,_finite_elements,set_finite_elements{});
        _function_space_cache.update();
        return this;
      }

      template<typename... Indices>
      const auto& functionSpace(Indices... indices) const
      {
        return child(_function_space,indices...);
      }

      template<typename... Indices>
      const auto& space(Indices... indices) const
      {
        return child(_function_space,indices...);
      }

      template<typename... Indices>
      auto& finiteElement(Indices... indices)
      {
        return TypeTree::child(_finite_elements,indices...);
      }

      template<typename... Indices>
      auto& basis(Indices... indices)
      {
        return finiteElement(indices...).basis();
      }

      const FunctionSpaceCache& functionSpaceCache() const
      {
        return _function_space_cache;
      }

      using Context::cache;

      Cache& cache()
      {
        return _function_space_cache;
      }

      FunctionSpaceData(Context&& ctx, const typename LFS::Traits::GridFunctionSpace& gfs, LFSCache&& cache)
        : Context(std::move(ctx))
        , _function_space(gfs)
        , _function_space_cache(std::move(cache))
        , _finite_elements(fe_transformation::transform(_function_space))
      {}

      Context* setup()
      {
        _function_space_cache.attach(_function_space);
        TypeTree::applyToTree(_finite_elements,set_context<Context>(*this));
        return this;
      }

    private:

      FunctionSpace _function_space;
      FunctionSpaceCache _function_space_cache;
      FiniteElements _finite_elements;

    };


    template<typename Context>
    struct TestSpaceData
      : public FunctionSpaceData<Context,CellFlavor::TestLocalSpace<Context>,CellFlavor::TestSpaceCache<Context>>
    {

      // avoid introducing name Context_, otherwise mayhem may occur
      using Context_ = FunctionSpaceData<Context,CellFlavor::TestLocalSpace<Context>,CellFlavor::TestSpaceCache<Context>>;

      using Test = Context_;

      Test& test()
      {
        return *this;
      }

      using Context_::cache;

      typename Context_::Cache& cache(Flavor::Test)
      {
        return Context_::cache();
      }

      TestSpaceData(Context&& ctx)
        : Context_(std::move(ctx),ctx.engine().testSpace(),ctx.engine().makeTestSpaceCache(typename Context::Flavor::Test{}))
      {}

    };

    template<typename Context>
    auto testSpaceData(Context&& ctx)
    {
      return TestSpaceData<Context>(std::move(ctx));
    }

    template<typename Context>
    struct GalerkinTrialSpaceData
      : public Context
    {

      using Trial = typename Context::Test;

      Trial& trial()
      {
        return *this;
      }

      GalerkinTrialSpaceData(Context&& ctx)
        : Context(std::move(ctx))
      {}

      using Context::cache;

      typename Context::Cache& cache(Flavor::Trial)
      {
        return Context::cache();
      }

    };

    template<typename Context>
    auto galerkinTrialSpaceData(Context&& ctx)
    {
      return GalerkinTrialSpaceData<Context>(std::move(ctx));
    }

    template<typename Context>
    struct NonGalerkinTrialSpaceData
      : public FunctionSpaceData<Context,CellFlavor::TrialLocalSpace<Context>,CellFlavor::TrialSpaceCache<Context>>
    {

      using Context_ = FunctionSpaceData<Context,CellFlavor::TrialLocalSpace<Context>,CellFlavor::TrialSpaceCache<Context>>;

      using Trial = Context_;

      Trial& trial()
      {
        return *this;
      }

      using Context_::cache;

      typename Context_::Cache& cache(Flavor::Trial)
      {
        return Context_::cache();
      }

      NonGalerkinTrialSpaceData(Context&& ctx)
        : Context_(std::move(ctx),ctx.engine().trialSpace(),ctx.engine().makeTrialSpaceCache(typename Context::Flavor::Trial{}))
      {}

    };

    template<typename Context>
    auto nonGalerkinTrialSpaceData(Context&& ctx)
    {
      return NonGalerkinTrialSpaceData<Context>(std::move(ctx));
    }

    template<typename Context>
    auto trialSpaceData(Context&& ctx, std::enable_if_t<Context::isGalerkin(),int> = 0)
    {
      return GalerkinTrialSpaceData<Context>(std::move(ctx));
    }

    template<typename Context>
    auto trialSpaceData(Context&& ctx, std::enable_if_t<not Context::isGalerkin(),int> = 0)
    {
      return NonGalerkinTrialSpaceData<Context>(std::move(ctx));
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_FUNCTIONSPACEDATA_HH
