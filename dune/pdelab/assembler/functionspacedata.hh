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

    template<typename Cell>
    struct lfs_to_finite_elements
    {};

    template<typename PowerLocalFunctionSpace, typename Cell>
    TypeTree::SimplePowerNodeTransformation<
      PowerLocalFunctionSpace,
      lfs_to_finite_elements<Cell>,
      TypeTree::PowerNode
      >
    registerNodeTransformation(PowerLocalFunctionSpace* plfs, lfs_to_finite_elements<Cell>* t, PowerLocalFunctionSpaceTag* tag);

    template<typename CompositeLocalFunctionSpace, typename Cell>
    TypeTree::SimpleCompositeNodeTransformation<
      CompositeLocalFunctionSpace,
      lfs_to_finite_elements<Cell>,
      TypeTree::CompositeNode
      >
    registerNodeTransformation(CompositeLocalFunctionSpace* plfs, lfs_to_finite_elements<Cell>* t, CompositeLocalFunctionSpaceTag* tag);

    template<typename LocalFunctionSpace, typename Cell>
    TypeTree::SimpleLeafNodeTransformation<
      LocalFunctionSpace,
      lfs_to_finite_elements<Cell>,
      FiniteElementWrapper<typename LocalFunctionSpace::Traits::FiniteElement,Cell>
      >
    registerNodeTransformation(LocalFunctionSpace* plfs, lfs_to_finite_elements<Cell>* t, LeafLocalFunctionSpaceTag* tag);


    template<typename LFS, typename LFSCache, typename Cell>
    struct FunctionSpaceData
    {

      using FunctionSpace      = LFS;
      using FunctionSpaceCache = LFSCache;
      using Cache              = LFSCache;
      using CellType           = Cell;

    private:

      using FiniteElements = typename TypeTree::TransformTree<LFS,lfs_to_finite_elements<Cell>>::transformed_type;

    public:

      template<std::size_t... I>
      using SubSpace = TypeTree::Child<LFS,I...>;

      template<std::size_t... I>
      using FiniteElement = TypeTree::Child<FiniteElements,I...>;

      template<std::size_t... I>
      using Basis = typename FiniteElement<I...>::Basis;

      template<typename Context, typename CellContext, typename Element, typename Index>
      void bind(Context&, CellContext&, const Element& element, Index, Index)
      {
        _function_space.bind(element);
        TypeTree::applyToTreePair(_function_space,_finite_elements,set_finite_elements{});
        _function_space_cache.update();
      }

      template<typename... Indices>
      const auto& functionSpace(Indices... indices) const
      {
        return child(_function_space,indices...);
      }

      template<typename... Indices>
      auto& finiteElement(Indices... indices)
      {
        //using FiniteElement = std::decay_t<decltype(child(_function_space,indices...).finiteElement())>;
        //return FiniteElementWrapper<FiniteElement>(child(_function_space,indices...).finiteElement());
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

      const Cache& cache() const
      {
        return _function_space_cache;
      }

      template<typename GFS, typename... Args>
      FunctionSpaceData(const GFS& gfs, Args&&... args)
        : _function_space(gfs)
        , _function_space_cache(_function_space,std::forward<Args>(args)...)
        , _finite_elements(TypeTree::TransformTree<LFS,lfs_to_finite_elements<Cell>>::transform(_function_space))
      {}

      template<typename Context, typename CellContext, typename Engine>
      void setup(Context&, CellContext&, Engine&)
      {
        _function_space_cache.setLocalFunctionSpace(_function_space);
      }

    private:

      FunctionSpace _function_space;
      FunctionSpaceCache _function_space_cache;
      FiniteElements _finite_elements;

    };


    template<typename LFS, typename LFSCache, typename Cell>
    struct TestSpaceData
      : public FunctionSpaceData<LFS,LFSCache,Cell>
    {

      using Test = FunctionSpaceData<LFS,LFSCache,Cell>;

      Test& test()
      {
        return *this;
      }

      template<typename GFS, typename... Args>
      TestSpaceData(const GFS& gfs, Args&&... args)
        : Test(gfs,std::forward<Args>(args)...)
      {}

    };

    template<typename LFS, typename LFSCache, typename Cell, typename... Args>
    auto testSpaceData(const typename LFS::Traits::GridFunctionSpace& gfs, Args&&... args)
    {
      return TestSpaceData<LFS,LFSCache,Cell>(gfs,std::forward<Args>(args)...);
    }


    template<typename TestSpaceData>
    struct GalerkinTrialSpaceData
      : public TestSpaceData
    {

      using Trial = typename TestSpaceData::Test;

      Trial& trial()
      {
        return *this;
      }

      GalerkinTrialSpaceData(TestSpaceData&& test_space_data)
        : TestSpaceData(std::move(test_space_data))
      {}

    };

    template<typename TestSpaceData>
    auto galerkinTrialSpaceData(TestSpaceData&& testSpaceData)
    {
      return GalerkinTrialSpaceData<TestSpaceData>(std::forward<TestSpaceData>(testSpaceData));
    }

    template<typename TestSpaceData, typename LFS, typename LFSCache>
    struct NonGalerkinTrialSpaceData
      : public TestSpaceData
    {

      using Trial = FunctionSpaceData<LFS,LFSCache,typename TestSpaceData::CellType>;

      Trial& trial()
      {
        return _trial;
      }

      template<typename Context, typename CellContext, typename Element, typename Index>
      void bind(Context& ctx, CellContext& cell_ctx, const Element& element, Index entity_index, Index unique_index)
      {
        TestSpaceData::bind(ctx,cell_ctx,element,entity_index,unique_index);
        _trial.bind(ctx,cell_ctx,element,entity_index,unique_index);
      }

      template<typename Context, typename CellContext, typename Engine>
      void setup(Context& ctx, CellContext& cell_ctx, Engine& engine)
      {
        TestSpaceData::setup(ctx,cell_ctx,engine);
        _trial.setup(ctx,cell_ctx,engine);
      }

      template<typename GFS, typename... Args>
      NonGalerkinTrialSpaceData(TestSpaceData&& test_space_data, const GFS& gfs, Args&&... args)
        : TestSpaceData(std::move(test_space_data))
        , _trial(gfs,std::forward<Args>(args)...)
      {}

    private:

      Trial _trial;

    };

    template<typename LFS, typename LFSCache, typename TestSpaceData, typename... Args>
    auto nonGalerkinTrialSpaceData(TestSpaceData&& testSpaceData, const typename LFS::Traits::GridFunctionSpace& gfs, Args&&... args)
    {
      return NonGalerkinTrialSpaceData<TestSpaceData,LFS,LFSCache>(std::forward<TestSpaceData>(testSpaceData),gfs,std::forward<Args>(args)...);
    }

    template<typename LFS, typename LFSCache, typename TestSpaceData, typename... Args>
    auto trialSpaceData(std::true_type, TestSpaceData&& testSpaceData, const typename LFS::Traits::GridFunctionSpace& gfs, Args&&... args)
    {
      return GalerkinTrialSpaceData<TestSpaceData>(std::forward<TestSpaceData>(testSpaceData));
    }

    template<typename LFS, typename LFSCache, typename TestSpaceData, typename... Args>
    auto trialSpaceData(std::false_type, TestSpaceData&& testSpaceData, const typename LFS::Traits::GridFunctionSpace& gfs, Args&&... args)
    {
      return NonGalerkinTrialSpaceData<TestSpaceData,LFS,LFSCache>(std::forward<TestSpaceData>(testSpaceData),gfs,std::forward<Args>(args)...);
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_FUNCTIONSPACEDATA_HH
