// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_ENGINEBASE_HH
#define DUNE_PDELAB_ASSEMBLER_ENGINEBASE_HH

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/flavor.hh>

namespace Dune {
  namespace PDELab {

    enum class Galerkin { disable, enable, automatic };

    static constexpr auto enableGalerkin = std::integral_constant<Galerkin,Galerkin::enable>{};
    static constexpr auto disableGalerkin = std::integral_constant<Galerkin,Galerkin::disable>{};
    static constexpr auto automaticGalerkin = std::integral_constant<Galerkin,Galerkin::automatic>{};

    template<typename TestSpace_,
             typename TrialSpace_,
             typename TrialConstraints_,
             typename TestConstraints_,
             bool enable_flavors>
    struct LocalFunctionSpaceTypes;

    template<
      typename TestSpace_,
      typename TrialSpace_,
      typename TrialConstraints_,
      typename TestConstraints_
      >
    struct LocalFunctionSpaceTypes<TestSpace_,TrialSpace_,TrialConstraints_,TestConstraints_,true>
    {

      using TestSpace        = TestSpace_;
      using TestConstraints  = TestConstraints_;

      template<typename Flavor = Flavor::Generic>
      using TestLocalSpace   = LocalFunctionSpace<TestSpace,Flavor>;

      template<typename Flavor = Flavor::Generic>
      using TestSpaceCache   = LFSIndexCache<TestLocalSpace<Flavor>,TestConstraints>;


      using TrialSpace       = TrialSpace_;
      using TrialConstraints = TrialConstraints_;

      template<typename Flavor = Flavor::Generic>
      using TrialLocalSpace = LocalFunctionSpace<TrialSpace,Flavor>;

      template<typename Flavor = Flavor::Generic>
      using TrialSpaceCache = LFSIndexCache<TrialLocalSpace<Flavor>,TrialConstraints>;

    };

    template<
      typename TestSpace_,
      typename TrialSpace_,
      typename TrialConstraints_,
      typename TestConstraints_
      >
    struct LocalFunctionSpaceTypes<TestSpace_,TrialSpace_,TrialConstraints_,TestConstraints_,false>
    {

      using TestSpace        = TestSpace_;
      using TestConstraints  = TestConstraints_;

      template<typename = Flavor::Generic>
      using TestLocalSpace   = LocalFunctionSpace<TestSpace,Flavor::Generic>;

      template<typename = Flavor::Generic>
      using TestSpaceCache   = LFSIndexCache<TestLocalSpace<Flavor::Generic>,TestConstraints>;


      using TrialSpace       = TrialSpace_;
      using TrialConstraints = TrialConstraints_;

      template<typename = Flavor::Generic>
      using TrialLocalSpace  = LocalFunctionSpace<TrialSpace,Flavor::Generic>;

      template<typename = Flavor::Generic>
      using TrialSpaceCache  = LFSIndexCache<TrialLocalSpace<Flavor::Generic>,TrialConstraints>;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_ENGINEBASE_HH
