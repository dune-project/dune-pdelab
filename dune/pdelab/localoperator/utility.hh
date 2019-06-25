// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_UTILITY_HH
#define DUNE_PDELAB_LOCALOPERATOR_UTILITY_HH

#include <cassert>
#include <type_traits>

#include <dune/typetree/childextraction.hh>

#include <dune/pdelab/assembler/context.hh>
#include <dune/pdelab/assembler/utility.hh>
#include <dune/pdelab/assembler/finiteelementwrapper.hh>

namespace Dune::PDELab {

  namespace LocalOperator {

    namespace Impl {

      template<typename LFS, typename = void>
      struct extract_finite_element_traits
      {
        using type = typename LFS::Traits::FiniteElement;
      };

      template<typename LFS>
      struct extract_finite_element_traits<LFS,std::void_t<TypeTree::Child<LFS,0>>>
        : public extract_finite_element_traits<TypeTree::Child<LFS,0>>
      {};

      template<typename T>
      struct extract_local_test_space
      {
        using type = typename T::Test::FunctionSpace;
      };

      template<typename Ctx>
      struct extract_local_test_space<Experimental::Context::Context<Ctx>>
      {
        using type = typename Ctx::Test::FunctionSpace;
      };

      template<typename Ctx>
      struct extract_local_test_space<Experimental::Context::CellContext<Ctx>>
      {
        using type = typename Ctx::Test::FunctionSpace;
      };

      template<typename Ctx>
      struct extract_local_test_space<Experimental::Context::IntersectionContext<Ctx>>
      {
        using type = typename Ctx::Test::FunctionSpace;
      };

    }

    template<typename T>
    using TestSpace = typename Impl::extract_local_test_space<T>::type;

    template<typename T>
    constexpr int dimension = TestSpace<T>::Traits::GridFunctionSpace::Traits::EntitySet::dimension;

    // TODO: Make these work with spaces as well

    template<typename T, std::size_t... I>
    using Range = typename Experimental::FiniteElementWrapper<typename Impl::extract_finite_element_traits<TypeTree::Child<TestSpace<T>,I...>>::type,T>::Basis::Range;

    template<typename T, std::size_t... I>
    using RangeField = typename Experimental::FiniteElementWrapper<typename Impl::extract_finite_element_traits<TypeTree::Child<TestSpace<T>,I...>>::type,T>::Basis::RangeField;

    template<typename T, std::size_t... I>
    using Gradient = typename Experimental::FiniteElementWrapper<typename Impl::extract_finite_element_traits<TypeTree::Child<TestSpace<T>,I...>>::type,T>::Basis::Gradient;

    template<typename T, std::size_t... I>
    using Values = typename Experimental::FiniteElementWrapper<typename Impl::extract_finite_element_traits<TypeTree::Child<TestSpace<T>,I...>>::type,T>::Basis::Values;

    template<typename T, std::size_t... I>
    using Gradients = typename Experimental::FiniteElementWrapper<typename Impl::extract_finite_element_traits<TypeTree::Child<TestSpace<T>,I...>>::type,T>::Basis::Gradients;

    template<typename T, std::size_t... I>
    using ReferenceGradients = typename Experimental::FiniteElementWrapper<typename Impl::extract_finite_element_traits<TypeTree::Child<TestSpace<T>,I...>>::type,T>::Basis::ReferenceGradients;

  } // namespace LocalOperator

} // namespace Dune::PDELab


#endif // DUNE_PDELAB_LOCALOPERATOR_UTILITY_HH
