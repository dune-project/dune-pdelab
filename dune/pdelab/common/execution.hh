// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_EXECUTION_HH
#define DUNE_PDELAB_COMMON_EXECUTION_HH

/**
 * \file
 * \brief PDELab execution policies.
 */

#include <type_traits>

namespace Dune {
  namespace PDELab {
    namespace Execution {

      class SequencedPolicy{};
      class ParallelPolicy{};
      class ParallelUnsequencedPolicy{};
      class UnsequencedPolicy{};

      inline constexpr SequencedPolicy seq{};
      inline constexpr ParallelPolicy par{};
      inline constexpr ParallelUnsequencedPolicy par_unseq{};
      inline constexpr UnsequencedPolicy unseq{};

      template <class T>
      struct is_execution_policy : ::std::false_type {};
      template <>
      struct is_execution_policy<Dune::PDELab::Execution::SequencedPolicy> : ::std::true_type {};
      template <>
      struct is_execution_policy<Dune::PDELab::Execution::ParallelPolicy> : ::std::true_type {};
      template <>
      struct is_execution_policy<Dune::PDELab::Execution::ParallelUnsequencedPolicy> : ::std::true_type {};
      template <>
      struct is_execution_policy<Dune::PDELab::Execution::UnsequencedPolicy> : ::std::true_type {};

      template <class T>
      inline constexpr bool is_execution_policy_v = Dune::PDELab::Execution::is_execution_policy<T>::value;

    } // namespace Execution
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_EXECUTION_HH
