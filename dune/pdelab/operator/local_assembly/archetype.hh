#ifndef DUNE_PDELAB_OPERATOR_LOCAL_ASSEMBLY_ARCHETYPE_HH
#define DUNE_PDELAB_OPERATOR_LOCAL_ASSEMBLY_ARCHETYPE_HH

#include <dune/pdelab/concepts/container.hh>
#include <dune/pdelab/concepts/local_basis.hh>

#include <dune/grid/concepts.hh>

#include <type_traits>

namespace Dune::PDELab::inline Experimental::LocalAssembly {

struct Archetype {

  constexpr static auto localAssembleIsLinear() noexcept {
    return std::false_type{};
  }

  constexpr static auto localAssembleSkipEntity(const Dune::Concept::Entity auto& entity) noexcept {
    return std::false_type{};
  }

  constexpr static auto localAssembleSkipIntersection(const Dune::Concept::Intersection auto& intersection) noexcept {
    return std::false_type{};
  }

  constexpr static auto localAssembleDoVolume() noexcept {
    return std::false_type{};
  }

  void localAssembleVolume(
                                                       auto  time,
    const PDELab::Concept::LocalBasis                  auto& ltrial,
    const PDELab::Concept::LocalConstContainer         auto& lcoefficients,
    const PDELab::Concept::LocalBasis                  auto& ltest,
    PDELab::Concept::LocalMutableContainer             auto& lresidual) noexcept;

  void localAssembleJacobianVolume(
                                                       auto  time,
    const PDELab::Concept::LocalBasis                  auto& ltrial,
    const PDELab::Concept::LocalConstContainer         auto& llin_point,
    const PDELab::Concept::LocalBasis                  auto& ltest,
    PDELab::Concept::LocalMutableMatrix                auto& ljacobian) noexcept;

  void localAssembleJacobianVolumeApply(
                                                       auto  time,
    const PDELab::Concept::LocalBasis                  auto& ltrial,
    const PDELab::Concept::LocalConstContainer         auto& llin_point,
    const PDELab::Concept::LocalConstContainer         auto& lapp_point,
    const PDELab::Concept::LocalBasis                  auto& ltest,
    PDELab::Concept::LocalMutableContainer             auto& ljacobian) noexcept;

  void localAssemblePatternVolume(
    const PDELab::Concept::LocalBasis                  auto& ltrial,
    const PDELab::Concept::LocalBasis                  auto& ltest,
                                                       auto& lpattern) noexcept;

  constexpr static auto localAssembleDoSkeleton() {
    return std::false_type{};
  }

  void localAssembleSkeleton(
    const Dune::Concept::Intersection                  auto& intersection,
                                                       auto  time,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer         auto& lcoefficients_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    const PDELab::Concept::LocalBasis                  auto& ltrial_out,
    const PDELab::Concept::LocalConstContainer         auto& lcoefficients_out,
    const PDELab::Concept::LocalBasis                  auto& ltest_out,
    PDELab::Concept::LocalMutableContainer             auto& lresidual_in,
    PDELab::Concept::LocalMutableContainer             auto& lresidual_out) noexcept;

  void localAssembleJacobianSkeleton(
    const Dune::Concept::Intersection               auto& intersection,
                                                    auto  time,
    const PDELab::Concept::LocalBasis               auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer      auto& llin_point_in,
    const PDELab::Concept::LocalBasis               auto& ltest_in,
    const PDELab::Concept::LocalBasis               auto& ltrial_out,
    const PDELab::Concept::LocalConstContainer      auto& llin_point_out,
    const PDELab::Concept::LocalBasis               auto& ltest_out,
    PDELab::Concept::LocalMutableMatrix             auto& ljacobian_in_in,
    PDELab::Concept::LocalMutableMatrix             auto& ljacobian_in_out,
    PDELab::Concept::LocalMutableMatrix             auto& ljacobian_out_in,
    PDELab::Concept::LocalMutableMatrix             auto& ljacobian_out_out) noexcept;

  void localAssembleJacobianSkeletonApply(
    const Dune::Concept::Intersection                  auto& intersection,
                                                       auto  time,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer         auto& llin_point_in,
    const PDELab::Concept::LocalConstContainer         auto& lapp_point_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    const PDELab::Concept::LocalBasis                  auto& ltrial_out,
    const PDELab::Concept::LocalConstContainer         auto& llin_point_out,
    const PDELab::Concept::LocalConstContainer         auto& lapp_point_out,
    const PDELab::Concept::LocalBasis                  auto& ltest_out,
    PDELab::Concept::LocalMutableContainer             auto& ljacobian_in,
    PDELab::Concept::LocalMutableContainer             auto& ljacobian_out) noexcept;

  void localAssemblePatternSkeleton(
    const Dune::Concept::Intersection                  auto& intersection,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    const PDELab::Concept::LocalBasis                  auto& ltrial_out,
    const PDELab::Concept::LocalBasis                  auto& ltest_out,
                                                       auto& lpattern_in_in,
                                                       auto& lpattern_in_out,
                                                       auto& lpattern_out_in,
                                                       auto& lpattern_out_out) noexcept;

  constexpr static auto localAssembleDoBoundary() {
    return std::false_type{};
  }

  void localAssembleBoundary(
    const Dune::Concept::Intersection                  auto& intersection,
                                                       auto  time,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer         auto& lcoefficients_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    PDELab::Concept::LocalMutableContainer             auto& lresidual_in) noexcept;

  void localAssembleJacobianBoundary(
    const Dune::Concept::Intersection                  auto& intersection,
                                                       auto  time,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer         auto& llin_point_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    PDELab::Concept::LocalMutableContainer             auto& ljacobian_in) noexcept;

  void localAssembleJacobianBoundaryApply(
    const Dune::Concept::Intersection                  auto& intersection,
                                                       auto  time,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer         auto& llin_point_in,
    const PDELab::Concept::LocalConstContainer         auto& lapp_point_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    PDELab::Concept::LocalMutableContainer             auto& ljacobian_in) noexcept;

  void localAssemblePatternBoundary(
    const Dune::Concept::Intersection                  auto& intersection,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
                                                       auto& lpattern_in) noexcept;

};

} // namespace Dune::PDELab::inline Experimental::LocalAssembly

#endif // DUNE_PDELAB_OPERATOR_LOCAL_ASSEMBLY_ARCHETYPE_HH
