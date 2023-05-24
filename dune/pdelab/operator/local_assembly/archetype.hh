#ifndef DUNE_PDELAB_OPERATOR_LOCAL_ASSEMBLY_ARCHETYPE_HH
#define DUNE_PDELAB_OPERATOR_LOCAL_ASSEMBLY_ARCHETYPE_HH

#include <dune/pdelab/concepts/container.hh>
#include <dune/pdelab/concepts/local_basis.hh>

#include <dune/grid/concepts.hh>

#include <type_traits>

namespace Dune::PDELab::inline Experimental::LocalAssembly {

struct Archetype {

  constexpr static auto localAssembleIsLinear() {
    return std::false_type{};
  }

  constexpr static auto localAssembleSkipEntity(const Dune::Concept::Entity auto& entity) {
    return std::false_type{};
  }

  constexpr static auto localAssemblySkipIntersection(const Dune::Concept::Intersection auto& intersection) {
    return std::false_type{};
  }

  constexpr static auto localAssembleDoVolume() {
    return std::false_type{};
  }

  void localAssembleVolume(
                                                       auto  time_point,
    const PDELab::Concept::LocalBasis                  auto& ltrial,
    const PDELab::Concept::LocalConstContainer         auto& lcoefficients,
    const PDELab::Concept::LocalBasis                  auto& ltest,
    PDELab::Concept::LocalMutableContainer             auto& lresidual);

  void localAssembleJacobianVolume(
                                                       auto  time_point,
    const PDELab::Concept::LocalBasis                  auto& ltrial,
    const PDELab::Concept::LocalConstContainer         auto& llin_point,
    const PDELab::Concept::LocalBasis                  auto& ltest,
    PDELab::Concept::LocalMutableMatrix                auto& ljacobian);

  void localAssembleJacobianVolumeApply(
                                                       auto  time_point,
    const PDELab::Concept::LocalBasis                  auto& ltrial,
    const PDELab::Concept::LocalConstContainer         auto& llin_point,
    const PDELab::Concept::LocalConstContainer         auto& lapp_point,
    const PDELab::Concept::LocalBasis                  auto& ltest,
    PDELab::Concept::LocalMutableContainer             auto& ljacobian);

  void localAssemblePatternVolume(
    const PDELab::Concept::LocalBasis                  auto& ltrial,
    const PDELab::Concept::LocalBasis                  auto& ltest,
                                                       auto& lpattern);

  constexpr static auto localAssembleDoSkeleton() {
    return std::false_type{};
  }

  void localAssembleSkeleton(
    const Dune::Concept::Intersection                  auto& intersection,
                                                       auto  time_point,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer         auto& lcoefficients_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    const PDELab::Concept::LocalBasis                  auto& ltrial_out,
    const PDELab::Concept::LocalConstContainer         auto& lcoefficients_out,
    const PDELab::Concept::LocalBasis                  auto& ltest_out,
    PDELab::Concept::LocalMutableContainer             auto& lresidual_in,
    PDELab::Concept::LocalMutableContainer             auto& lresidual_out);

  void localAssembleJacobianSkeleton(
    const Dune::Concept::Intersection                  auto& intersection,
                                                       auto  time_point,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer         auto& llin_point_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    const PDELab::Concept::LocalBasis                  auto& ltrial_out,
    const PDELab::Concept::LocalConstContainer         auto& llin_point_out,
    const PDELab::Concept::LocalBasis                  auto& ltest_out,
    PDELab::Concept::LocalMutableContainer             auto& ljacobian_in,
    PDELab::Concept::LocalMutableContainer             auto& ljacobian_out);

  void localAssembleJacobianSkeletonApply(
    const Dune::Concept::Intersection                  auto& intersection,
                                                       auto  time_point,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer         auto& llin_point_in,
    const PDELab::Concept::LocalConstContainer         auto& lapp_point_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    const PDELab::Concept::LocalBasis                  auto& ltrial_out,
    const PDELab::Concept::LocalConstContainer         auto& llin_point_out,
    const PDELab::Concept::LocalConstContainer         auto& lapp_point_out,
    const PDELab::Concept::LocalBasis                  auto& ltest_out,
    PDELab::Concept::LocalMutableContainer             auto& ljacobian_in,
    PDELab::Concept::LocalMutableContainer             auto& ljacobian_out);

  void localAssemblePatternSkeleton(
    const Dune::Concept::Intersection                  auto& intersection,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    const PDELab::Concept::LocalBasis                  auto& ltrial_out,
    const PDELab::Concept::LocalBasis                  auto& ltest_out,
                                                       auto& lpattern_in_in,
                                                       auto& lpattern_in_out,
                                                       auto& lpattern_out_in,
                                                       auto& lpattern_out_out);

  constexpr static auto localAssembleDoBoundary() {
    return std::false_type{};
  }

  void localAssembleBoundary(
    const Dune::Concept::Intersection                  auto& intersection,
                                                       auto  time_point,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer         auto& lcoefficients_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    PDELab::Concept::LocalMutableContainer             auto& lresidual_in);

  void localAssembleJacobianBoundary(
    const Dune::Concept::Intersection                  auto& intersection,
                                                       auto  time_point,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer         auto& llin_point_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    PDELab::Concept::LocalMutableContainer             auto& ljacobian_in);

  void localAssembleJacobianBoundaryApply(
    const Dune::Concept::Intersection                  auto& intersection,
                                                       auto  time_point,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer         auto& llin_point_in,
    const PDELab::Concept::LocalConstContainer         auto& lapp_point_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
    PDELab::Concept::LocalMutableContainer             auto& ljacobian_in);

  void localAssemblePatternBoundary(
    const Dune::Concept::Intersection                  auto& intersection,
    const PDELab::Concept::LocalBasis                  auto& ltrial_in,
    const PDELab::Concept::LocalBasis                  auto& ltest_in,
                                                       auto& lpattern_in);

};

} // namespace Dune::PDELab::inline Experimental::LocalAssembly

#endif // DUNE_PDELAB_OPERATOR_LOCAL_ASSEMBLY_ARCHETYPE_HH
