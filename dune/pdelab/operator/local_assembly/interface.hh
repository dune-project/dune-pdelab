#ifndef DUNE_ASSEMBLER_OPERATOR_LOCAL_ASSEMLBY_INTERFACE_HH
#define DUNE_ASSEMBLER_OPERATOR_LOCAL_ASSEMLBY_INTERFACE_HH

#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/concepts/container.hh>

#include <dune/grid/concepts/intersection.hh>
#include <dune/grid/concepts/entity.hh>

namespace Dune::PDELab::inline Experimental::LocalAssembly {

namespace Impl {

inline namespace Default {

constexpr auto localAssembleSkipEntity(const auto& lop, const Dune::Concept::Entity auto& entity) {
  if constexpr (requires { { lop.localAssembleSkipEntity(entity) } -> std::convertible_to<bool>; })
    return lop.localAssembleSkipEntity(entity);
  else
    return std::false_type();
}

constexpr auto localAssembleSkipIntersection(const auto& lop, const Dune::Concept::Intersection auto& intersection) {
  if constexpr (requires { { lop.localAssembleSkipIntersection(intersection) } -> std::convertible_to<bool>; })
    return lop.localAssembleSkipIntersection(intersection);
  else
    return std::false_type();
}

constexpr auto localAssembleIsLinear(const auto& lop) {
  if constexpr (requires { { lop.localAssembleIsLinear() } -> std::convertible_to<bool>; })
    return lop.localAssembleIsLinear();
  else
    return std::false_type();
}

constexpr auto localAssembleDoVolume(const auto& lop) {
  if constexpr (requires { { lop.localAssembleDoVolume() } -> std::convertible_to<bool>; })
    return lop.localAssembleDoVolume();
  else
    return std::false_type();
}

void localAssembleVolume(
                                                  auto& lop,
                                                  auto  time_point,
  const PDELab::Concept::LocalBasis               auto& ltrial,
  const PDELab::Concept::LocalConstContainer      auto& lcoefficients,
  const PDELab::Concept::LocalBasis               auto& ltest,
  PDELab::Concept::LocalMutableContainer          auto& lresidual)
{
  lop.localAssembleVolume(time_point, ltrial, lcoefficients, ltest, lresidual);
}

void localAssemblePatternVolume(
                                                  auto& lop,
  const PDELab::Concept::LocalBasis               auto& ltrial,
  const PDELab::Concept::LocalBasis               auto& ltest,
                                                  auto& lpattern)
{
  lop.localAssemblePatternVolume(ltrial, ltest, lpattern);
}

void localAssembleJacobianVolume(
                                                  auto& lop,
                                                  auto  time_point,
  const PDELab::Concept::LocalBasis               auto& ltrial,
  const PDELab::Concept::LocalConstContainer      auto& llin_point,
  const PDELab::Concept::LocalBasis               auto& ltest,
  PDELab::Concept::LocalMutableMatrix             auto& ljacobian)
{
  lop.localAssembleJacobianVolume(time_point, ltrial, llin_point, ltest, ljacobian);
}

void localAssembleJacobianVolumeApply(
                                                  auto& lop,
                                                  auto  time_point,
  const PDELab::Concept::LocalBasis               auto& ltrial,
  const PDELab::Concept::LocalConstContainer      auto& llin_point,
  const PDELab::Concept::LocalConstContainer      auto& lpoint,
  const PDELab::Concept::LocalBasis               auto& ltest,
  PDELab::Concept::LocalMutableContainer          auto& ljacobian)
{
  lop.localAssembleJacobianVolumeApply(time_point, ltrial, llin_point, lpoint, ltest, ljacobian);
}

constexpr auto localAssembleDoSkeleton(const auto& lop) {
  if constexpr (requires { { lop.localAssembleDoSkeleton() } -> std::convertible_to<bool>; })
    return lop.localAssembleDoSkeleton();
  else
    return std::false_type();
}

void localAssembleSkeleton(
                                                  auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
                                                  auto  time_point,
  const PDELab::Concept::LocalBasis               auto& ltrial_in,
  const PDELab::Concept::LocalConstContainer      auto& lcoefficients_in,
  const PDELab::Concept::LocalBasis               auto& ltest_in,
  const PDELab::Concept::LocalBasis               auto& ltrial_out,
  const PDELab::Concept::LocalConstContainer      auto& lcoefficients_out,
  const PDELab::Concept::LocalBasis               auto& ltest_out,
  PDELab::Concept::LocalMutableContainer          auto& lresidual_in,
  PDELab::Concept::LocalMutableContainer          auto& lresidual_out)
{
  lop.localAssembleSkeleton(
    intersection, time_point,
    ltrial_in, lcoefficients_in, ltest_in,
    ltrial_out, lcoefficients_out, ltest_out,
    lresidual_in,
    lresidual_out);
}

void localAssemblePatternSkeleton(
                                          auto& lop,
  const Dune::Concept::Intersection       auto& intersection,
  const PDELab::Concept::LocalBasis       auto& ltrial_in,
  const PDELab::Concept::LocalBasis       auto& ltest_in,
  const PDELab::Concept::LocalBasis       auto& ltrial_out,
  const PDELab::Concept::LocalBasis       auto& ltest_out,
                                          auto& lpattern_in_in,
                                          auto& lpattern_in_out,
                                          auto& lpattern_out_in,
                                          auto& lpattern_out_out)
{
  lop.localAssemblePatternSkeleton(
    intersection,
    ltrial_in, ltest_in,
    ltrial_out, ltest_out,
    lpattern_in_in,
    lpattern_in_out,
    lpattern_out_in,
    lpattern_out_out);
}

void localAssembleJacobianSkeleton(
                                                  auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
                                                  auto  time_point,
  const PDELab::Concept::LocalBasis               auto& ltrial_in,
  const PDELab::Concept::LocalConstContainer      auto& llin_point_in,
  const PDELab::Concept::LocalBasis               auto& ltest_in,
  const PDELab::Concept::LocalBasis               auto& ltrial_out,
  const PDELab::Concept::LocalConstContainer      auto& llin_point_out,
  const PDELab::Concept::LocalBasis               auto& ltest_out,
  PDELab::Concept::LocalMutableMatrix             auto& ljacobian_in_in,
  PDELab::Concept::LocalMutableMatrix             auto& ljacobian_in_out,
  PDELab::Concept::LocalMutableMatrix             auto& ljacobian_out_in,
  PDELab::Concept::LocalMutableMatrix             auto& ljacobian_out_out)
{
  lop.localAssembleJacobianSkeleton(
    intersection, time_point,
    ltrial_in, llin_point_in, ltest_in,
    ltrial_out, llin_point_out, ltest_out,
    ljacobian_in_in, ljacobian_in_out,
    ljacobian_out_in, ljacobian_out_out);
}

void localAssembleJacobianSkeletonApply(
                                                  auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
                                                  auto  time_point,
  const PDELab::Concept::LocalBasis               auto& ltrial_in,
  const PDELab::Concept::LocalConstContainer      auto& llin_point_in,
  const PDELab::Concept::LocalConstContainer      auto& lpoint_in,
  const PDELab::Concept::LocalBasis               auto& ltest_in,
  const PDELab::Concept::LocalBasis               auto& ltrial_out,
  const PDELab::Concept::LocalConstContainer      auto& llin_point_out,
  const PDELab::Concept::LocalConstContainer      auto& lpoint_out,
  const PDELab::Concept::LocalBasis               auto& ltest_out,
  PDELab::Concept::LocalMutableContainer          auto& ljacobian_in,
  PDELab::Concept::LocalMutableContainer          auto& ljacobian_out)
{
  lop.localAssembleJacobianSkeletonApply(
    intersection, time_point,
    ltrial_in, llin_point_in, lpoint_in, ltest_in,
    ltrial_out, llin_point_out, lpoint_out, ltest_out,
    ljacobian_in,
    ljacobian_out);
}

constexpr auto localAssembleDoBoundary(const auto& lop) {
  if constexpr (requires { { lop.localAssembleDoBoundary() } -> std::convertible_to<bool>; })
    return lop.localAssembleDoBoundary();
  else
    return std::false_type();
}

void localAssembleBoundary(
                                                  auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
                                                  auto  time_point,
  const PDELab::Concept::LocalBasis               auto& ltrial_in,
  const PDELab::Concept::LocalConstContainer      auto& lcoefficients_in,
  const PDELab::Concept::LocalBasis               auto& ltest_in,
                                                  auto& lresidual_in)
{
  lop.localAssembleBoundary(
    intersection, time_point,
    ltrial_in, lcoefficients_in, ltest_in,
    lresidual_in);
}

void localAssemblePatternBoundary(
                                                  auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
  const PDELab::Concept::LocalBasis               auto& ltrial_in,
  const PDELab::Concept::LocalBasis               auto& ltest_in,
                                                  auto& lpattern_in_in)
{
  lop.localAssemblePatternBoundary(
    intersection,
    ltrial_in, ltest_in,
    lpattern_in_in);
}

void localAssembleJacobianBoundary(
                                                  auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
                                                  auto  time_point,
  const PDELab::Concept::LocalBasis               auto& ltrial_in,
  const PDELab::Concept::LocalConstContainer      auto& llin_point_in,
  const PDELab::Concept::LocalBasis               auto& ltest_in,
  PDELab::Concept::LocalMutableContainer          auto& ljacobian_in)
{
  lop.localAssembleJacobianBoundary(
    intersection, time_point,
    ltrial_in, llin_point_in, ltest_in,
    ljacobian_in);
}

void localAssembleJacobianBoundaryApply(
                                                  auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
                                                  auto  time_point,
  const PDELab::Concept::LocalBasis               auto& ltrial_in,
  const PDELab::Concept::LocalConstContainer      auto& llin_point_in,
  const PDELab::Concept::LocalConstContainer      auto& lpoint_in,
  const PDELab::Concept::LocalBasis               auto& ltest_in,
  PDELab::Concept::LocalMutableContainer          auto& ljacobian_in)
{
  lop.localAssembleJacobianBoundaryApply(
    intersection, time_point,
    ltrial_in, llin_point_in, lpoint_in, ltest_in,
    ljacobian_in);
}

} // namespace Default

struct SkipEntity {
  constexpr std::convertible_to<bool> auto operator()(const auto& lop, const Dune::Concept::Entity auto& entity) const {
    return localAssembleSkipEntity(lop, entity);
  }
};

struct SkipIntersection {
  constexpr std::convertible_to<bool> auto operator()(const auto& lop, const Dune::Concept::Intersection auto& intersection) const {
    return localAssembleSkipIntersection(lop, intersection);
  }
};

struct IsLinear {
  constexpr std::convertible_to<bool> auto operator()(const auto& lop) const {
    return localAssembleIsLinear(lop);
  }
};

struct DoVolume {
  constexpr std::convertible_to<bool> auto operator()(const auto& lop) const {
    return localAssembleDoVolume(lop);
  }

  static void apply(auto& lop, auto&& f) {
    // static/dynamic switch
    const auto do_assemble = localAssembleDoVolume(lop);
    if constexpr (requires { std::bool_constant<do_assemble>{}; } ) {
      if constexpr (do_assemble)
        f(lop);
    } else if (do_assemble) {
      f(lop);
    }
  }
};

struct Volume {
  void operator()(
                                                    auto& lop,
                                                    auto  time_point,
    const PDELab::Concept::LocalBasis               auto& ltrial,
    const PDELab::Concept::LocalConstContainer      auto& lcoefficients,
    const PDELab::Concept::LocalBasis               auto& ltest,
    PDELab::Concept::LocalMutableContainer          auto& lresidual) const
  requires requires { { localAssembleDoVolume(lop) } -> std::convertible_to<bool>; }
  {
    DoVolume::apply(lop, [&](auto& _lop){
      localAssembleVolume(_lop, time_point, ltrial, lcoefficients, ltest, lresidual);
    });
  }
};

struct PatternVolume {
  void operator()(
                                                    auto& lop,
    const PDELab::Concept::LocalBasis               auto& ltrial,
    const PDELab::Concept::LocalBasis               auto& ltest,
                                                    auto& lpattern) const
  requires requires { { localAssembleDoVolume(lop) } -> std::convertible_to<bool>; }
  {
    DoVolume::apply(lop, [&](auto& _lop){
      localAssemblePatternVolume(_lop, ltrial, ltest, lpattern);
    });
  }
};

struct JacobianVolume {
  void operator()(
                                                    auto& lop,
                                                    auto  time_point,
    const PDELab::Concept::LocalBasis               auto& ltrial,
    const PDELab::Concept::LocalConstContainer      auto& llin_point,
    const PDELab::Concept::LocalBasis               auto& ltest,
    PDELab::Concept::LocalMutableMatrix             auto& ljacobian) const
  {
    DoVolume::apply(lop, [&](auto& _lop){
      localAssembleJacobianVolume(_lop, time_point, ltrial, llin_point, ltest, ljacobian);
    });
  }
};

struct JacobianVolumeApply {
  void operator()(
                                                    auto& lop,
                                                    auto  time_point,
    const PDELab::Concept::LocalBasis               auto& ltrial,
    const PDELab::Concept::LocalConstContainer      auto& llin_point,
    const PDELab::Concept::LocalConstContainer      auto& lpoint,
    const PDELab::Concept::LocalBasis               auto& ltest,
    PDELab::Concept::LocalMutableContainer          auto& ljacobian) const
  {
    DoVolume::apply(lop, [&](auto& _lop){
      localAssembleJacobianVolumeApply(_lop, time_point, ltrial, llin_point, lpoint, ltest, ljacobian);
    });
  }
};

struct DoSkeleton {
  constexpr std::convertible_to<bool> auto operator()(const auto& lop) const {
    return localAssembleDoSkeleton(lop);
  }

  static void apply(auto& lop, auto&& f) {
    // static/dynamic switch
    const auto do_assemble = localAssembleDoSkeleton(lop);
    if constexpr (requires { std::bool_constant<do_assemble>{}; } ) {
      if constexpr (do_assemble)
        f(lop);
    } else if (do_assemble) {
      f(lop);
    }
  }
};

struct Skeleton {
  void operator()(
                                                    auto& lop,
    const Dune::Concept::Intersection               auto& intersection,
                                                    auto  time_point,
    const PDELab::Concept::LocalBasis               auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer      auto& lcoefficients_in,
    const PDELab::Concept::LocalBasis               auto& ltest_in,
    const PDELab::Concept::LocalBasis               auto& ltrial_out,
    const PDELab::Concept::LocalConstContainer      auto& lcoefficients_out,
    const PDELab::Concept::LocalBasis               auto& ltest_out,
    PDELab::Concept::LocalMutableContainer          auto& lresidual_in,
    PDELab::Concept::LocalMutableContainer          auto& lresidual_out) const
  {
    DoSkeleton::apply(lop, [&](auto& _lop){
      localAssembleSkeleton(_lop, intersection, time_point,
              ltrial_in,  lcoefficients_in, ltest_in,
              ltrial_out, lcoefficients_out, ltest_out,
              lresidual_in, lresidual_out);
    });
  }
};

struct PatternSkeleton {
  void operator()(
                                                    auto& lop,
    const Dune::Concept::Intersection               auto& intersection,
    const PDELab::Concept::LocalBasis               auto& ltrial_in,
    const PDELab::Concept::LocalBasis               auto& ltest_in,
    const PDELab::Concept::LocalBasis               auto& ltrial_out,
    const PDELab::Concept::LocalBasis               auto& ltest_out,
                                                    auto& lresidual_in_in,
                                                    auto& lresidual_in_out,
                                                    auto& lresidual_out_in,
                                                    auto& lresidual_out_out) const
  {
    DoSkeleton::apply(lop, [&](auto& _lop){
      localAssemblePatternSkeleton(_lop, intersection,
                ltrial_in, ltest_in,
                ltrial_out, ltest_out,
                lresidual_in_in, lresidual_in_out,
                lresidual_out_in, lresidual_out_out);
    });
  }
};

struct JacobianSkeleton {
  void operator()(
                                                    auto& lop,
    const Dune::Concept::Intersection               auto& intersection,
                                                    auto  time_point,
    const PDELab::Concept::LocalBasis               auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer      auto& llin_point_in,
    const PDELab::Concept::LocalBasis               auto& ltest_in,
    const PDELab::Concept::LocalBasis               auto& ltrial_out,
    const PDELab::Concept::LocalConstContainer      auto& llin_point_out,
    const PDELab::Concept::LocalBasis               auto& ltest_out,
    PDELab::Concept::LocalMutableMatrix             auto& ljacobian_in_in,
    PDELab::Concept::LocalMutableMatrix             auto& ljacobian_in_out,
    PDELab::Concept::LocalMutableMatrix             auto& ljacobian_out_in,
    PDELab::Concept::LocalMutableMatrix             auto& ljacobian_out_out) const
  {
    DoSkeleton::apply(lop, [&](auto& _lop){
      localAssembleJacobianSkeleton(lop, intersection, time_point,
                  ltrial_in, llin_point_in, ltest_in,
                  ltrial_out, llin_point_out, ltest_out,
                  ljacobian_in_in, ljacobian_in_out,
                  ljacobian_out_in, ljacobian_out_out);
    });
  }
};

struct JacobianSkeletonApply {
  void operator()(
                                                    auto& lop,
    const Dune::Concept::Intersection               auto& intersection,
                                                    auto  time_point,
    const PDELab::Concept::LocalBasis               auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer      auto& llin_point_in,
    const PDELab::Concept::LocalConstContainer      auto& lpoint_in,
    const PDELab::Concept::LocalBasis               auto& ltest_in,
    const PDELab::Concept::LocalBasis               auto& ltrial_out,
    const PDELab::Concept::LocalConstContainer      auto& llin_point_out,
    const PDELab::Concept::LocalConstContainer      auto& lpoint_out,
    const PDELab::Concept::LocalBasis               auto& ltest_out,
    PDELab::Concept::LocalMutableContainer          auto& ljacobian_in,
    PDELab::Concept::LocalMutableContainer          auto& ljacobian_out) const
  {
    DoSkeleton::apply(lop, [&](auto& _lop){
      localAssembleJacobianSkeletonApply(_lop, intersection, time_point,
                  ltrial_in, llin_point_in, lpoint_in, ltest_in,
                  ltrial_out, llin_point_out, lpoint_out, ltest_out,
                  ljacobian_in, ljacobian_out);
    });
  }
};

struct DoBoundary {
  constexpr std::convertible_to<bool> auto operator()(const auto& lop) const {
    return localAssembleDoBoundary(lop);
  }

  static void apply(auto& lop, auto&& f) {
    // static/dynamic switch
    const auto do_assemble = localAssembleDoBoundary(lop);
    if constexpr (requires { std::bool_constant<do_assemble>{}; } ) {
      if constexpr (do_assemble)
        f(lop);
    } else if (do_assemble) {
      f(lop);
    }
  }
};


struct Boundary {
  void operator()(
                                                    auto& lop,
    const Dune::Concept::Intersection               auto& intersection,
                                                    auto  time_point,
    const PDELab::Concept::LocalBasis               auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer      auto& lcoefficients_in,
    const PDELab::Concept::LocalBasis               auto& ltest_in,
    PDELab::Concept::LocalMutableContainer          auto& lresidual_in) const
  {
    DoBoundary::apply(lop, [&](auto& _lop){
      localAssembleBoundary(_lop, intersection, time_point, ltrial_in, lcoefficients_in, ltest_in, lresidual_in);
    });
  }
};

struct PatternBoundary {
  void operator()(
                                                    auto& lop,
    const Dune::Concept::Intersection               auto& intersection,
    const PDELab::Concept::LocalBasis               auto& ltrial_in,
    const PDELab::Concept::LocalBasis               auto& ltest_in,
                                                    auto& lpattern_in_in) const
  {
    DoBoundary::apply(lop, [&](auto& _lop){
      localAssemblePatternBoundary(_lop, intersection, ltrial_in, ltest_in, lpattern_in_in);
    });
  }
};


struct JacobianBoundary {
  void operator()(
                                                    auto& lop,
    const Dune::Concept::Intersection               auto& intersection,
                                                    auto  time_point,
    const PDELab::Concept::LocalBasis               auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer      auto& llin_point_in,
    const PDELab::Concept::LocalBasis               auto& ltest_in,
    PDELab::Concept::LocalMutableMatrix             auto& ljacobian_in) const
  {
    DoBoundary::apply(lop, [&](auto& _lop){
      localAssembleJacobianBoundary(_lop, intersection, time_point, ltrial_in, llin_point_in, ltest_in, ljacobian_in);
    });
  }
};

struct JacobianBoundaryApply {
  void operator()(
                                                    auto& lop,
    const Dune::Concept::Intersection               auto& intersection,
                                                    auto  time_point,
    const PDELab::Concept::LocalBasis               auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer      auto& llin_point_in,
    const PDELab::Concept::LocalConstContainer      auto& lpoint_in,
    const PDELab::Concept::LocalBasis               auto& ltest_in,
    PDELab::Concept::LocalMutableContainer          auto& ljacobian_in) const
  {
    DoBoundary::apply(lop, [&](auto& _lop){
      localAssembleJacobianBoundaryApply(lop, intersection, time_point, ltrial_in, llin_point_in, lpoint_in, ltest_in, ljacobian_in);
    });
  }
};

} // namespace Impl

namespace Default = Impl::Default;

//! Customization point objects
inline const Impl::SkipEntity             skipEntity{};
inline const Impl::SkipIntersection       skipIntersection{};
inline const Impl::IsLinear               isLinear{};

inline const Impl::DoVolume               doVolume{};
inline const Impl::Volume                 volume{};
inline const Impl::PatternVolume          patternVolume{};
inline const Impl::JacobianVolume         jacobianVolume{};
inline const Impl::JacobianVolumeApply    jacobianVolumeApply{};

inline const Impl::DoSkeleton             doSkeleton{};
inline const Impl::Skeleton               skeleton{};
inline const Impl::PatternSkeleton        patternSkeleton{};
inline const Impl::JacobianSkeleton       jacobianSkeleton{};
inline const Impl::JacobianSkeletonApply  jacobianSkeletonApply{};

inline const Impl::DoBoundary             doBoundary{};
inline const Impl::Boundary               boundary{};
inline const Impl::PatternBoundary        patternBoundary{};
inline const Impl::JacobianBoundary       jacobianBoundary{};
inline const Impl::JacobianBoundaryApply  jacobianBoundaryApply{};


} // namespace Dune::PDELab::inline Experimental::LocalAssembly

#endif // DUNE_ASSEMBLER_OPERATOR_LOCAL_ASSEMLBY_INTERFACE_HH
