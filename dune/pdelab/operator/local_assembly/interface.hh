#ifndef DUNE_ASSEMBLER_OPERATOR_LOCAL_ASSEMLBY_INTERFACE_HH
#define DUNE_ASSEMBLER_OPERATOR_LOCAL_ASSEMLBY_INTERFACE_HH

#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/concepts/container.hh>

#include <dune/grid/concepts/intersection.hh>
#include <dune/grid/concepts/entity.hh>

namespace Dune::PDELab::inline Experimental::LocalAssembly {

namespace Impl {

inline namespace Default {

constexpr auto localAssembleSkipEntity(const auto& lop, const Dune::Concept::Entity auto& entity)
requires requires { { lop.localAssembleSkipEntity(entity) } -> std::convertible_to<bool>; }
{
  return lop.localAssembleSkipEntity(entity);
}

constexpr auto localAssembleSkipIntersection(const auto& lop, const Dune::Concept::Intersection auto& intersection)
requires requires { { lop.localAssembleSkipIntersection(intersection) } -> std::convertible_to<bool>; }
{
  return lop.localAssembleSkipIntersection(intersection);
}

constexpr auto localAssembleIsLinear(const auto& lop)
requires requires { { lop.localAssembleIsLinear() } -> std::convertible_to<bool>; }
{
  return lop.localAssembleIsLinear();
}

constexpr auto localAssembleDoVolume(const auto& lop)
requires requires { { lop.localAssembleDoVolume() } -> std::convertible_to<bool>; }
{
  return lop.localAssembleDoVolume();
}

void localAssembleVolume(                         auto& lop,
                                                  auto  time,
  const PDELab::Concept::LocalBasis               auto& ltrial,
  const PDELab::Concept::LocalConstContainer      auto& lcoefficients,
  const PDELab::Concept::LocalBasis               auto& ltest,
  PDELab::Concept::LocalMutableContainer          auto& lresidual)
requires requires { lop.localAssembleVolume(time, ltrial, lcoefficients, ltest, lresidual); }
{
  lop.localAssembleVolume(time, ltrial, lcoefficients, ltest, lresidual);
}

void localAssemblePatternVolume(                  auto& lop,
  const PDELab::Concept::LocalBasis               auto& ltrial,
  const PDELab::Concept::LocalBasis               auto& ltest,
                                                  auto& lpattern)
requires requires { lop.localAssemblePatternVolume(ltrial, ltest, lpattern); }
{
  lop.localAssemblePatternVolume(ltrial, ltest, lpattern);
}

void localAssembleJacobianVolume(                 auto& lop,
                                                  auto  time,
  const PDELab::Concept::LocalBasis               auto& ltrial,
  const PDELab::Concept::LocalConstContainer      auto& llin_point,
  const PDELab::Concept::LocalBasis               auto& ltest,
  PDELab::Concept::LocalMutableMatrix             auto& ljacobian)
requires requires { lop.localAssembleJacobianVolume(time, ltrial, llin_point, ltest, ljacobian); }
{
  lop.localAssembleJacobianVolume(time, ltrial, llin_point, ltest, ljacobian);
}

void localAssembleJacobianVolumeApply(            auto& lop,
                                                  auto  time,
  const PDELab::Concept::LocalBasis               auto& ltrial,
  const PDELab::Concept::LocalConstContainer      auto& llin_point,
  const PDELab::Concept::LocalConstContainer      auto& lpoint,
  const PDELab::Concept::LocalBasis               auto& ltest,
  PDELab::Concept::LocalMutableContainer          auto& ljacobian)
requires requires { lop.localAssembleJacobianVolumeApply(time, ltrial, llin_point, lpoint, ltest, ljacobian); }
{
  lop.localAssembleJacobianVolumeApply(time, ltrial, llin_point, lpoint, ltest, ljacobian);
}

constexpr auto localAssembleDoSkeleton(const auto& lop)
requires requires { { lop.localAssembleDoSkeleton() } -> std::convertible_to<bool>; }
{
  return lop.localAssembleDoSkeleton();
}

void localAssembleSkeleton(                       auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
                                                  auto  time,
  const PDELab::Concept::LocalBasis               auto& ltrial_in,
  const PDELab::Concept::LocalConstContainer      auto& lcoefficients_in,
  const PDELab::Concept::LocalBasis               auto& ltest_in,
  const PDELab::Concept::LocalBasis               auto& ltrial_out,
  const PDELab::Concept::LocalConstContainer      auto& lcoefficients_out,
  const PDELab::Concept::LocalBasis               auto& ltest_out,
  PDELab::Concept::LocalMutableContainer          auto& lresidual_in,
  PDELab::Concept::LocalMutableContainer          auto& lresidual_out)
requires requires { lop.localAssembleSkeleton(intersection, time, ltrial_in, lcoefficients_in, ltest_in, ltrial_out, lcoefficients_out, ltest_out, lresidual_in, lresidual_out); }
{
  lop.localAssembleSkeleton(intersection, time, ltrial_in, lcoefficients_in, ltest_in, ltrial_out, lcoefficients_out, ltest_out, lresidual_in, lresidual_out);
}

void localAssemblePatternSkeleton(        auto& lop,
  const Dune::Concept::Intersection       auto& intersection,
  const PDELab::Concept::LocalBasis       auto& ltrial_in,
  const PDELab::Concept::LocalBasis       auto& ltest_in,
  const PDELab::Concept::LocalBasis       auto& ltrial_out,
  const PDELab::Concept::LocalBasis       auto& ltest_out,
                                          auto& lpattern_in_in,
                                          auto& lpattern_in_out,
                                          auto& lpattern_out_in,
                                          auto& lpattern_out_out)
requires requires {  lop.localAssemblePatternSkeleton(intersection, ltrial_in, ltest_in, ltrial_out, ltest_out, lpattern_in_in, lpattern_in_out, lpattern_out_in, lpattern_out_out); }
{
  lop.localAssemblePatternSkeleton(intersection, ltrial_in, ltest_in, ltrial_out, ltest_out, lpattern_in_in, lpattern_in_out, lpattern_out_in, lpattern_out_out);
}

void localAssembleJacobianSkeleton(                auto& lop,
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
  PDELab::Concept::LocalMutableMatrix             auto& ljacobian_out_out)
requires requires { lop.localAssembleJacobianSkeleton(intersection, time, ltrial_in, llin_point_in, ltest_in, ltrial_out, llin_point_out, ltest_out, ljacobian_in_in, ljacobian_in_out, ljacobian_out_in, ljacobian_out_out); }
{
  lop.localAssembleJacobianSkeleton(intersection, time, ltrial_in, llin_point_in, ltest_in, ltrial_out, llin_point_out, ltest_out, ljacobian_in_in, ljacobian_in_out, ljacobian_out_in, ljacobian_out_out);
}

void localAssembleJacobianSkeletonApply(          auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
                                                  auto  time,
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
requires requires { lop.localAssembleJacobianSkeletonApply(intersection, time, ltrial_in, llin_point_in, lpoint_in, ltest_in, ltrial_out, llin_point_out, lpoint_out, ltest_out, ljacobian_in, ljacobian_out); }
{
  lop.localAssembleJacobianSkeletonApply(intersection, time, ltrial_in, llin_point_in, lpoint_in, ltest_in, ltrial_out, llin_point_out, lpoint_out, ltest_out, ljacobian_in, ljacobian_out);
}

constexpr auto localAssembleDoBoundary(const auto& lop)
requires requires { { lop.localAssembleDoBoundary() } -> std::convertible_to<bool>; }
{
  return lop.localAssembleDoBoundary();
}

void localAssembleBoundary(                       auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
                                                  auto  time,
  const PDELab::Concept::LocalBasis               auto& ltrial_in,
  const PDELab::Concept::LocalConstContainer      auto& lcoefficients_in,
  const PDELab::Concept::LocalBasis               auto& ltest_in,
                                                  auto& lresidual_in)
requires requires { lop.localAssembleBoundary(intersection, time, ltrial_in, lcoefficients_in, ltest_in, lresidual_in); }
{
  lop.localAssembleBoundary(intersection, time, ltrial_in, lcoefficients_in, ltest_in, lresidual_in);
}

void localAssemblePatternBoundary(                auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
  const PDELab::Concept::LocalBasis               auto& ltrial_in,
  const PDELab::Concept::LocalBasis               auto& ltest_in,
                                                  auto& lpattern_in_in)
requires requires { lop.localAssemblePatternBoundary(intersection, ltrial_in, ltest_in, lpattern_in_in); }
{
  lop.localAssemblePatternBoundary(intersection, ltrial_in, ltest_in, lpattern_in_in);
}

void localAssembleJacobianBoundary(               auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
                                                  auto  time,
  const PDELab::Concept::LocalBasis               auto& ltrial_in,
  const PDELab::Concept::LocalConstContainer      auto& llin_point_in,
  const PDELab::Concept::LocalBasis               auto& ltest_in,
  PDELab::Concept::LocalMutableMatrix             auto& ljacobian_in)
requires requires { lop.localAssembleJacobianBoundary(intersection, time, ltrial_in, llin_point_in, ltest_in, ljacobian_in); }
{
  lop.localAssembleJacobianBoundary(intersection, time, ltrial_in, llin_point_in, ltest_in, ljacobian_in);
}

void localAssembleJacobianBoundaryApply(          auto& lop,
  const Dune::Concept::Intersection               auto& intersection,
                                                  auto  time,
  const PDELab::Concept::LocalBasis               auto& ltrial_in,
  const PDELab::Concept::LocalConstContainer      auto& llin_point_in,
  const PDELab::Concept::LocalConstContainer      auto& lpoint_in,
  const PDELab::Concept::LocalBasis               auto& ltest_in,
  PDELab::Concept::LocalMutableContainer          auto& ljacobian_in)
requires requires { lop.localAssembleJacobianBoundaryApply(intersection, time, ltrial_in, llin_point_in, lpoint_in, ltest_in, ljacobian_in); }
{
  lop.localAssembleJacobianBoundaryApply(intersection, time, ltrial_in, llin_point_in, lpoint_in, ltest_in, ljacobian_in);
}

} // namespace Default


struct SkipEntity {
  constexpr auto operator()(const auto& lop, const Dune::Concept::Entity auto& entity) const noexcept
  requires requires { { localAssembleSkipEntity(lop, entity) } -> std::convertible_to<bool>; } {
    return localAssembleSkipEntity(lop, entity);
  }

  constexpr std::false_type operator()(const auto& lop, const Dune::Concept::Entity auto& entity) const noexcept
  requires (!requires { localAssembleSkipEntity(lop, entity); }) {
    return {};
  }
};

struct SkipIntersection {
  constexpr auto operator()(const auto& lop, const Dune::Concept::Intersection auto& intersection) const noexcept
  requires requires { { localAssembleSkipIntersection(lop, intersection) } -> std::convertible_to<bool>; } {
    return localAssembleSkipIntersection(lop, intersection);
  }

  constexpr std::false_type operator()(const auto& lop, const Dune::Concept::Intersection auto& intersection) const noexcept
  requires (!requires { localAssembleSkipIntersection(lop, intersection); }) {
    return {};
  }
};

struct IsLinear {
  constexpr auto operator()(const auto& lop) const noexcept
  requires requires { { localAssembleIsLinear(lop) } -> std::convertible_to<bool>; } {
    return localAssembleIsLinear(lop);
  }

  constexpr std::false_type operator()(const auto& lop) const noexcept
  requires (!requires { localAssembleIsLinear(lop); }) {
    return {};
  }
};

struct DoVolume {
  constexpr auto operator()(const auto& lop) const noexcept
  requires requires { { localAssembleDoVolume(lop) } -> std::convertible_to<bool>; }
  {
    return localAssembleDoVolume(lop);
  }

  constexpr std::false_type operator()(const auto& lop) const noexcept
  requires (!requires { localAssembleDoVolume(lop); }) {
    return {};
  }

  static void and_then(auto& lop, auto&& f)
  {
    // static/dynamic switch
    const auto do_assemble = DoVolume{}(lop);
    if constexpr (requires { std::bool_constant<do_assemble>{}; } ) {
      if constexpr (do_assemble)
        f(lop);
    } else if (do_assemble) {
      f(lop);
    }
  }
};

struct Volume {
  template<class... Args>
  void operator()(auto& lop, Args&&... args) const noexcept
  {
    DoVolume::and_then(lop, [&](auto& _lop){
      localAssembleVolume(_lop, std::forward<Args>(args)...);
    });
  }
};

struct PatternVolume {
  template<class... Args>
  void operator()(auto& lop, Args&&... args) const noexcept
  {
    DoVolume::and_then(lop, [&](auto& _lop){
      localAssemblePatternVolume(_lop, std::forward<Args>(args)...);
    });
  }
};

struct JacobianVolume {
  template<class... Args>
  void operator()(auto& lop, Args&&... args) const noexcept
  {
    DoVolume::and_then(lop, [&](auto& _lop){
      localAssembleJacobianVolume(_lop, std::forward<Args>(args)...);
    });
  }
};

struct JacobianVolumeApply {
  template<class... Args>
  void operator()(auto& lop, Args&&... args) const noexcept
  {
    DoVolume::and_then(lop, [&](auto& _lop){
      localAssembleJacobianVolumeApply(_lop, std::forward<Args>(args)...);
    });
  }
};

struct DoSkeleton {
  constexpr auto operator()(const auto& lop) const noexcept
  requires requires { { localAssembleDoSkeleton(lop) } -> std::convertible_to<bool>; }
  {
    return localAssembleDoSkeleton(lop);
  }

  constexpr std::false_type operator()(const auto& lop) const noexcept
  requires (!requires { localAssembleDoSkeleton(lop); }) {
    return {};
  }

  static void and_then(auto& lop, auto&& f) {
    // static/dynamic switch
    const auto do_assemble = DoSkeleton{}(lop);
    if constexpr (requires { std::bool_constant<do_assemble>{}; } ) {
      if constexpr (do_assemble)
        f(lop);
    } else if (do_assemble) {
      f(lop);
    }
  }
};

struct Skeleton {
  template<class... Args>
  void operator()(auto& lop, Args&&... args) const noexcept
  {
    DoSkeleton::and_then(lop, [&](auto& _lop){
      localAssembleSkeleton(_lop, std::forward<Args>(args)...);
    });
  }
};

struct PatternSkeleton {
  template<class... Args>
  void operator()(auto& lop, Args&&... args) const noexcept
  {
    DoSkeleton::and_then(lop, [&](auto& _lop){
      localAssemblePatternSkeleton(_lop, std::forward<Args>(args)...);
    });
  }
};

struct JacobianSkeleton {
  template<class... Args>
  void operator()(auto& lop, Args&&... args) const noexcept
  {
    DoSkeleton::and_then(lop, [&](auto& _lop){
      localAssembleJacobianSkeleton(_lop, std::forward<Args>(args)...);
    });
  }
};

struct JacobianSkeletonApply {
  template<class... Args>
  void operator()(auto& lop, Args&&... args) const noexcept
  {
    DoSkeleton::and_then(lop, [&](auto& _lop){
      localAssembleJacobianSkeletonApply(_lop, std::forward<Args>(args)...);
    });
  }
};

struct DoBoundary {

  constexpr auto operator()(const auto& lop) const noexcept
  requires requires { { localAssembleDoBoundary(lop) } -> std::convertible_to<bool>; } {
    return localAssembleDoBoundary(lop);
  }

  constexpr std::false_type operator()(const auto& lop) const noexcept
  requires (!requires { localAssembleDoBoundary(lop); }) {
    return {};
  }

  static void and_then(auto& lop, auto&& f) {
    // static/dynamic switch
    const auto do_assemble = DoBoundary{}(lop);
    if constexpr (requires { std::bool_constant<do_assemble>{}; } ) {
      if constexpr (do_assemble)
        f(lop);
    } else if (do_assemble) {
      f(lop);
    }
  }
};


struct Boundary {
  template<class... Args>
  void operator()(auto& lop, Args&&... args) const noexcept
  {
    DoBoundary::and_then(lop, [&](auto& _lop){
      localAssembleBoundary(_lop, std::forward<Args>(args)...);
    });
  }
};

struct PatternBoundary {
  template<class... Args>
  void operator()(auto& lop, Args&&... args) const noexcept
  {
    DoBoundary::and_then(lop, [&](auto& _lop){
      localAssemblePatternBoundary(_lop, std::forward<Args>(args)...);
    });
  }
};


struct JacobianBoundary {
  template<class... Args>
  void operator()(auto& lop, Args&&... args) const noexcept
  {
    DoBoundary::and_then(lop, [&](auto& _lop){
      localAssembleJacobianBoundary(_lop, std::forward<Args>(args)...);
    });
  }
};

struct JacobianBoundaryApply {
  template<class... Args>
  void operator()(auto& lop, Args&&... args) const noexcept
{
    DoBoundary::and_then(lop, [&](auto& _lop){
      localAssembleJacobianBoundaryApply(_lop, std::forward<Args>(args)...);
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
