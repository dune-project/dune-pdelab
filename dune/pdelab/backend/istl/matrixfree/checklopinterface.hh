// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_CHECKLOPINTERFACE_HH
#define DUNE_PDELAB_BACKEND_ISTL_MATRIXFREE_CHECKLOPINTERFACE_HH


namespace Dune {
  namespace PDELab {
    namespace impl{
      // Some of the matrix free solvers require that the local operator describing
      // the model uses the new interface for nonlinear Jacobian apply methods where
      // the current solution and the linearization point may have different
      // type. This is best explained by an example:
      //
      //
      // old: void jacobian_apply_volume (const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y) const {}
      // new: void jacobian_apply_volume (const EG& eg, const LFSU& lfsu, const X& x, const Z& z, const LFSV& lfsv, Y& y) const {}
      //
      // In the new interface the coefficient z has type Z which may be different
      // than X.
      //
      // The purpose of this file is to provide a helper function to check whether a
      // local operator complies to the new interface. I'm all in favor of removing
      // this code as soon as possible.
      //
      // The current implementation is based on C++-14 SFINAE. It is entirely
      // possible that this could be replaced with features from the current
      // concepts implementation in dune-common.

      // SFINAE Infrastructure
      template <typename T>
      struct Result
      {
      private:
        template <typename A> constexpr auto test(int)
          -> decltype(std::declval<T>()(std::declval<A>()), std::true_type())
        {
          return std::true_type();
        }
        template <typename A> constexpr std::false_type test(...)
        {
          return std::false_type();
        }

      public:
        template <typename A> constexpr auto operator()(const A& p)
        {
          return test<A>(int());
        }
      };

      template <typename T> constexpr auto lambdaToTemplate(const T& t)
      {
        return Result<T>();
      }

      // The actual tests if the new method exists
      //
      // Note: The types in the declval don't really matter. The important part
      // is that they are all different. This way we check that we have the new
      // interface since it would not be valid to pass this to the old one.
      //
      // Note: Don't forget the & for the last argument. This is necessary
      // since the last argument of jacobian_apply_volume is passed by
      // reference (without being const).
      auto hasNewJacobianApplyVolume =
        lambdaToTemplate([](auto&& lop) -> decltype(lop.jacobian_apply_volume(std::declval<bool>(),
                                                                              std::declval<int>(),
                                                                              std::declval<double>(),
                                                                              std::declval<short int>(),
                                                                              std::declval<long int>(),
                                                                              std::declval<long long int&>()
                                                      )) {});
      auto hasNewJacobianApplyVolumePostSkeleton =
        lambdaToTemplate([](auto&& lop) -> decltype(lop.jacobian_apply_volume_post_skeleton(std::declval<bool>(),
                                                                                            std::declval<int>(),
                                                                                            std::declval<double>(),
                                                                                            std::declval<short int>(),
                                                                                            std::declval<long int>(),
                                                                                            std::declval<long long int&>()
                                                      )) {});
      auto hasNewJacobianApplySkeleton =
        lambdaToTemplate([](auto&& lop) -> decltype(lop.jacobian_apply_skeleton(std::declval<bool>(),
                                                                                std::declval<int>(),
                                                                                std::declval<double>(),
                                                                                std::declval<short int>(),
                                                                                std::declval<long int>(),
                                                                                std::declval<int>(),
                                                                                std::declval<double>(),
                                                                                std::declval<short int>(),
                                                                                std::declval<long int>(),
                                                                                std::declval<long long int&>(),
                                                                                std::declval<long long int&>()
                                                      )) {});
      auto hasNewJacobianApplyBoundary =
        lambdaToTemplate([](auto&& lop) -> decltype(lop.jacobian_apply_boundary(std::declval<bool>(),
                                                                                std::declval<int>(),
                                                                                std::declval<double>(),
                                                                                std::declval<short int>(),
                                                                                std::declval<long int>(),
                                                                                std::declval<long long int&>()
                                                      )) {});

      // Similar tests that show if the old or the new methods exist
      //
      // Note: Since the calls to the old interface are also valid in the new
      // interface we only get the information if the method exists in one of
      // the two implementations.
      auto hasOldOrNewJacobianApplyVolume =
        lambdaToTemplate([](auto&& lop) -> decltype(lop.jacobian_apply_volume(std::declval<bool>(),
                                                                              std::declval<int>(),
                                                                              std::declval<double>(),
                                                                              std::declval<double>(),
                                                                              std::declval<long int>(),
                                                                              std::declval<long long int&>()
                                                      )) {});
      auto hasOldOrNewJacobianApplyVolumePostSkeleton =
        lambdaToTemplate([](auto&& lop) -> decltype(lop.jacobian_apply_volume_post_skeleton(std::declval<bool>(),
                                                                                            std::declval<int>(),
                                                                                            std::declval<double>(),
                                                                                            std::declval<double>(),
                                                                                            std::declval<long int>(),
                                                                                            std::declval<long long int&>()
                                                      )) {});
      auto hasOldOrNewJacobianApplySkeleton =
        lambdaToTemplate([](auto&& lop) -> decltype(lop.jacobian_apply_skeleton(std::declval<bool>(),
                                                                                std::declval<int>(),
                                                                                std::declval<double>(),
                                                                                std::declval<double>(),
                                                                                std::declval<long int>(),
                                                                                std::declval<int>(),
                                                                                std::declval<double>(),
                                                                                std::declval<double>(),
                                                                                std::declval<long int>(),
                                                                                std::declval<long long int&>(),
                                                                                std::declval<long long int&>()
                                                      )) {});
      auto hasOldOrNewJacobianApplyBoundary =
        lambdaToTemplate([](auto&& lop) -> decltype(lop.jacobian_apply_boundary(std::declval<bool>(),
                                                                                std::declval<int>(),
                                                                                std::declval<double>(),
                                                                                std::declval<double>(),
                                                                                std::declval<long int>(),
                                                                                std::declval<long long int&>()
                                                      )) {});

      // Now comes the test that should actually be called. If a localoperator
      // implements the OldOrNew interface for one of the nonlinear jacobian
      // apply methods but not the new interface we can conlude that only the
      // old one is implemented.
      template <typename T>
      constexpr auto hasOldLOPInterface(T& t) -> typename std::enable_if<
        (decltype(hasOldOrNewJacobianApplyVolume(t))::value && !decltype(hasNewJacobianApplyVolume(t))::value)
        || (decltype(hasOldOrNewJacobianApplyVolumePostSkeleton(t))::value && !decltype(hasNewJacobianApplyVolumePostSkeleton(t))::value)
        || (decltype(hasOldOrNewJacobianApplySkeleton(t))::value && !decltype(hasNewJacobianApplySkeleton(t))::value)
        || (decltype(hasOldOrNewJacobianApplyBoundary(t))::value && !decltype(hasNewJacobianApplyBoundary(t))::value),
        std::true_type>::type
      {
        return std::true_type();
      }
      template <typename T>
      constexpr auto hasOldLOPInterface(T& t) -> typename std::enable_if<
        !((decltype(hasOldOrNewJacobianApplyVolume(t))::value && !decltype(hasNewJacobianApplyVolume(t))::value)
          || (decltype(hasOldOrNewJacobianApplyVolumePostSkeleton(t))::value && !decltype(hasNewJacobianApplyVolumePostSkeleton(t))::value)
          || (decltype(hasOldOrNewJacobianApplySkeleton(t))::value && !decltype(hasNewJacobianApplySkeleton(t))::value)
          || (decltype(hasOldOrNewJacobianApplyBoundary(t))::value && !decltype(hasNewJacobianApplyBoundary(t))::value)),
        std::false_type>::type
      {
        return std::false_type();
      }
    }
  }
}


#endif
