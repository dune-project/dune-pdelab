// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_CALLSWITCH_HH
#define DUNE_PDELAB_LOCALOPERATOR_CALLSWITCH_HH

#include <dune/common/concept.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
  namespace PDELab {

    namespace Impl {

    /** \internal */

    // compile time switching of function call
    template<typename LOP, bool doIt, bool isLinear = LOP::isLinear>
    struct LocalAssemblerCallSwitchHelper
    {
      //================
      // Pattern methods
      //================
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_volume (const LOP& lop, const LFSU& lfsu, const LFSV& lfsv, LocalPattern& pattern)
      {
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_volume_post_skeleton
      ( const LOP& lop,
        const LFSU& lfsu, const LFSV& lfsv,
        LocalPattern& pattern)
      {
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_skeleton (const LOP& lop, const LFSU& lfsu_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const LFSV& lfsv_n,
        LocalPattern& pattern_sn,
        LocalPattern& pattern_ns)
      {
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_boundary(const LOP& lop,
        const LFSU& lfsu_s, const LFSV& lfsv_s,
        LocalPattern& pattern_ss)
      {
      }

      //==============
      // Alpha methods
      //==============
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_volume (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
      {
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_volume_post_skeleton (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
      {
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_skeleton (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        R& r_s, R& r_n)
      {
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_boundary (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        R& r_s)
      {
      }

      //===============
      // Lambda methods
      //===============
      template<typename EG, typename LFSV, typename R>
      static void lambda_volume (const LOP& lop, const EG& eg, const LFSV& lfsv, R& r)
      {
      }
      template<typename EG, typename LFSV, typename R>
      static void lambda_volume_post_skeleton (const LOP& lop, const EG& eg, const LFSV& lfsv, R& r)
      {
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_skeleton(const LOP& lop, const IG& ig,
        const LFSV& lfsv_s, const LFSV& lfsv_n,
        R& r_s, R& r_n)
      {
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_boundary (const LOP& lop, const IG& ig, const LFSV& lfsv, R& r)
      {
      }

      //=======================
      // Jacobian apply methods
      //=======================
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_volume (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y)
      {
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_volume_post_skeleton (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y)
      {
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_skeleton (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const X& z_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n)
      {
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_boundary (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        Y& y_s)
      {
      }

      //=================
      // Jacobian methods
      //=================
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_volume (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat)
      {
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_volume_post_skeleton (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M& mat)
      {
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_skeleton (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        M & mat_ss, M & mat_sn,
        M & mat_ns, M & mat_nn)
      {
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_boundary (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        M & mat_ss)
      {
      }
    };


    template<typename LOP>
    struct LocalAssemblerCallSwitchHelper<LOP,true,true>
    {
      //================
      // Pattern methods
      //================
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_volume (const LOP& lop, const LFSU& lfsu, const LFSV& lfsv, LocalPattern& pattern)
      {
        lop.pattern_volume(lfsu,lfsv,pattern);
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_volume_post_skeleton
      ( const LOP& lop,
        const LFSU& lfsu, const LFSV& lfsv,
        LocalPattern& pattern)
      {
        lop.pattern_volume_post_skeleton(lfsu,lfsv,pattern);
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_skeleton (const LOP& lop, const LFSU& lfsu_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const LFSV& lfsv_n,
        LocalPattern& pattern_sn,
        LocalPattern& pattern_ns)
      {
        lop.pattern_skeleton(lfsu_s,lfsv_s,lfsu_n,lfsv_n,
          pattern_sn, pattern_ns);
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_boundary(const LOP& lop,
        const LFSU& lfsu_s, const LFSV& lfsv_s,
        LocalPattern& pattern_ss)
      {
        lop.pattern_boundary(lfsu_s,lfsv_s,pattern_ss);
      }

      //==============
      // Alpha methods
      //==============
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_volume (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
      {
        lop.alpha_volume(eg,lfsu,x,lfsv,r);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_volume_post_skeleton (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
      {
        lop.alpha_volume_post_skeleton(eg,lfsu,x,lfsv,r);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_skeleton (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        R& r_s, R& r_n)
      {
        lop.alpha_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,r_s,r_n);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_boundary (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        R& r_s)
      {
        lop.alpha_boundary(ig,lfsu_s,x_s,lfsv_s,r_s);
      }

      //===============
      // Lambda methods
      //===============
      template<typename EG, typename LFSV, typename R>
      static void lambda_volume (const LOP& lop, const EG& eg, const LFSV& lfsv, R& r)
      {
        lop.lambda_volume(eg,lfsv,r);
      }
      template<typename EG, typename LFSV, typename R>
      static void lambda_volume_post_skeleton (const LOP& lop, const EG& eg, const LFSV& lfsv, R& r)
      {
        lop.lambda_volume_post_skeleton(eg,lfsv,r);
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_skeleton(const LOP& lop, const IG& ig,
        const LFSV& lfsv_s, const LFSV& lfsv_n,
        R& r_s, R& r_n)
      {
        lop.lambda_skeleton(ig, lfsv_s, lfsv_n, r_s, r_n);
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_boundary (const LOP& lop, const IG& ig, const LFSV& lfsv, R& r)
      {
        lop.lambda_boundary(ig,lfsv,r);
      }

      //=======================
      // Jacobian apply methods
      //=======================
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      static auto jacobian_apply_volume (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y)
      {
        lop.jacobian_apply_volume(eg,lfsu,z,lfsv,y);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_volume_post_skeleton (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y)
      {
        lop.jacobian_apply_volume_post_skeleton(eg,lfsu,z,lfsv,y);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_skeleton (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const X& z_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n)
      {
        lop.jacobian_apply_skeleton(ig,lfsu_s,z_s,lfsv_s,lfsu_n,z_n,lfsv_n,y_s,y_n);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_boundary (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        Y& y_s)
      {
        lop.jacobian_apply_boundary(ig,lfsu_s,z_s,lfsv_s,y_s);
      }

      //=================
      // Jacobian methods
      //=================
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_volume (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat)
      {
        lop.jacobian_volume(eg,lfsu,x,lfsv,mat);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_volume_post_skeleton (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat)
      {
        lop.jacobian_volume_post_skeleton(eg,lfsu,x,lfsv,mat);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_skeleton (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        M & mat_ss, M & mat_sn,
        M & mat_ns, M & mat_nn)
      {
        lop.jacobian_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,
          mat_ss, mat_sn, mat_ns, mat_nn);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_boundary (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        M & mat_ss)
      {
        lop.jacobian_boundary(ig,lfsu_s,x_s,lfsv_s,mat_ss);
      }
    };

    /** \brief Specialization of LocalAssemblerCallSwitchHelper for the non-linear case

        We derive from the linear case and only replace those calls which differ, namely all jacobian_apply* calls.
        Compared to the linear case these methods now also pass the linearization point x to the local operator.
     */
    template<typename LOP>
    struct LocalAssemblerCallSwitchHelper<LOP,true,false> :
      public LocalAssemblerCallSwitchHelper<LOP,true,true>
    {
      //=======================
      // Jacobian apply methods
      //=======================
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      static auto jacobian_apply_volume (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y)
      {
        lop.jacobian_apply_volume(eg,lfsu,x,z,lfsv,y);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_volume_post_skeleton (const LOP& lop, const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y)
      {
        lop.jacobian_apply_volume_post_skeleton(eg,lfsu,x,z,lfsv,y);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_skeleton (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const X& z_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n)
      {
        lop.jacobian_apply_skeleton(ig,lfsu_s,x_s,z_s,lfsv_s,lfsu_n,x_n,z_n,lfsv_n,y_s,y_n);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_boundary (const LOP& lop, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        Y& y_s)
      {
        lop.jacobian_apply_boundary(ig,lfsu_s,x_s,z_s,lfsv_s,y_s);
      }
    };

    } // namespace Impl

    /** compile time switching of function call */
    template<typename LOP, bool doIt>
    using LocalAssemblerCallSwitch =
      Impl::LocalAssemblerCallSwitchHelper<LOP,doIt>;

    namespace LocalOperatorApply {

      auto patternVolume = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doPatternVolume>::
          pattern_volume(lop, args...);
      };

      auto patternVolumePostSkeleton = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doPatternVolumePostSkeleton>::
          pattern_volume_post_skeleton(lop, args...);
      };

      auto patternSkeleton = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doPatternSkeleton>::
          pattern_skeleton(lop, args...);
      };

      auto patternBoundary = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doPatternBoundary>::
          pattern_boundary(lop, args...);
      };

      //////////////////////
      auto alphaVolume = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doAlphaVolume>::
          alpha_volume(lop, args...);
      };

      auto alphaVolumePostSkeleton = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doAlphaVolumePostSkeleton>::
          alpha_volume_post_skeleton(lop, args...);
      };

      auto alphaSkeleton = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doAlphaSkeleton>::
          alpha_skeleton(lop, args...);
      };

      auto alphaBoundary = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doAlphaBoundary>::
          alpha_boundary(lop, args...);
      };


      //////////////////////
      auto lambdaVolume = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doLambdaVolume>::
          lambda_volume(lop, args...);
      };

      auto lambdaVolumePostSkeleton = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doLambdaVolumePostSkeleton>::
          lambda_volume_post_skeleton(lop, args...);
      };

      auto lambdaSkeleton = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doLambdaSkeleton>::
          lambda_skeleton(lop, args...);
      };

      auto lambdaBoundary = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doLambdaBoundary>::
          lambda_boundary(lop, args...);
      };


      //////////////////////
      auto jacobianVolume = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doAlphaVolume>::
          jacobian_volume(lop, args...);
      };

      auto jacobianVolumePostSkeleton = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doAlphaVolumePostSkeleton>::
          jacobian_volume_post_skeleton(lop, args...);
      };

      auto jacobianSkeleton = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doAlphaSkeleton>::
          jacobian_skeleton(lop, args...);
      };

      auto jacobianBoundary = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doAlphaBoundary>::
          jacobian_boundary(lop, args...);
      };


      //////////////////////
      auto jacobianApplyVolume = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doAlphaVolume>::
          jacobian_apply_volume(lop, args...);
      };

      auto jacobianApplyVolumePostSkeleton = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doAlphaVolumePostSkeleton>::
          jacobian_apply_volume_post_skeleton(lop, args...);
      };

      auto jacobianApplySkeleton = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doAlphaSkeleton>::
          jacobian_apply_skeleton(lop, args...);
      };

      auto jacobianApplyBoundary = [](const auto& lop, auto&... args)
      {
        using LOP = std::decay_t<decltype(lop)>;
        Impl::LocalAssemblerCallSwitchHelper<LOP,LOP::doAlphaBoundary>::
          jacobian_apply_boundary(lop, args...);
      };

    } // namespace LocalOperatorApply

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_CALLSWITCH_HH
