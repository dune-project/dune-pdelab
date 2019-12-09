// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_CALLSWITCH_HH
#define DUNE_PDELAB_LOCALOPERATOR_CALLSWITCH_HH

#include <dune/common/concept.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
  namespace PDELab {


    // compile time switching of function call
    template<typename LA, bool doIt, bool isLinear>
    struct LocalAssemblerCallSwitch
    {
      //================
      // Pattern methods
      //================
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_volume (const LA& la, const LFSU& lfsu, const LFSV& lfsv, LocalPattern& pattern)
      {
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_volume_post_skeleton
      ( const LA& la,
        const LFSU& lfsu, const LFSV& lfsv,
        LocalPattern& pattern)
      {
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_skeleton (const LA& la, const LFSU& lfsu_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const LFSV& lfsv_n,
        LocalPattern& pattern_sn,
        LocalPattern& pattern_ns)
      {
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_boundary(const LA& la,
        const LFSU& lfsu_s, const LFSV& lfsv_s,
        LocalPattern& pattern_ss)
      {
      }

      //==============
      // Alpha methods
      //==============
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
      {
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
      {
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_skeleton (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        R& r_s, R& r_n)
      {
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_boundary (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        R& r_s)
      {
      }

      //===============
      // Lambda methods
      //===============
      template<typename EG, typename LFSV, typename R>
      static void lambda_volume (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
      }
      template<typename EG, typename LFSV, typename R>
      static void lambda_volume_post_skeleton (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_skeleton(const LA& la, const IG& ig,
        const LFSV& lfsv_s, const LFSV& lfsv_n,
        R& r_s, R& r_n)
      {
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_boundary (const LA& la, const IG& ig, const LFSV& lfsv, R& r)
      {
      }

      //=======================
      // Jacobian apply methods
      //=======================
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y)
      {
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y)
      {
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_skeleton (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const X& z_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n)
      {
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_boundary (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        Y& y_s)
      {
      }

      //=================
      // Jacobian methods
      //=================
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat)
      {
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M& mat)
      {
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_skeleton (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        M & mat_ss, M & mat_sn,
        M & mat_ns, M & mat_nn)
      {
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_boundary (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        M & mat_ss)
      {
      }
    };


    template<typename LA>
    struct LocalAssemblerCallSwitch<LA,true,true>
    {
      //================
      // Pattern methods
      //================
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_volume (const LA& la, const LFSU& lfsu, const LFSV& lfsv, LocalPattern& pattern)
      {
        la.pattern_volume(lfsu,lfsv,pattern);
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_volume_post_skeleton
      ( const LA& la,
        const LFSU& lfsu, const LFSV& lfsv,
        LocalPattern& pattern)
      {
        la.pattern_volume_post_skeleton(lfsu,lfsv,pattern);
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_skeleton (const LA& la, const LFSU& lfsu_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const LFSV& lfsv_n,
        LocalPattern& pattern_sn,
        LocalPattern& pattern_ns)
      {
        la.pattern_skeleton(lfsu_s,lfsv_s,lfsu_n,lfsv_n,
          pattern_sn, pattern_ns);
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_boundary(const LA& la,
        const LFSU& lfsu_s, const LFSV& lfsv_s,
        LocalPattern& pattern_ss)
      {
        la.pattern_boundary(lfsu_s,lfsv_s,pattern_ss);
      }

      //==============
      // Alpha methods
      //==============
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
      {
        la.alpha_volume(eg,lfsu,x,lfsv,r);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
      {
        la.alpha_volume_post_skeleton(eg,lfsu,x,lfsv,r);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_skeleton (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        R& r_s, R& r_n)
      {
        la.alpha_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,r_s,r_n);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_boundary (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        R& r_s)
      {
        la.alpha_boundary(ig,lfsu_s,x_s,lfsv_s,r_s);
      }

      //===============
      // Lambda methods
      //===============
      template<typename EG, typename LFSV, typename R>
      static void lambda_volume (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
        la.lambda_volume(eg,lfsv,r);
      }
      template<typename EG, typename LFSV, typename R>
      static void lambda_volume_post_skeleton (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
        la.lambda_volume_post_skeleton(eg,lfsv,r);
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_skeleton(const LA& la, const IG& ig,
        const LFSV& lfsv_s, const LFSV& lfsv_n,
        R& r_s, R& r_n)
      {
        la.lambda_skeleton(ig, lfsv_s, lfsv_n, r_s, r_n);
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_boundary (const LA& la, const IG& ig, const LFSV& lfsv, R& r)
      {
        la.lambda_boundary(ig,lfsv,r);
      }

      //=======================
      // Jacobian apply methods
      //=======================
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      static auto jacobian_apply_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y)
      {
        la.jacobian_apply_volume(eg,lfsu,z,lfsv,y);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y)
      {
        la.jacobian_apply_volume_post_skeleton(eg,lfsu,z,lfsv,y);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_skeleton (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const X& z_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n)
      {
        la.jacobian_apply_skeleton(ig,lfsu_s,z_s,lfsv_s,lfsu_n,z_n,lfsv_n,y_s,y_n);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_boundary (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        Y& y_s)
      {
        la.jacobian_apply_boundary(ig,lfsu_s,z_s,lfsv_s,y_s);
      }

      //=================
      // Jacobian methods
      //=================
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat)
      {
        la.jacobian_volume(eg,lfsu,x,lfsv,mat);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat)
      {
        la.jacobian_volume_post_skeleton(eg,lfsu,x,lfsv,mat);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_skeleton (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        M & mat_ss, M & mat_sn,
        M & mat_ns, M & mat_nn)
      {
        la.jacobian_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,
          mat_ss, mat_sn, mat_ns, mat_nn);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_boundary (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        M & mat_ss)
      {
        la.jacobian_boundary(ig,lfsu_s,x_s,lfsv_s,mat_ss);
      }
    };


    template<typename LA>
    struct LocalAssemblerCallSwitch<LA,true,false>
    {
      //================
      // Pattern methods
      //================
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_volume (const LA& la, const LFSU& lfsu, const LFSV& lfsv, LocalPattern& pattern)
      {
        la.pattern_volume(lfsu,lfsv,pattern);
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_volume_post_skeleton
      ( const LA& la,
        const LFSU& lfsu, const LFSV& lfsv,
        LocalPattern& pattern)
      {
        la.pattern_volume_post_skeleton(lfsu,lfsv,pattern);
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_skeleton (const LA& la, const LFSU& lfsu_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const LFSV& lfsv_n,
        LocalPattern& pattern_sn,
        LocalPattern& pattern_ns)
      {
        la.pattern_skeleton(lfsu_s,lfsv_s,lfsu_n,lfsv_n,
          pattern_sn, pattern_ns);
      }
      template<typename LFSU, typename LFSV, typename LocalPattern>
      static void pattern_boundary(const LA& la,
        const LFSU& lfsu_s, const LFSV& lfsv_s,
        LocalPattern& pattern_ss)
      {
        la.pattern_boundary(lfsu_s,lfsv_s,pattern_ss);
      }

      //==============
      // Alpha methods
      //==============
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
      {
        la.alpha_volume(eg,lfsu,x,lfsv,r);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
      {
        la.alpha_volume_post_skeleton(eg,lfsu,x,lfsv,r);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_skeleton (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        R& r_s, R& r_n)
      {
        la.alpha_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,r_s,r_n);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      static void alpha_boundary (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        R& r_s)
      {
        la.alpha_boundary(ig,lfsu_s,x_s,lfsv_s,r_s);
      }

      //===============
      // Lambda methods
      //===============
      template<typename EG, typename LFSV, typename R>
      static void lambda_volume (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
        la.lambda_volume(eg,lfsv,r);
      }
      template<typename EG, typename LFSV, typename R>
      static void lambda_volume_post_skeleton (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
        la.lambda_volume_post_skeleton(eg,lfsv,r);
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_skeleton(const LA& la, const IG& ig,
        const LFSV& lfsv_s, const LFSV& lfsv_n,
        R& r_s, R& r_n)
      {
        la.lambda_skeleton(ig, lfsv_s, lfsv_n, r_s, r_n);
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_boundary (const LA& la, const IG& ig, const LFSV& lfsv, R& r)
      {
        la.lambda_boundary(ig,lfsv,r);
      }

      //=======================
      // Jacobian apply methods
      //=======================
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      static auto jacobian_apply_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y)
      {
        la.jacobian_apply_volume(eg,lfsu,x,z,lfsv,y);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv, Y& y)
      {
        la.jacobian_apply_volume_post_skeleton(eg,lfsu,x,z,lfsv,y);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_skeleton (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const X& z_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n)
      {
        la.jacobian_apply_skeleton(ig,lfsu_s,x_s,z_s,lfsv_s,lfsu_n,x_n,z_n,lfsv_n,y_s,y_n);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      static void jacobian_apply_boundary (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        Y& y_s)
      {
        la.jacobian_apply_boundary(ig,lfsu_s,x_s,z_s,lfsv_s,y_s);
      }

      //=================
      // Jacobian methods
      //=================
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat)
      {
        la.jacobian_volume(eg,lfsu,x,lfsv,mat);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat)
      {
        la.jacobian_volume_post_skeleton(eg,lfsu,x,lfsv,mat);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_skeleton (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        M & mat_ss, M & mat_sn,
        M & mat_ns, M & mat_nn)
      {
        la.jacobian_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,
          mat_ss, mat_sn, mat_ns, mat_nn);
      }
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      static void jacobian_boundary (const LA& la, const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        M & mat_ss)
      {
        la.jacobian_boundary(ig,lfsu_s,x_s,lfsv_s,mat_ss);
      }
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_CALLSWITCH_HH
