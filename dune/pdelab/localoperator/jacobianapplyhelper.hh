// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=4 sts=4:

#ifndef DUNE_PDELAB_LOCALOPERATOR_JACOBIANAPPLYHELPER_HH
#define DUNE_PDELAB_LOCALOPERATOR_JACOBIANAPPLYHELPER_HH

// This file implements helper functions to call the jacobian_apply_* methods
// on a local operator. Since the interface of linear and nonlinear operators
// looks different this helper functions are useful in some places.
//
// Note: This should be removed as soon as a good replacement is implemented in
// callswitch.hh! For this reason it is put into an impl namespace and it is
// not intended to be used in user code.

namespace Dune{
  namespace PDELab{
    namespace impl{
      //======================
      // Jacobian Apply Volume
      //======================

      template <typename LOP, typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      std::enable_if_t<LOP::isLinear> jacobianApplyVolume(
          const LOP& lop,
          const EG& eg,
          const LFSU& lfsu, const X& z, const LFSV& lfsv,
          Y& y)
      {
          lop.jacobian_apply_volume(eg, lfsu, z, lfsv, y);
      }
        template <typename LOP, typename EG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
        std::enable_if_t<LOP::isLinear> jacobianApplyVolume(
            const LOP& lop,
            const EG& eg,
            const LFSU& lfsu, const X& x, const Z& z, const LFSV& lfsv,
            Y& y)
        {
            DUNE_THROW(Dune::Exception, "You try to call a nonlinear method on a linear operator");
        }
        template <typename LOP, typename EG, typename LFSU, typename X, typename LFSV, typename Y>
        std::enable_if_t<not LOP::isLinear> jacobianApplyVolume(
            const LOP& lop,
            const EG& eg,
            const LFSU& lfsu, const X& z, const LFSV& lfsv,
            Y& y)
        {
            DUNE_THROW(Dune::Exception, "You try to call a linear method on a nonlinear operator");
        }
        template <typename LOP, typename EG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
        std::enable_if_t<not LOP::isLinear> jacobianApplyVolume(
            const LOP& lop,
            const EG& eg,
            const LFSU& lfsu, const X& x, const Z& z, const LFSV& lfsv,
            Y& y)
        {
            lop.jacobian_apply_volume(eg, lfsu, x, z, lfsv, y);
        }

        //========================
        // Jacobian Apply Skeleton
        //========================

        template <typename LOP, typename IG, typename LFSU, typename X, typename LFSV, typename Y>
        std::enable_if_t<LOP::isLinear> jacobianApplySkeleton(
            const LOP& lop,
            const IG& ig,
            const LFSU& lfsu_s, const X& z_s, const LFSV& lfsv_s,
            const LFSU& lfsu_n, const X& z_n, const LFSV& lfsv_n,
            Y& y_s, Y& y_n)
        {
            lop.jacobian_apply_skeleton(ig, lfsu_s, z_s, lfsv_s, lfsu_n, z_n, lfsv_n, y_s, y_n);
        }
        template <typename LOP, typename IG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
        std::enable_if_t<LOP::isLinear> jacobianApplySkeleton(
            const LOP& lop,
            const IG& ig,
            const LFSU& lfsu_s, const X& x_s, const Z& z_s, const LFSV& lfsv_s,
            const LFSU& lfsu_n, const X& x_n, const Z& z_n, const LFSV& lfsv_n,
            Y& y_s, Y& y_n)
        {
            DUNE_THROW(Dune::Exception, "You try to call a nonlinear method on a linear operator");
        }
        template <typename LOP, typename IG, typename LFSU, typename X, typename LFSV, typename Y>
        std::enable_if_t<not LOP::isLinear> jacobianApplySkeleton(
            const LOP& lop,
            const IG& ig,
            const LFSU& lfsu_s, const X& z_s, const LFSV& lfsv_s,
            const LFSU& lfsu_n, const X& z_n, const LFSV& lfsv_n,
            Y& y_s, Y& y_n)
        {
            DUNE_THROW(Dune::Exception, "You try to call a nonlinear method on a linear operator");
        }
        template <typename LOP, typename IG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
        std::enable_if_t<not LOP::isLinear> jacobianApplySkeleton(
            const LOP& lop,
            const IG& ig,
            const LFSU& lfsu_s, const X& x_s, const Z& z_s, const LFSV& lfsv_s,
            const LFSU& lfsu_n, const X& x_n, const Z& z_n, const LFSV& lfsv_n,
            Y& y_s, Y& y_n)
        {
            lop.jacobian_apply_skeleton(ig, lfsu_s, x_s, z_s, lfsv_s, lfsu_n, x_n, z_n, lfsv_n, y_s, y_n);
        }

        //========================
        // Jacobian Apply Boundary
        //========================

        template<typename LOP, typename IG, typename LFSU, typename X, typename LFSV, typename Y>
        std::enable_if_t<LOP::isLinear> jacobianApplyBoundary(
            const LOP& lop,
            const IG& ig,
            const LFSU& lfsu_s, const X& z_s, const LFSV& lfsv_s,
            Y& y_s)
        {
            lop.jacobian_apply_boundary(ig, lfsu_s, z_s, lfsv_s, y_s);
        }
        template<typename LOP, typename IG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
        std::enable_if_t<LOP::isLinear> jacobianApplyBoundary(
            const LOP& lop,
            const IG& ig,
            const LFSU& lfsu_s, const X& x_s, const Z& z_s, const LFSV& lfsv_s,
            Y& y_s)
        {
            DUNE_THROW(Dune::Exception, "You try to call a nonlinear method on a linear operator");
        }
        template<typename LOP, typename IG, typename LFSU, typename X, typename LFSV, typename Y>
        std::enable_if_t<not LOP::isLinear> jacobianApplyBoundary(
            const LOP& lop,
            const IG& ig,
            const LFSU& lfsu_s, const X& z_s, const LFSV& lfsv_s,
            Y& y_s)
        {
            DUNE_THROW(Dune::Exception, "You try to call a nonlinear method on a linear operator");
        }
        template<typename LOP, typename IG, typename LFSU, typename X, typename Z, typename LFSV, typename Y>
        std::enable_if_t<not LOP::isLinear> jacobianApplyBoundary(
            const LOP& lop,
            const IG& ig,
            const LFSU& lfsu_s, const X& x_s, const Z& z_s, const LFSV& lfsv_s,
            Y& y_s)
        {
            lop.jacobian_apply_boundary(ig, lfsu_s, x_s, z_s, lfsv_s, y_s);
        }

    } //namespace impl
  } // namespace PDELab
} // namespace Dune


#endif
