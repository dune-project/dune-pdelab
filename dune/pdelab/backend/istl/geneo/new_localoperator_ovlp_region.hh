#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_NEW_LOCALOPERATOR_OVLP_REGION_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_NEW_LOCALOPERATOR_OVLP_REGION_HH

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

namespace Dune {
  namespace PDELab {

    /*!
     * \brief Wrapper for LocalOperators restricting their action to areas overlapping with other subdomains
     * Any entity not covered by multiple subdomains will be ignored by not evaluating the underlying LocalOperator.
     * This wrapper can be used to generate the overlap matrix needed for the GenEO coarse space.
     * \tparam LocalOperatorBase The LocalOperator to be wrapped
     * \tparam GFS The corresponding grid function space
     */
    template<typename LocalOperatorBase, typename Vector>
    class NewLocalOperatorOvlpRegion {

      //typedef Dune::PDELab::Backend::Vector<GFS,int> V;

    public:
      NewLocalOperatorOvlpRegion (LocalOperatorBase& base_, const Vector& partUnity)
       : baseop(base_), partUnity_(partUnity) {
      }

      // pattern assembly flags
      enum { doPatternVolume = LocalOperatorBase::doPatternVolume };
      enum { doPatternSkeleton = LocalOperatorBase::doPatternSkeleton };
      enum { doPatternVolumePostSkeleton = LocalOperatorBase::doPatternVolumePostSkeleton };
      enum { doPatternBoundary = LocalOperatorBase::doPatternBoundary };

      // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doAlphaVolumePostSkeleton = LocalOperatorBase::doAlphaVolumePostSkeleton };
      //enum { doAlphaPostSkeletonVolume = LocalOperatorBase::doAlphaPostSkeletonVolume };
      enum { doAlphaBoundary = true };

      enum { doAlphaSkeleton = LocalOperatorBase::doAlphaSkeleton };
      enum { doLambdaVolume  = LocalOperatorBase::doLambdaVolume };
      enum { doLambdaBoundary = LocalOperatorBase::doLambdaBoundary };
      enum { doLambdaSkeleton = LocalOperatorBase::doLambdaSkeleton };
      enum { doSkeletonTwoSided = LocalOperatorBase::doSkeletonTwoSided };

      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_volume (const LFSU& lfsu, const LFSV& lfsv, LocalPattern& pattern) const
      {
        if (entity_is_interior(lfsu))
          return;
        baseop.pattern_volume(lfsu, lfsv, pattern);
      }

      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_skeleton (const LFSU& lfsu_s, const LFSV& lfsv_s, const LFSU& lfsu_n, const LFSV& lfsv_n, LocalPattern& pattern_sn, LocalPattern& pattern_ns) const
      {
        if (entity_is_interior(lfsu_s))
          return;
        baseop.pattern_skeleton(lfsu_s, lfsv_s, lfsu_n, lfsv_n, pattern_sn, pattern_ns);
      }

      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_volume_post_skeleton (const LFSU& lfsu, const LFSV& lfsv, LocalPattern& pattern) const
      {
        if (entity_is_interior(lfsu))
          return;
        baseop.pattern_volume(lfsu, lfsv, pattern);
      }

      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_boundary (const LFSU& lfsu_s, const LFSV& lfsv_s, LocalPattern& pattern_ss) const
      {
        if (entity_is_interior(lfsu_s))
          return;
        baseop.pattern_boundary (lfsu_s, lfsv_s, pattern_ss);
      }

      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        if (entity_is_interior(lfsu))
          return;
        baseop.alpha_volume(eg, lfsu, x, lfsv, r);
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M& mat) const
      {
        if (entity_is_interior(lfsu))
          return;
        baseop.jacobian_volume (eg, lfsu, x, lfsv, mat);
      }

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        if (entity_is_interior(lfsv))
          return;
        baseop.lambda_volume (eg, lfsv, r);
      }

      // post skeleton: compute time step allowable for cell; to be done later
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                      const LFSV& lfsv, R& r) const
      {
        if (entity_is_interior(lfsu))
          return;
        baseop.alpha_volume_post_skeleton (eg, lfsu, x, lfsv, r);
      }

      // boundary integral
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        if (entity_is_interior(lfsu_s))
          return;
        baseop.alpha_boundary (ig, lfsu_s, x_s, lfsv_s, r_s);
      }

      // jacobian contribution from boundary
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              M& mat_s) const
      {
        if (entity_is_interior(lfsu_s))
          return;
        baseop.jacobian_boundary (ig, lfsu_s, x_s, lfsv_s, mat_s);
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        if (entity_is_interior(lfsu_s))
          return;
        baseop.alpha_skeleton (ig, lfsu_s, x_s, lfsv_s, lfsu_n, x_n, lfsv_n, r_s, r_n);
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                              M& mat_ss, M& mat_sn,
                              M& mat_ns, M& mat_nn) const
      {
        if (entity_is_interior(lfsu_s))
          return;
        baseop.jacobian_skeleton (ig, lfsu_s, x_s, lfsv_s, lfsu_n, x_n, lfsv_n, mat_ss, mat_sn, mat_ns, mat_nn);
      }


      void setTime (double t)
      {
        baseop.setTime(t);
      }

    private:

      template <typename LFS>
      bool entity_is_interior (const LFS& lfs) const {
        LFSIndexCache<LFS> cache(lfs);
        cache.update();
        for (std::size_t i = 0; i < cache.size(); i++)
          {
            if (partUnity_[cache.containerIndex(i)[0]] > 0.0 &&
                partUnity_[cache.containerIndex(i)[0]] < 1.0)
              return false;
          }
        return true;
      }

      LocalOperatorBase& baseop;
      const Vector& partUnity_;
    };

  }
}

#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_LOCALOPERATOR_OVLP_REGION_HH
