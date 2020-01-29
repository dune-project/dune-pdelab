// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_WEIGHTEDSUM_HH
#define DUNE_PDELAB_LOCALOPERATOR_WEIGHTEDSUM_HH

#include <cstddef>

#include <dune/common/forloop.hh>
#include <dune/common/fvector.hh>
#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/typetraits.hh>

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>

#include <dune/pdelab/localoperator/sum.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    //! A local operator to take the weighted sum of other local operators
    /**
     * \nosubgrouping
     *
     * If the weight for one summand is zero, calls to that local operators
     * evaluation and pattern methods are eliminated at run-time.
     *
     * \tparam K    Type of the scaling factors.
     * \tparam Args Tuple of local operators.  Must fulfill \c
     *              tuple_size<Args>::value>=1.
     */
    template<typename K, typename Args>
    class WeightedSumLocalOperator
    {
      static const std::size_t size = std::tuple_size<Args>::value;

      typedef typename ForEachType<AddPtrTypeEvaluator, Args>::Type ArgPtrs;
      typedef typename ForEachType<AddRefTypeEvaluator, Args>::Type ArgRefs;

      ArgPtrs lops;
      typedef FieldVector<K, size> Weights;
      Weights weights;

      template<typename F, typename... Args2>
      void applyLops(F && f, Args2 &&... args)
      {
        Hybrid::forEach(Std::make_index_sequence<size-1>{},
          [&](auto i){
            if(weights[i] != K(0))
              f(*Hybrid::elementAt(lops, i), std::forward<Args2>(args)...);});
      }

    public:
      //////////////////////////////////////////////////////////////////////
      //
      //! \name Construction and modification
      //! \{
      //

      //! construct a WeightedSumLocalOperator
      /**
       * No summand local operators are set.  They must be initialized with
       * setSummand() before the constructed object is used.
       */
      WeightedSumLocalOperator
      ( const Weights& weights_ = Weights(1))
        : weights(weights_)
      { }

      //! construct a WeightedSumLocalOperator from a tuple of local operators
      WeightedSumLocalOperator
      ( const ArgRefs& lops_,
        const Weights& weights_ = Weights(1))
        : lops(transformTuple<AddPtrTypeEvaluator>(lops_)), weights(weights_)
      { }

      //! set the i'th component of the sum
      template<std::size_t i>
      void setSummand(typename std::tuple_element<i,Args>::type& summand)
      { std::get<i>(lops) = &summand; }

      //! get the i'th component of the sum
      template<std::size_t i>
      typename std::tuple_element<i,Args>::type& getSummand()
      { return *std::get<i>(lops); }

      //! set the weight for the i'th component of the sum
      void setWeight(K w, std::size_t i)
      { weights[i] = w; }

      //! get the weight for the i'th component of the sum
      K getWeight(std::size_t i)
      { return weights[i]; }

      //! \} Construction and modification

      ////////////////////////////////////////////////////////////////////////
      //
      //! \name Control flags
      //! \{
      //

    private:
      template<typename T>
      using PatternVolumeValue = std::integral_constant<bool, T::doPatternVolume>;
      template<typename T>
      using PatternVolumePostSkeletonValue = std::integral_constant<bool, T::doPatternVolumePostSkeleton>;
      template<typename T>
      using PatternSkeletonValue = std::integral_constant<bool, T::doPatternSkeleton>;
      template<typename T>
      using PatternBoundaryValue = std::integral_constant<bool, T::doPatternBoundary>;

      template<typename T>
      using AlphaVolumeValue = std::integral_constant<bool, T::doAlphaVolume>;
      template<typename T>
      using AlphaVolumePostSkeletonValue = std::integral_constant<bool, T::doAlphaVolumePostSkeleton>;
      template<typename T>
      using AlphaSkeletonValue = std::integral_constant<bool, T::doAlphaSkeleton>;
      template<typename T>
      using AlphaBoundaryValue = std::integral_constant<bool, T::doAlphaBoundary>;

      template<typename T>
      using LambdaVolumeValue = std::integral_constant<bool, T::doLambdaVolume>;
      template<typename T>
      using LambdaVolumePostSkeletonValue = std::integral_constant<bool, T::doLambdaVolumePostSkeleton>;
      template<typename T>
      using LambdaSkeletonValue = std::integral_constant<bool, T::doLambdaSkeleton>;
      template<typename T>
      using LambdaBoundaryValue = std::integral_constant<bool, T::doLambdaBoundary>;

      template<typename T>
      using OneSidedSkeletonRequiredValue = std::integral_constant
      < bool, ( ( T::doAlphaSkeleton || T::doLambdaSkeleton) && ! T::doSkeletonTwoSided)>;
      template<typename T>
      using TwoSidedSkeletonRequiredValue = std::integral_constant
      < bool, ( ( T::doAlphaSkeleton || T::doLambdaSkeleton) && T::doSkeletonTwoSided)>;

      template<typename T>
      using IsLinearValue = std::integral_constant<bool, T::isLinear>;

    public:
      //! \brief Whether to assemble the pattern on the elements, i.e. whether
      //!        or not pattern_volume() should be called.
      enum { doPatternVolume             =
             combineOr<Args,PatternVolumeValue>()             };

      //! \brief Whether to assemble the pattern on the elements after the
      //!        skeleton has been handled, i.e. whether or not
      //!        pattern_volume_post_skeleton() should be called.
      enum { doPatternVolumePostSkeleton =
             combineOr<Args,PatternVolumePostSkeletonValue>() };
      //! \brief Whether to assemble the pattern on the interior
      //!        intersections, i.e. whether or not pattern_skeleton() should
      //!        be called.
      enum { doPatternSkeleton           =
             combineOr<Args,PatternSkeletonValue>()           };
      //! \brief Whether to assemble the pattern on the boundary
      //!        intersections, i.e. whether or not pattern_boundary() should
      //!        be called.
      enum { doPatternBoundary           =
             combineOr<Args,PatternBoundaryValue>()           };

      //! \brief Whether to call the local operator's alpha_volume(),
      //!        jacobian_apply_volume() and jacobian_volume().
      enum { doAlphaVolume               =
             combineOr<Args,AlphaVolumeValue>()               };
      //! \brief Whether to call the local operator's
      //!        alpha_volume_post_skeleton(),
      //!        jacobian_apply_volume_post_skeleton() and
      //!        jacobian_volume_post_skeleton().
      enum { doAlphaVolumePostSkeleton   =
             combineOr<Args,AlphaVolumePostSkeletonValue>()   };
      //! \brief Whether to call the local operator's alpha_skeleton(),
      //!        jacobian_apply_skeleton() and jacobian_skeleton().
      enum { doAlphaSkeleton             =
             combineOr<Args,AlphaSkeletonValue>()             };
      //! \brief Whether to call the local operator's alpha_boundary(),
      //!        jacobian_apply_boundary() and jacobian_boundary().
      enum { doAlphaBoundary             =
             combineOr<Args,AlphaBoundaryValue>()             };

      //! \brief Whether to call the local operator's lambda_volume().
      enum { doLambdaVolume              =
             combineOr<Args,LambdaVolumeValue>()              };
      //! \brief Whether to call the local operator's
      //!        lambda_volume_post_skeleton().
      enum { doLambdaVolumePostSkeleton  =
             combineOr<Args,LambdaVolumePostSkeletonValue>()  };
      //! \brief Whether to call the local operator's lambda_skeleton().
      enum { doLambdaSkeleton            =
             combineOr<Args,LambdaSkeletonValue>()            };
      //! \brief Whether to call the local operator's lambda_boundary().
      enum { doLambdaBoundary            =
             combineOr<Args,LambdaBoundaryValue>()            };

      //! \brief Whether to visit the skeleton methods from both sides
      enum { doSkeletonTwoSided          =
             combineOr<Args,TwoSidedSkeletonRequiredValue>()  };
      static_assert(!(combineOr<Args,OneSidedSkeletonRequiredValue>() &&
                      combineOr<Args,TwoSidedSkeletonRequiredValue>()),
                    "Some summands require a one-sided skelton, others a "
                    "two-sided skeleton.  This is not supported.");

      //! \brief Whether this is a linear operator
      enum { isLinear = combineAnd<Args, IsLinearValue>() };

      //! \} Control flags

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the sparsity pattern
      //! \{
      //

      //! get an element's contribution to the sparsity pattern
      /**
       * \note Summands with zero weight don't contribute to the sparsity
       *       pattern, and the calls to the pattern methods are eliminated at
       *       run-time.
       */
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_volume
      ( const LFSU& lfsu, const LFSV& lfsv,
        LocalPattern& pattern) const
      {
        applyLops(LocalOperatorApply::patternVolume, lfsu, lfsv, pattern);
      }

      //! \brief get an element's contribution to the sparsity pattern after
      //!        the intersections have been handled
      /**
       * \note Summands with zero weight don't contribute to the sparsity
       *       pattern, and the calls to the pattern methods are eliminated at
       *       run-time.
       */
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_volume_post_skeleton
      ( const LFSU& lfsu, const LFSV& lfsv,
        LocalPattern& pattern) const
      {
        applyLops(LocalOperatorApply::patternVolumePostSkeleton, lfsu, lfsv, pattern);
      }

      //! get an internal intersection's contribution to the sparsity pattern
      /**
       * \note Summands with zero weight don't contribute to the sparsity
       *       pattern, and the calls to the pattern methods are eliminated at
       *       run-time.
       */
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_skeleton
      ( const LFSU& lfsu_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const LFSV& lfsv_n,
        LocalPattern& pattern_sn,
        LocalPattern& pattern_ns) const
      {
        applyLops(LocalOperatorApply::patternSkeleton,
          lfsu_s, lfsv_s, lfsu_n, lfsv_n,
                pattern_sn, pattern_ns);
      }

      //! get a boundary intersection's contribution to the sparsity pattern
      /**
       * \note Summands with zero weight don't contribute to the sparsity
       *       pattern, and the calls to the pattern methods are eliminated at
       *       run-time.
       */
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_boundary
      ( const LFSU& lfsu_s, const LFSV& lfsv_s,
        LocalPattern& pattern_ss) const
      {
        applyLops(LocalOperatorApply::patternBoundary, lfsu_s, lfsv_s, pattern_ss);
      }

      //! \} Methods for the sparsity pattern

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the residual -- non-constant parts
      //! \{
      //

      //! get an element's contribution to alpha
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename R>
      void alpha_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        R& r) const
      {
        applyLops(LocalOperatorApply::alphaVolume, eg, lfsu, x, lfsv, r);
      }

      //! \brief get an element's contribution to alpha after the
      //!        intersections have been handled
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename R>
      void alpha_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        R& r) const
      {
        applyLops(LocalOperatorApply::alphaVolumePostSkeleton, eg, lfsu, x, lfsv, r);
      }

      //! get an internal intersections's contribution to alpha
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename R>
      void alpha_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        R& r_s, R& r_n) const
      {
        applyLops(LocalOperatorApply::alphaSkeleton, ig,
                lfsu_s, x_s, lfsv_s,
                lfsu_n, x_n, lfsv_n,
                r_s, r_n);
      }

      //! get a boundary intersections's contribution to alpha
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename R>
      void alpha_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        R& r_s) const
      {
        applyLops(LocalOperatorApply::alphaVolumePostSkeleton, ig, lfsu_s, x_s, lfsv_s, r_s);
      }

      //! \} Methods for the residual -- non-constant parts

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the residual -- constant parts
      //! \{
      //

      //! get an element's contribution to lambda
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSV, typename R>
      void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
      {
        applyLops(LocalOperatorApply::lambdaVolume, eg, lfsv, r);
      }

      //! \brief get an element's contribution to lambda after the
      //!        intersections have been handled
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSV, typename R>
      void lambda_volume_post_skeleton(const EG& eg,
                                       const LFSV& lfsv,
                                       R& r) const
      {
        applyLops(LocalOperatorApply::lambdaVolumePostSkeleton, eg, lfsv, r);
      }

      //! get an internal intersections's contribution to lambda
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSV, typename R>
      void lambda_skeleton(const IG& ig,
                           const LFSV& lfsv_s, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        applyLops(LocalOperatorApply::lambdaSkeleton, ig, lfsv_s, lfsv_n, r_s, r_n);
      }

      //! get a boundary intersections's contribution to lambda
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary(const IG& ig, const LFSV& lfsv_s, R& r_s) const
      {
        applyLops(LocalOperatorApply::lambdaBoundary, ig, lfsv_s, r_s);
      }

      //! \} Methods for the residual -- constant parts

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the application of the jacobian
      //! \{
      //

      //! apply an element's jacobian
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        Y& y) const
      {
        applyLops(LocalOperatorApply::jacobianApplyVolume, eg, lfsu, x, lfsv, y);
      }

      //! \brief apply an element's jacobian after the intersections have been
      //!        handled
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        Y& y) const
      {
        applyLops(LocalOperatorApply::jacobianApplyVolumePostSkeleton, eg, lfsu, x, lfsv, y);
      }

      //! apply an internal intersections's jacobians
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n) const
      {
        applyLops(LocalOperatorApply::jacobianApplySkeleton, ig,
                lfsu_s, x_s, lfsv_s,
                lfsu_n, x_n, lfsv_n,
                y_s, y_n);
      }

      //! apply a boundary intersections's jacobian
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        Y& y_s) const
      {
        applyLops(LocalOperatorApply::jacobianApplyBoundary, ig, lfsu_s, x_s, lfsv_s, y_s);
      }

      //! \} Methods for the application of the jacobian

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods to extract the jacobian
      //! \{
      //

      //! get an element's jacobian
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        LocalMatrix& mat) const
      {
        applyLops(LocalOperatorApply::jacobianVolume, eg, lfsu, x, lfsv, mat);
      }

      //! get an element's jacobian after the intersections have been handled
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        LocalMatrix& mat) const
      {
        applyLops(LocalOperatorApply::jacobianVolumePostSkeleton, eg, lfsu, x, lfsv, mat);
      }

      //! apply an internal intersections's jacobians
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        LocalMatrix& mat_ss, LocalMatrix& mat_sn,
        LocalMatrix& mat_ns, LocalMatrix& mat_nn) const
      {
        applyLops(LocalOperatorApply::jacobianSkeleton, ig,
                lfsu_s, x_s, lfsv_s,
                lfsu_n, x_n, lfsv_n,
                mat_ss, mat_sn, mat_ns, mat_nn);
      }

      //! get a boundary intersections's jacobian
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        LocalMatrix& mat_ss) const
      {
        applyLops(LocalOperatorApply::jacobianBoundary, ig, lfsu_s, x_s, lfsv_s, mat_ss);
      }

      //! \} Methods to extract the jacobian

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for instationary problems
      //! \{
      //

      //! Export type used for time values
      typedef typename std::tuple_element<0, Args>::type::RealType RealType;

    private:
      // template meta program helpers for the methods related to instationary
      // stuff

      struct Apply {
      template<typename LOP>
      static void setTime(const LOP& lop, RealType t)
      {
        lop.setTime(t);
      }

      template<typename LOP>
      static void preStep(const LOP& lop,
                          RealType time, RealType dt, int stages)
      {
        lop.preStep(time, dt, stages);
      }

      template<typename LOP>
      static void postStep(const LOP& lop)
      {
        lop.postStep();
      }

      template<typename LOP>
      static void preStage(const LOP& lop, RealType time, int r)
      {
        lop.preStage(time, r);
      }

      template<typename LOP>
      static void postStage(const LOP& lop)
      {
        lop.postStage();
      }

      template<typename LOP>
      static RealType suggestTimestep(const LOP& lop, RealType & dt)
      {
        dt = std::min(dt,lop.suggestTimestep(dt));
      }
      };

    public:
      //! set time for subsequent evaluation
      void setTime (RealType t)
      {
        applyLops(Apply::setTime, t);
      }

      //! get current time
      RealType getTime () const
      {
        return get<0>(lops)->getTime();
      }

      //! to be called once before each time step
      void preStep (RealType time, RealType dt, int stages)
      {
        applyLops(Apply::preStep, time, dt, stages);
      }

      //! to be called once at the end of each time step
      void postStep ()
      {
        applyLops(Apply::postStep, lops);
      }

      //! to be called once before each stage
      void preStage (RealType time, int r)
      {
        applyLops(Apply::preStage, time, r);
      }

      //! get current stage
      int getStage () const
      {
        return get<0>(lops)->getStage();
      }

      //! to be called once at the end of each stage
      void postStage ()
      {
        applyLops(Apply::postStage, lops);
      }

      //! to be called after stage 1
      /**
       * \note This operator simply chains suggestTimestep() methods of all
       *       the component local operators together and hopes that the
       *       result will be meaningful.
       */
      RealType suggestTimestep (RealType dt) const
      {
        applyLops(Apply::suggestTimestep, dt);
        return dt;
      }

      //! \} Methods for instationary problems
    };

  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_WEIGHTEDSUM_HH
