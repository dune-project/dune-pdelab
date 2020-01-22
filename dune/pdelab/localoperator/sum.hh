// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_SUM_HH
#define DUNE_PDELAB_LOCALOPERATOR_SUM_HH

#include <cstddef>

#include <tuple>
#include <dune/pdelab/common/forloop.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/typetraits.hh>

#include <dune/pdelab/localoperator/callswitch.hh>

namespace Dune {
  namespace PDELab {

template<template<typename> typename F, typename... Args>
constexpr auto staticCombineOr(const std::tuple<Args...> & t)
{
    return std::disjunction<F<Args>...>{};
}

template<typename Tuple, template<typename> typename F>
constexpr auto combineOr() {
    return decltype(staticCombineOr<F>(std::declval<Tuple>())){};
}

    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    //! A local operator to take the sum of other local operators
    /**
     * \nosubgrouping
     *
     * \tparam Args Tuple of local operators.  Must fulfill \c
     *              tuple_size<Args>::value>=1.
     */
    template<typename Args>
    class InstationarySumLocalOperator
    {
      static const std::size_t size = std::tuple_size<Args>::value;

      typedef typename ForEachType<AddPtrTypeEvaluator, Args>::Type ArgPtrs;
      typedef typename ForEachType<AddRefTypeEvaluator, Args>::Type ArgRefs;

      ArgPtrs lops;

    public:
      //////////////////////////////////////////////////////////////////////
      //
      //! \name Construction and modification
      //! \{
      //

      //! \brief construct a InstationarySumLocalOperator from a tuple of
      //!        local operators
      InstationarySumLocalOperator(const ArgRefs& lops_)
        : lops(transformTuple<AddPtrTypeEvaluator>(lops_))
      { }

      //! set the i'th component of the sum
      template<std::size_t i>
      void setSummand(typename std::tuple_element<i,Args>::type& summand)
      { std::get<i>(lops) = &summand; }

      //! get the i'th component of the sum
      template<std::size_t i>
      typename std::tuple_element<i,Args>::type& getSummand()
      { return *std::get<i>(lops); }

      //! \} Construction and modification

      ////////////////////////////////////////////////////////////////////////
      //
      //! \name Control flags
      //! \{
      //

    private:
      template<typename T>
      using PatternVolumeValue = std::integral_constant<bool, T::PatternVolumeValue>;
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

    public:
      //! \brief Whether to assemble the pattern on the elements, i.e. whether
      //!        or not pattern_volume() should be called.
      enum { doPatternVolume             =
             combineOr<Args,PatternVolumeValue>             };
      //! \brief Whether to assemble the pattern on the elements after the
      //!        skeleton has been handled, i.e. whether or not
      //!        pattern_volume_post_skeleton() should be called.
      enum { doPatternVolumePostSkeleton =
             combineOr<Args,PatternVolumePostSkeletonValue> };
      //! \brief Whether to assemble the pattern on the interior
      //!        intersections, i.e. whether or not pattern_skeleton() should
      //!        be called.
      enum { doPatternSkeleton           =
             combineOr<Args,PatternSkeletonValue>           };
      //! \brief Whether to assemble the pattern on the boundary
      //!        intersections, i.e. whether or not pattern_boundary() should
      //!        be called.
      enum { doPatternBoundary           =
             combineOr<Args,PatternBoundaryValue>           };

      //! \brief Whether to call the local operator's alpha_volume(),
      //!        jacobian_apply_volume() and jacobian_volume().
      enum { doAlphaVolume               =
             combineOr<Args,AlphaVolumeValue>               };
      //! \brief Whether to call the local operator's
      //!        alpha_volume_post_skeleton(),
      //!        jacobian_apply_volume_post_skeleton() and
      //!        jacobian_volume_post_skeleton().
      enum { doAlphaVolumePostSkeleton   =
             combineOr<Args,AlphaVolumePostSkeletonValue>   };
      //! \brief Whether to call the local operator's alpha_skeleton(),
      //!        jacobian_apply_skeleton() and jacobian_skeleton().
      enum { doAlphaSkeleton             =
             combineOr<Args,AlphaSkeletonValue>             };
      //! \brief Whether to call the local operator's alpha_boundary(),
      //!        jacobian_apply_boundary() and jacobian_boundary().
      enum { doAlphaBoundary             =
             combineOr<Args,AlphaBoundaryValue>             };

      //! \brief Whether to call the local operator's lambda_volume().
      enum { doLambdaVolume              =
             combineOr<Args,LambdaVolumeValue>              };
      //! \brief Whether to call the local operator's
      //!        lambda_volume_post_skeleton().
      enum { doLambdaVolumePostSkeleton  =
             combineOr<Args,LambdaVolumePostSkeletonValue>  };
      //! \brief Whether to call the local operator's lambda_skeleton().
      enum { doLambdaSkeleton            =
             combineOr<Args,LambdaSkeletonValue>            };
      //! \brief Whether to call the local operator's lambda_boundary().
      enum { doLambdaBoundary            =
             combineOr<Args,LambdaBoundaryValue>            };

      //! \brief Whether to visit the skeleton methods from both sides
      enum { doSkeletonTwoSided          =
             combineOr<Args,TwoSidedSkeletonRequiredValue>  };
      static_assert(!(combineOr<Args,OneSidedSkeletonRequiredValue> &&
                      combineOr<Args,TwoSidedSkeletonRequiredValue>),
                    "Some summands require a one-sided skelton, others a "
                    "two-sided skeleton.  This is not supported.");

      //! \} Control flags

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the sparsity pattern
      //! \{
      //

    private:
      // template meta program helpers for the pattern_* methods

      template<int i>
      struct PatternVolumeOperation {
        template<typename LFSU, typename LFSV, typename LocalPattern>
        static void apply(const ArgPtrs& lops,
                          const LFSU& lfsu, const LFSV& lfsv,
                          LocalPattern& pattern)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doPatternVolume>::
            pattern_volume(*get<i>(lops), lfsu, lfsv, pattern);
        }
      };

      template<int i>
      struct PatternVolumePostSkeletonOperation {
        template<typename LFSU, typename LFSV, typename LocalPattern>
        static void apply(const ArgPtrs& lops,
                          const LFSU& lfsu, const LFSV& lfsv,
                          LocalPattern& pattern)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doPatternVolumePostSkeleton>::
            pattern_volume_post_skeleton(*get<i>(lops), lfsu, lfsv, pattern);
        }
      };

      template<int i>
      struct PatternSkeletonOperation {
        template<typename LFSU, typename LFSV, typename LocalPattern>
        static void apply(const ArgPtrs& lops,
                          const LFSU& lfsu_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const LFSV& lfsv_n,
                          LocalPattern& pattern_sn,
                          LocalPattern& pattern_ns)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doPatternSkeleton>::
            pattern_skeleton(*get<i>(lops),
                             lfsu_s, lfsv_s, lfsu_n, lfsv_n,
                             pattern_sn, pattern_ns);
        }
      };

      template<int i>
      struct PatternBoundaryOperation {
        template<typename LFSU, typename LFSV, typename LocalPattern>
        static void apply(const ArgPtrs& lops,
                          const LFSU& lfsu_s, const LFSV& lfsv_s,
                          LocalPattern& pattern_ss)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doPatternBoundary>::
            pattern_boundary(*get<i>(lops), lfsu_s, lfsv_s, pattern_ss);
        }
      };

    public:
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
        ForLoop<PatternVolumeOperation, 0, size-1>::
          apply(lops, lfsu, lfsv, pattern);
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
        ForLoop<PatternVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, lfsu, lfsv, pattern);
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
        ForLoop<PatternSkeletonOperation, 0, size-1>::
          apply(lops, lfsu_s, lfsv_s, lfsu_n, lfsv_n,
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
        ForLoop<PatternBoundaryOperation, 0, size-1>::
          apply(lops, lfsu_s, lfsv_s, pattern_ss);
      }

      //! \} Methods for the sparsity pattern

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the residual -- non-constant parts
      //! \{
      //

    private:
      // template meta program helpers for the alpha_* methods

      template<int i>
      struct AlphaVolumeOperation {
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename R>
        static void apply(const ArgPtrs& lops, const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          R& r)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doAlphaVolume>::
          alpha_volume(*get<i>(lops), eg, lfsu, x, lfsv, r);
        }
      };

      template<int i>
      struct AlphaVolumePostSkeletonOperation {
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename R>
        static void apply(const ArgPtrs& lops, const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          R& r)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doAlphaVolumePostSkeleton>::
            alpha_volume_post_skeleton(*get<i>(lops), eg,
                                       lfsu, x, lfsv,
                                       r);
        }
      };

      template<int i>
      struct AlphaSkeletonOperation {
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename R>
        static void apply(const ArgPtrs& lops, const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          R& r_s, R& r_n)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doAlphaSkeleton>::
            alpha_skeleton(*get<i>(lops), ig,
                           lfsu_s, x_s, lfsv_s,
                           lfsu_n, x_n, lfsv_n,
                           r_s, r_n);
        }
      };

      template<int i>
      struct AlphaBoundaryOperation {
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename R>
        static void apply(const ArgPtrs& lops, const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          R& r_s)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doAlphaBoundary>::
            alpha_boundary(*get<i>(lops), ig, lfsu_s, x_s, lfsv_s, r_s);
        }
      };

    public:
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
        ForLoop<AlphaVolumeOperation, 0, size-1>::
          apply(lops, eg, lfsu, x, lfsv, r);
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
        ForLoop<AlphaVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, eg, lfsu, x, lfsv, r);
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
        ForLoop<AlphaSkeletonOperation, 0, size-1>::
          apply(lops, ig,
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
        ForLoop<AlphaVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, ig, lfsu_s, x_s, lfsv_s, r_s);
      }

      //! \} Methods for the residual -- non-constant parts

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the residual -- constant parts
      //! \{
      //

    private:
      // template meta program helpers for the lambda_* methods

      template<int i>
      struct LambdaVolumeOperation {
        template<typename EG, typename LFSV, typename R>
        static void apply(const ArgPtrs& lops, const EG& eg,
                          const LFSV& lfsv,
                          R& r)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doLambdaVolume>::
            lambda_volume(*get<i>(lops), eg, lfsv, r);
        }
      };

      template<int i>
      struct LambdaVolumePostSkeletonOperation {
        template<typename EG, typename LFSV, typename R>
        static void apply(const ArgPtrs& lops, const EG& eg,
                          const LFSV& lfsv,
                          R& r)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doLambdaVolumePostSkeleton>::
            lambda_volume_post_skeleton(*get<i>(lops), eg, lfsv, r);
        }
      };

      template<int i>
      struct LambdaSkeletonOperation {
        template<typename IG, typename LFSV, typename R>
        static void apply(const ArgPtrs& lops, const IG& ig,
                          const LFSV& lfsv_s,
                          const LFSV& lfsv_n,
                          R& r_s, R& r_n)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doLambdaSkeleton>::
            lambda_skeleton(*get<i>(lops), ig,
                            lfsv_s, lfsv_n,
                            r_s, r_n);
        }
      };

      template<int i>
      struct LambdaBoundaryOperation {
        template<typename IG, typename LFSV, typename R>
        static void apply(const ArgPtrs& lops, const IG& ig,
                          const LFSV& lfsv_s,
                          R& r_s)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doLambdaBoundary>::
            lambda_boundary(*get<i>(lops), ig, lfsv_s, r_s);
        }
      };

    public:
      //! get an element's contribution to lambda
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSV, typename R>
      void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
      {
        ForLoop<LambdaVolumeOperation, 0, size-1>::
          apply(lops, eg, lfsv, r);
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
        ForLoop<LambdaVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, eg, lfsv, r);
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
        ForLoop<LambdaSkeletonOperation, 0, size-1>::
          apply(lops, ig, lfsv_s, lfsv_n, r_s, r_n);
      }

      //! get a boundary intersections's contribution to lambda
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary(const IG& ig, const LFSV& lfsv_s, R& r_s) const
      {
        ForLoop<LambdaBoundaryOperation, 0, size-1>::
          apply(lops, ig, lfsv_s, r_s);
      }

      //! \} Methods for the residual -- constant parts

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the application of the jacobian
      //! \{
      //

    private:
      // template meta program helpers for the jacobian_apply_* methods

      template<int i>
      struct JacobianApplyVolumeOperation {
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename Y>
        static void apply(const ArgPtrs& lops, const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          Y& y)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doAlphaVolume>::
            jacobian_apply_volume(*get<i>(lops), eg, lfsu, x, lfsv, y);
        }
      };

      template<int i>
      struct JacobianApplyVolumePostSkeletonOperation {
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename Y>
        static void apply(const ArgPtrs& lops, const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          Y& y)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doAlphaVolumePostSkeleton>::
            jacobian_apply_volume_post_skeleton(*get<i>(lops), eg,
                                                lfsu, x, lfsv,
                                                y);
        }
      };

      template<int i>
      struct JacobianApplySkeletonOperation {
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename Y>
        static void apply(const ArgPtrs& lops, const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          Y& y_s, Y& y_n)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doAlphaSkeleton>::
          jacobian_apply_skeleton(*get<i>(lops), ig,
                                  lfsu_s, x_s, lfsv_s,
                                  lfsu_n, x_n, lfsv_n,
                                  y_s, y_n);
        }
      };

      template<int i>
      struct JacobianApplyBoundaryOperation {
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename Y>
        static void apply(const ArgPtrs& lops, const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          Y& y_s)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doAlphaBoundary>::
            jacobian_apply_boundary(*get<i>(lops), ig,
                                    lfsu_s, x_s, lfsv_s,
                                    y_s);
        }
      };

    public:
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
        ForLoop<JacobianApplyVolumeOperation, 0, size-1>::
          apply(lops, eg, lfsu, x, lfsv, y);
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
        ForLoop<JacobianApplyVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, eg, lfsu, x, lfsv, y);
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
        ForLoop<JacobianApplySkeletonOperation, 0, size-1>::
          apply(lops, ig,
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
        ForLoop<JacobianApplyBoundaryOperation, 0, size-1>::
          apply(lops, ig, lfsu_s, x_s, lfsv_s, y_s);
      }

      //! \} Methods for the application of the jacobian

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods to extract the jacobian
      //! \{
      //

    private:
      // template meta program helpers for the jacobian_apply_* methods

      template<int i>
      struct JacobianVolumeOperation {
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename LocalMatrix>
        static void apply(const ArgPtrs& lops, const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          LocalMatrix& mat)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doAlphaVolume>::
            jacobian_volume(*get<i>(lops), eg, lfsu, x, lfsv, mat);
        }
      };

      template<int i>
      struct JacobianVolumePostSkeletonOperation {
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename LocalMatrix>
        static void apply(const ArgPtrs& lops, const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          LocalMatrix& mat)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doAlphaVolumePostSkeleton>::
            jacobian_volume_post_skeleton(*get<i>(lops), eg,
                                          lfsu, x, lfsv,
                                          mat);
        }
      };

      template<int i>
      struct JacobianSkeletonOperation {
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename LocalMatrix>
        static void apply(const ArgPtrs& lops, const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          LocalMatrix& mat_ss, LocalMatrix& mat_sn,
                          LocalMatrix& mat_ns, LocalMatrix& mat_nn)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doAlphaSkeleton>::
            jacobian_skeleton(*get<i>(lops), ig,
                              lfsu_s, x_s, lfsv_s,
                              lfsu_n, x_n, lfsv_n,
                              mat_ss, mat_sn, mat_ns, mat_nn);
        }
      };

      template<int i>
      struct JacobianBoundaryOperation {
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename LocalMatrix>
        static void apply(const ArgPtrs& lops, const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          LocalMatrix& mat_ss)
        {
          LocalAssemblerCallSwitch<typename std::tuple_element<i,Args>::type,
            std::tuple_element<i,Args>::type::doAlphaBoundary>::
            jacobian_boundary(*get<i>(lops), ig,
                              lfsu_s, x_s, lfsv_s, mat_ss);
        }
      };

    public:
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
        ForLoop<JacobianVolumeOperation, 0, size-1>::
          apply(lops, eg, lfsu, x, lfsv, mat);
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
        ForLoop<JacobianVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, eg, lfsu, x, lfsv, mat);
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
        ForLoop<JacobianSkeletonOperation, 0, size-1>::
          apply(lops, ig,
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
        ForLoop<JacobianBoundaryOperation, 0, size-1>::
          apply(lops, ig, lfsu_s, x_s, lfsv_s, mat_ss);
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

      template<int i> struct SetTimeOperation {
        static void apply(ArgPtrs& lops, RealType t)
        { get<i>(lops)->setTime(t); }
      };

      template<int i> struct PreStepOperation {
        static void apply(ArgPtrs& lops,
                          RealType time, RealType dt, int stages)
        { get<i>(lops)->preStep(time, dt, stages); }
      };

      template<int i> struct PostStepOperation {
        static void apply(ArgPtrs& lops)
        { get<i>(lops)->postStep(); }
      };

      template<int i> struct PreStageOperation {
        static void apply(ArgPtrs& lops, RealType time, int r)
        { get<i>(lops)->preStage(time, r); }
      };

      template<int i> struct PostStageOperation {
        static void apply(ArgPtrs& lops)
        { get<i>(lops)->postStage(); }
      };

      template<int i> struct SuggestTimestepOperation {
        static void apply(ArgPtrs& lops, RealType& dt)
        { dt = get<i>(lops)->suggestTimestep(dt); }
      };

    public:
      //! set time for subsequent evaluation
      void setTime (RealType t)
      {
        ForLoop<SetTimeOperation, 0, size-1>::apply(lops, t);
      }

      //! get current time
      RealType getTime () const
      {
        return get<0>(lops)->getTime();
      }

      //! to be called once before each time step
      void preStep (RealType time, RealType dt, int stages)
      {
        ForLoop<PreStepOperation, 0, size-1>::apply(lops, time, dt, stages);
      }

      //! to be called once at the end of each time step
      void postStep ()
      {
        ForLoop<PostStepOperation, 0, size-1>::apply(lops);
      }

      //! to be called once before each stage
      void preStage (RealType time, int r)
      {
        ForLoop<PreStageOperation, 0, size-1>::apply(lops, time, r);
      }

      //! get current stage
      int getStage () const
      {
        return get<0>(lops)->getStage();
      }

      //! to be called once at the end of each stage
      void postStage ()
      {
        ForLoop<PostStageOperation, 0, size-1>::apply(lops);
      }

      //! to be called after stage 1
      /**
       * \note This operator simply chains suggestTimestep() methods of all
       *       the component local operators together and hopes that the
       *       result will be meaningful.
       */
      RealType suggestTimestep (RealType dt) const
      {
        ForLoop<SuggestTimestepOperation, 0, size-1>::apply(lops, dt);
        return dt;
      }

      //! \} Methods for instationary problems
    };

    //! \} group LocalOperatorDefaultImp
  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_SUM_HH
