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

#include <dune/pdelab/localoperator/callswitch.hh>

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
      static const std::size_t size = tuple_size<Args>::value;

      typedef typename ForEachType<AddPtrTypeEvaluator, Args>::Type ArgPtrs;
      typedef typename ForEachType<AddRefTypeEvaluator, Args>::Type ArgRefs;

      ArgPtrs lops;
      typedef FieldVector<K, size> Weights;
      Weights weights;

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
      void setSummand(typename tuple_element<i,Args>::type& summand)
      { get<i>(lops) = &summand; }

      //! get the i'th component of the sum
      template<std::size_t i>
      typename tuple_element<i,Args>::type& getSummand()
      { return *get<i>(lops); }

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
      template<typename T1, typename T2>
      struct OrOperation
        : public std::integral_constant<bool, T1::value || T2:: value>
      { };
      template<template<int> class Value>
      struct AccFlag : public GenericForLoop<OrOperation, Value, 0, size-1>
      { };

      template<int i>
      struct PatternVolumeValue : public std::integral_constant
      < bool, tuple_element<i, Args>::type::doPatternVolume>
      { };
      template<int i>
      struct PatternVolumePostSkeletonValue : public std::integral_constant
      < bool, tuple_element<i, Args>::type::doPatternVolumePostSkeleton>
      { };
      template<int i>
      struct PatternSkeletonValue : public std::integral_constant
      < bool, tuple_element<i, Args>::type::doPatternSkeleton>
      { };
      template<int i>
      struct PatternBoundaryValue : public std::integral_constant
      < bool, tuple_element<i, Args>::type::doPatternBoundary>
      { };

      template<int i>
      struct AlphaVolumeValue : public std::integral_constant
      < bool, tuple_element<i, Args>::type::doAlphaVolume>
      { };
      template<int i>
      struct AlphaVolumePostSkeletonValue : public std::integral_constant
      < bool, tuple_element<i, Args>::type::doAlphaVolumePostSkeleton>
      { };
      template<int i>
      struct AlphaSkeletonValue : public std::integral_constant
      < bool, tuple_element<i, Args>::type::doAlphaSkeleton>
      { };
      template<int i>
      struct AlphaBoundaryValue : public std::integral_constant
      < bool, tuple_element<i, Args>::type::doAlphaBoundary>
      { };

      template<int i>
      struct LambdaVolumeValue : public std::integral_constant
      < bool, tuple_element<i, Args>::type::doLambdaVolume>
      { };
      template<int i>
      struct LambdaVolumePostSkeletonValue : public std::integral_constant
      < bool, tuple_element<i, Args>::type::doLambdaVolumePostSkeleton>
      { };
      template<int i>
      struct LambdaSkeletonValue : public std::integral_constant
      < bool, tuple_element<i, Args>::type::doLambdaSkeleton>
      { };
      template<int i>
      struct LambdaBoundaryValue : public std::integral_constant
      < bool, tuple_element<i, Args>::type::doLambdaBoundary>
      { };

      template<int i>
      struct OneSidedSkeletonRequiredValue : public std::integral_constant
      < bool, ( ( tuple_element<i, Args>::type::doAlphaSkeleton ||
                  tuple_element<i, Args>::type::doLambdaSkeleton) &&
                ! tuple_element<i, Args>::type::doSkeletonTwoSided)>
      { };
      template<int i>
      struct TwoSidedSkeletonRequiredValue : public std::integral_constant
      < bool, ( ( tuple_element<i, Args>::type::doAlphaSkeleton ||
                  tuple_element<i, Args>::type::doLambdaSkeleton) &&
                tuple_element<i, Args>::type::doSkeletonTwoSided)>
      { };

    public:
      //! \brief Whether to assemble the pattern on the elements, i.e. whether
      //!        or not pattern_volume() should be called.
      enum { doPatternVolume             =
             AccFlag<PatternVolumeValue>::value             };
      //! \brief Whether to assemble the pattern on the elements after the
      //!        skeleton has been handled, i.e. whether or not
      //!        pattern_volume_post_skeleton() should be called.
      enum { doPatternVolumePostSkeleton =
             AccFlag<PatternVolumePostSkeletonValue>::value };
      //! \brief Whether to assemble the pattern on the interior
      //!        intersections, i.e. whether or not pattern_skeleton() should
      //!        be called.
      enum { doPatternSkeleton           =
             AccFlag<PatternSkeletonValue>::value           };
      //! \brief Whether to assemble the pattern on the boundary
      //!        intersections, i.e. whether or not pattern_boundary() should
      //!        be called.
      enum { doPatternBoundary           =
             AccFlag<PatternBoundaryValue>::value           };

      //! \brief Whether to call the local operator's alpha_volume(),
      //!        jacobian_apply_volume() and jacobian_volume().
      enum { doAlphaVolume               =
             AccFlag<AlphaVolumeValue>::value               };
      //! \brief Whether to call the local operator's
      //!        alpha_volume_post_skeleton(),
      //!        jacobian_apply_volume_post_skeleton() and
      //!        jacobian_volume_post_skeleton().
      enum { doAlphaVolumePostSkeleton   =
             AccFlag<AlphaVolumePostSkeletonValue>::value   };
      //! \brief Whether to call the local operator's alpha_skeleton(),
      //!        jacobian_apply_skeleton() and jacobian_skeleton().
      enum { doAlphaSkeleton             =
             AccFlag<AlphaSkeletonValue>::value             };
      //! \brief Whether to call the local operator's alpha_boundary(),
      //!        jacobian_apply_boundary() and jacobian_boundary().
      enum { doAlphaBoundary             =
             AccFlag<AlphaBoundaryValue>::value             };

      //! \brief Whether to call the local operator's lambda_volume().
      enum { doLambdaVolume              =
             AccFlag<LambdaVolumeValue>::value              };
      //! \brief Whether to call the local operator's
      //!        lambda_volume_post_skeleton().
      enum { doLambdaVolumePostSkeleton  =
             AccFlag<LambdaVolumePostSkeletonValue>::value  };
      //! \brief Whether to call the local operator's lambda_skeleton().
      enum { doLambdaSkeleton            =
             AccFlag<LambdaSkeletonValue>::value            };
      //! \brief Whether to call the local operator's lambda_boundary().
      enum { doLambdaBoundary            =
             AccFlag<LambdaBoundaryValue>::value            };

      //! \brief Whether to visit the skeleton methods from both sides
      enum { doSkeletonTwoSided          =
             AccFlag<TwoSidedSkeletonRequiredValue>::value  };
      static_assert(!(AccFlag<OneSidedSkeletonRequiredValue>::value &&
                      AccFlag<TwoSidedSkeletonRequiredValue>::value),
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
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const LFSU& lfsu, const LFSV& lfsv,
                          LocalPattern& pattern)
        {
          if(weights[i] != K(0))
            LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
              tuple_element<i,Args>::type::doPatternVolume>::
              pattern_volume(*get<i>(lops), lfsu, lfsv, pattern);
        }
      };

      template<int i>
      struct PatternVolumePostSkeletonOperation {
        template<typename LFSU, typename LFSV, typename LocalPattern>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const LFSU& lfsu, const LFSV& lfsv,
                          LocalPattern& pattern)
        {
          if(weights[i] != K(0))
            LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
              tuple_element<i,Args>::type::doPatternVolumePostSkeleton>::
              pattern_volume_post_skeleton(*get<i>(lops), lfsu, lfsv, pattern);
        }
      };

      template<int i>
      struct PatternSkeletonOperation {
        template<typename LFSU, typename LFSV, typename LocalPattern>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const LFSU& lfsu_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const LFSV& lfsv_n,
                          LocalPattern& pattern_sn,
                          LocalPattern& pattern_ns)
        {
          if(weights[i] != K(0))
            LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
              tuple_element<i,Args>::type::doPatternSkeleton>::
              pattern_skeleton(*get<i>(lops),
                               lfsu_s, lfsv_s, lfsu_n, lfsv_n,
                               pattern_sn, pattern_ns);
        }
      };

      template<int i>
      struct PatternBoundaryOperation {
        template<typename LFSU, typename LFSV, typename LocalPattern>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const LFSU& lfsu_s, const LFSV& lfsv_s,
                          LocalPattern& pattern_ss)
        {
          if(weights[i] != K(0))
            LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
              tuple_element<i,Args>::type::doPatternBoundary>::
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
          apply(lops, weights, lfsu, lfsv, pattern);
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
          apply(lops, weights, lfsu, lfsv, pattern);
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
          apply(lops, weights, lfsu_s, lfsv_s, lfsu_n, lfsv_n,
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
          apply(lops, weights, lfsu_s, lfsv_s, pattern_ss);
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
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          WeightedVectorAccumulationView<C>& r)
        {
          apply(lops, weights[i]*r.weight(), eg, lfsu, x, lfsv, r.container());
        }
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, typename C::weight_type weight,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          C& r)
        {
          if(weight != K(0)) {
            WeightedVectorAccumulationView<C> view(r, weight);
            LocalAssemblerCallSwitch<Arg, Arg::doAlphaVolume>::
              alpha_volume(*get<i>(lops), eg, lfsu, x, lfsv, view);
          }
        }
      };

      template<int i>
      struct AlphaVolumePostSkeletonOperation {
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          WeightedVectorAccumulationView<C>& r)
        {
          apply(lops, weights[i]*r.weight(), eg, lfsu, x, lfsv, r.container());
        }
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, typename C::weight_type weight,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          C& r)
        {
          if(weight != K(0)) {
            WeightedVectorAccumulationView<C> view(r, weight);
            LocalAssemblerCallSwitch<Arg, Arg::doAlphaVolumePostSkeleton>::
              alpha_volume_post_skeleton(*get<i>(lops), eg,
                                         lfsu, x, lfsv,
                                         view);
          }
        }
      };

      template<int i>
      struct AlphaSkeletonOperation {
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          WeightedVectorAccumulationView<C>& r_s,
                          WeightedVectorAccumulationView<C>& r_n)
        {
          apply(lops, weights[i]*r_s.weight(), weights[i]*r_n.weight(),
                ig,
                lfsu_s, x_s, lfsv_s,
                lfsu_n, x_n, lfsv_n,
                r_s.container(), r_n.container());
        }
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops,
                          typename C::weight_type weight_s,
                          typename C::weight_type weight_n,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          C& r_s, C& r_n)
        {
          if(weight_s != K(0) || weight_n != K(0)) {
            WeightedVectorAccumulationView<C> view_s(r_s, weight_s);
            WeightedVectorAccumulationView<C> view_n(r_n, weight_n);
            LocalAssemblerCallSwitch<Arg, Arg::doAlphaSkeleton>::
              alpha_skeleton(*get<i>(lops), ig,
                             lfsu_s, x_s, lfsv_s,
                             lfsu_n, x_n, lfsv_n,
                             view_s, view_n);
          }
        }
      };

      template<int i>
      struct AlphaBoundaryOperation {
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          WeightedVectorAccumulationView<C>& r_s)
        {
          apply(lops, weights[i]*r_s.weight(), ig,
                lfsu_s, x_s, lfsv_s,
                r_s.container());
        }
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops,
                          typename C::weight_type weight_s,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          C& r_s)
        {
          if(weight_s != K(0)) {
            WeightedVectorAccumulationView<C> view_s(r_s, weight_s);
            LocalAssemblerCallSwitch<Arg, Arg::doAlphaBoundary>::
              alpha_boundary(*get<i>(lops), ig, lfsu_s, x_s, lfsv_s, view_s);
          }
        }
      };

    public:
      //! get an element's contribution to alpha
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename C>
      void alpha_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        C& r) const
      {
        ForLoop<AlphaVolumeOperation, 0, size-1>::
          apply(lops, weights, eg, lfsu, x, lfsv, r);
      }

      //! \brief get an element's contribution to alpha after the
      //!        intersections have been handled
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename C>
      void alpha_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        C& r) const
      {
        ForLoop<AlphaVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, weights, eg, lfsu, x, lfsv, r);
      }

      //! get an internal intersections's contribution to alpha
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename C>
      void alpha_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        C& r_s, C& r_n) const
      {
        ForLoop<AlphaSkeletonOperation, 0, size-1>::
          apply(lops, weights, ig,
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
               typename C>
      void alpha_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        C& r_s) const
      {
        ForLoop<AlphaBoundaryOperation, 0, size-1>::
          apply(lops, weights, ig, lfsu_s, x_s, lfsv_s, r_s);
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
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename EG, typename LFSV, typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const EG& eg,
                          const LFSV& lfsv,
                          WeightedVectorAccumulationView<C>& r)
        {
          apply(lops, weights[i]*r.weight(), eg, lfsv, r.container());
        }
        template<typename EG, typename LFSV, typename C>
        static void apply(const ArgPtrs& lops, typename C::weight_type weight,
                          const EG& eg,
                          const LFSV& lfsv,
                          C& r)
        {
          if(weight != K(0)) {
            WeightedVectorAccumulationView<C> view(r, weight);
            LocalAssemblerCallSwitch<Arg, Arg::doLambdaVolume>::
              lambda_volume(*get<i>(lops), eg, lfsv, view);
          }
        }
      };

      template<int i>
      struct LambdaVolumePostSkeletonOperation {
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename EG, typename LFSV, typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const EG& eg,
                          const LFSV& lfsv,
                          WeightedVectorAccumulationView<C>& r)
        {
          apply(lops, weights[i]*r.weight(), eg, lfsv, r.container());
        }
        template<typename EG, typename LFSV, typename C>
        static void apply(const ArgPtrs& lops, typename C::weight_type weight,
                          const EG& eg,
                          const LFSV& lfsv,
                          C& r)
        {
          if(weight != K(0)) {
            WeightedVectorAccumulationView<C> view(r, weight);
            LocalAssemblerCallSwitch<Arg, Arg::doLambdaVolumePostSkeleton>::
              lambda_volume_post_skeleton(*get<i>(lops), eg, lfsv, view);
          }
        }
      };

      template<int i>
      struct LambdaSkeletonOperation {
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename IG, typename LFSV, typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const IG& ig,
                          const LFSV& lfsv_s, const LFSV& lfsv_n,
                          WeightedVectorAccumulationView<C>& r_s,
                          WeightedVectorAccumulationView<C>& r_n)
        {
          apply(lops, weights[i]*r_s.weight(), weights[i]*r_n.weight(),
                ig,
                lfsv_s, lfsv_n,
                r_s.container(), r_n.container());
        }
        template<typename IG, typename LFSV, typename C>
        static void apply(const ArgPtrs& lops,
                          typename C::weight_type weight_s,
                          typename C::weight_type weight_n,
                          const IG& ig,
                          const LFSV& lfsv_s, const LFSV& lfsv_n,
                          C& r_s, C& r_n)
        {
          if(weight_s != K(0) || weight_n != K(0)) {
            WeightedVectorAccumulationView<C> view_s(r_s, weight_s);
            WeightedVectorAccumulationView<C> view_n(r_n, weight_n);
            LocalAssemblerCallSwitch<Arg, Arg::doLambdaSkeleton>::
              lambda_skeleton(*get<i>(lops), ig,
                              lfsv_s, lfsv_n,
                              view_s, view_n);
          }
        }
      };

      template<int i>
      struct LambdaBoundaryOperation {
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename IG, typename LFSV, typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const IG& ig,
                          const LFSV& lfsv_s,
                          WeightedVectorAccumulationView<C>& r_s)
        {
          apply(lops, weights[i]*r_s.weight(), ig, lfsv_s, r_s.container());
        }
        template<typename IG, typename LFSV, typename C>
        static void apply(const ArgPtrs& lops,
                          typename C::weight_type weight_s,
                          const IG& ig,
                          const LFSV& lfsv_s,
                          C& r_s)
        {
          if(weight_s != K(0)) {
            WeightedVectorAccumulationView<C> view_s(r_s, weight_s);
            LocalAssemblerCallSwitch<Arg, Arg::doLambdaBoundary>::
              lambda_boundary(*get<i>(lops), ig, lfsv_s, view_s);
          }
        }
      };

    public:
      //! get an element's contribution to lambda
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSV, typename C>
      void lambda_volume(const EG& eg, const LFSV& lfsv, C& r) const
      {
        ForLoop<LambdaVolumeOperation, 0, size-1>::
          apply(lops, weights, eg, lfsv, r);
      }

      //! \brief get an element's contribution to lambda after the
      //!        intersections have been handled
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSV, typename C>
      void lambda_volume_post_skeleton(const EG& eg,
                                       const LFSV& lfsv,
                                       C& r) const
      {
        ForLoop<LambdaVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, weights, eg, lfsv, r);
      }

      //! get an internal intersections's contribution to lambda
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSV, typename C>
      void lambda_skeleton(const IG& ig,
                           const LFSV& lfsv_s, const LFSV& lfsv_n,
                           C& r_s, C& r_n) const
      {
        ForLoop<LambdaSkeletonOperation, 0, size-1>::
          apply(lops, weights, ig, lfsv_s, lfsv_n, r_s, r_n);
      }

      //! get a boundary intersections's contribution to lambda
      /**
       * \note Summands with zero weight don't contribute to the residual, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSV, typename C>
      void lambda_boundary(const IG& ig, const LFSV& lfsv_s, C& r_s) const
      {
        ForLoop<LambdaBoundaryOperation, 0, size-1>::
          apply(lops, weights, ig, lfsv_s, r_s);
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
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          WeightedVectorAccumulationView<C>& r)
        {
          apply(lops, weights[i]*r.weight(), eg, lfsu, x, lfsv, r.container());
        }
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, typename C::weight_type weight,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          C& r)
        {
          if(weight != K(0)) {
            WeightedVectorAccumulationView<C> view(r, weight);
            LocalAssemblerCallSwitch<Arg, Arg::doAlphaVolume>::
              jacobian_apply_volume(*get<i>(lops), eg, lfsu, x, lfsv, view);
          }
        }
      };

      template<int i>
      struct JacobianApplyVolumePostSkeletonOperation {
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          WeightedVectorAccumulationView<C>& r)
        {
          apply(lops, weights[i]*r.weight(), eg, lfsu, x, lfsv, r.container());
        }
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, typename C::weight_type weight,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          C& r)
        {
          if(weight != K(0)) {
            WeightedVectorAccumulationView<C> view(r, weight);
            LocalAssemblerCallSwitch<Arg, Arg::doAlphaVolumePostSkeleton>::
              jacobian_apply_volume_post_skeleton(*get<i>(lops), eg,
                                                  lfsu, x, lfsv,
                                                  view);
          }
        }
      };

      template<int i>
      struct JacobianApplySkeletonOperation {
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          WeightedVectorAccumulationView<C>& r_s,
                          WeightedVectorAccumulationView<C>& r_n)
        {
          apply(lops, weights[i]*r_s.weight(), weights[i]*r_n.weight(),
                ig,
                lfsu_s, x_s, lfsv_s,
                lfsu_n, x_n, lfsv_n,
                r_s.container(), r_n.container());
        }
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops,
                          typename C::weight_type weight_s,
                          typename C::weight_type weight_n,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          C& r_s, C& r_n)
        {
          if(weight_s != K(0) || weight_n != K(0)) {
            WeightedVectorAccumulationView<C> view_s(r_s, weight_s);
            WeightedVectorAccumulationView<C> view_n(r_n, weight_n);
            LocalAssemblerCallSwitch<Arg, Arg::doAlphaSkeleton>::
              jacobian_apply_skeleton(*get<i>(lops), ig,
                                      lfsu_s, x_s, lfsv_s,
                                      lfsu_n, x_n, lfsv_n,
                                      view_s, view_n);
          }
        }
      };

      template<int i>
      struct JacobianApplyBoundaryOperation {
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          WeightedVectorAccumulationView<C>& r_s)
        {
          apply(lops, weights[i]*r_s.weight(), ig,
                lfsu_s, x_s, lfsv_s,
                r_s.container());
        }
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops,
                          typename C::weight_type weight_s,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          C& r_s)
        {
          if(weight_s != K(0)) {
            WeightedVectorAccumulationView<C> view_s(r_s, weight_s);
            LocalAssemblerCallSwitch<Arg, Arg::doAlphaBoundary>::
              jacobian_apply_boundary(*get<i>(lops), ig,
                                      lfsu_s, x_s, lfsv_s,
                                      view_s);
          }
        }
      };

    public:
      //! apply an element's jacobian
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename C>
      void jacobian_apply_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        C& r) const
      {
        ForLoop<JacobianApplyVolumeOperation, 0, size-1>::
          apply(lops, weights, eg, lfsu, x, lfsv, r);
      }

      //! \brief apply an element's jacobian after the intersections have been
      //!        handled
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename C>
      void jacobian_apply_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        C& r) const
      {
        ForLoop<JacobianApplyVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, weights, eg, lfsu, x, lfsv, r);
      }

      //! apply an internal intersections's jacobians
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename C>
      void jacobian_apply_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        C& r_s, C& r_n) const
      {
        ForLoop<JacobianApplySkeletonOperation, 0, size-1>::
          apply(lops, weights, ig,
                lfsu_s, x_s, lfsv_s,
                lfsu_n, x_n, lfsv_n,
                r_s, r_n);
      }

      //! apply a boundary intersections's jacobian
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename C>
      void jacobian_apply_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        C& r_s) const
      {
        ForLoop<JacobianApplyBoundaryOperation, 0, size-1>::
          apply(lops, weights, ig, lfsu_s, x_s, lfsv_s, r_s);
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
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          WeightedMatrixAccumulationView<C>& m)
        {
          apply(lops, weights[i]*m.weight(), eg, lfsu, x, lfsv, m.container());
        }
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, typename C::weight_type weight,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          C& m)
        {
          if(weight != K(0)) {
            WeightedMatrixAccumulationView<C> view(m, weight);
            LocalAssemblerCallSwitch<Arg, Arg::doAlphaVolume>::
              jacobian_volume(*get<i>(lops), eg, lfsu, x, lfsv, view);
          }
        }
      };

      template<int i>
      struct JacobianVolumePostSkeletonOperation {
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          WeightedMatrixAccumulationView<C>& m)
        {
          apply(lops, weights[i]*m.weight(), eg, lfsu, x, lfsv, m.container());
        }
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, typename C::weight_type weight,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          C& m)
        {
          if(weight != K(0)) {
            WeightedMatrixAccumulationView<C> view(m, weight);
            LocalAssemblerCallSwitch<Arg, Arg::doAlphaVolumePostSkeleton>::
              jacobian_volume_post_skeleton(*get<i>(lops), eg,
                                            lfsu, x, lfsv,
                                            view);
          }
        }
      };

      template<int i>
      struct JacobianSkeletonOperation {
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          WeightedMatrixAccumulationView<C>& m_ss,
                          WeightedMatrixAccumulationView<C>& m_sn,
                          WeightedMatrixAccumulationView<C>& m_ns,
                          WeightedMatrixAccumulationView<C>& m_nn)
        {
          apply(lops,
                weights[i]*m_ss.weight(), weights[i]*m_sn.weight(),
                weights[i]*m_ns.weight(), weights[i]*m_nn.weight(),
                ig,
                lfsu_s, x_s, lfsv_s,
                lfsu_n, x_n, lfsv_n,
                m_ss.container(), m_sn.container(),
                m_ns.container(), m_nn.container());
        }
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops,
                          typename C::weight_type weight_ss,
                          typename C::weight_type weight_sn,
                          typename C::weight_type weight_ns,
                          typename C::weight_type weight_nn,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          C& m_ss, C& m_sn, C& m_ns, C& m_nn)
        {
          if(weight_ss != K(0) || weight_sn != K(0) ||
             weight_ns != K(0) || weight_nn != K(0))
          {
            WeightedMatrixAccumulationView<C> view_ss(m_ss, weight_ss);
            WeightedMatrixAccumulationView<C> view_sn(m_sn, weight_sn);
            WeightedMatrixAccumulationView<C> view_ns(m_ns, weight_ns);
            WeightedMatrixAccumulationView<C> view_nn(m_nn, weight_nn);
            LocalAssemblerCallSwitch<Arg, Arg::doAlphaSkeleton>::
              jacobian_skeleton(*get<i>(lops), ig,
                                lfsu_s, x_s, lfsv_s,
                                lfsu_n, x_n, lfsv_n,
                                view_ss, view_sn, view_ns, view_nn);
          }
        }
      };

      template<int i>
      struct JacobianBoundaryOperation {
        typedef typename tuple_element<i,Args>::type Arg;
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops, const Weights& weights,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          WeightedMatrixAccumulationView<C>& m_ss)
        {
          apply(lops, weights[i]*m_ss.weight(), ig,
                lfsu_s, x_s, lfsv_s,
                m_ss.container());
        }
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename C>
        static void apply(const ArgPtrs& lops,
                          typename C::weight_type weight_ss,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          C& m_ss)
        {
          if(weight_ss != K(0))
          {
            WeightedMatrixAccumulationView<C> view_ss(m_ss, weight_ss);
            LocalAssemblerCallSwitch<Arg, Arg::doAlphaBoundary>::
              jacobian_boundary(*get<i>(lops), ig,
                                lfsu_s, x_s, lfsv_s, view_ss);
          }
        }
      };

    public:
      //! get an element's jacobian
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename C>
      void jacobian_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        C& m) const
      {
        ForLoop<JacobianVolumeOperation, 0, size-1>::
          apply(lops, weights, eg, lfsu, x, lfsv, m);
      }

      //! get an element's jacobian after the intersections have been handled
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename C>
      void jacobian_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        C& m) const
      {
        ForLoop<JacobianVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, weights, eg, lfsu, x, lfsv, m);
      }

      //! apply an internal intersections's jacobians
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename C>
      void jacobian_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        C& m_ss, C& m_sn, C& m_ns, C& m_nn) const
      {
        ForLoop<JacobianSkeletonOperation, 0, size-1>::
          apply(lops, weights, ig,
                lfsu_s, x_s, lfsv_s,
                lfsu_n, x_n, lfsv_n,
                m_ss, m_sn, m_ns, m_nn);
      }

      //! get a boundary intersections's jacobian
      /**
       * \note Summands with zero weight don't contribute to the jacobian, and
       *       the calls to the evaluation methods are eliminated at run-time.
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename C>
      void jacobian_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        C& m_ss) const
      {
        ForLoop<JacobianBoundaryOperation, 0, size-1>::
          apply(lops, weights, ig, lfsu_s, x_s, lfsv_s, m_ss);
      }

      //! \} Methods to extract the jacobian

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for instationary problems
      //! \{
      //

      //! Export type used for time values
      typedef typename tuple_element<0, Args>::type::RealType RealType;

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

#endif // DUNE_PDELAB_LOCALOPERATOR_WEIGHTEDSUM_HH
