// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_SUM_HH
#define DUNE_PDELAB_LOCALOPERATOR_SUM_HH

#include <cstddef>

#include <dune/common/forloop.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/tuples.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/utility.hh>

#include <dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    //! A local operator to take the sum of other local operators
    template<typename L0,
             typename L1 = void, typename L2 = void, typename L3 = void,
             typename L4 = void, typename L5 = void, typename L6 = void,
             typename L7 = void, typename L8 = void, typename L9 = void>
    class SumLocalOperator
    {
      static const std::size_t size =
        is_same<L1, void>::value ? 1 :
        is_same<L2, void>::value ? 2 :
        is_same<L3, void>::value ? 3 :
        is_same<L4, void>::value ? 4 :
        is_same<L5, void>::value ? 5 :
        is_same<L6, void>::value ? 6 :
        is_same<L7, void>::value ? 7 :
        is_same<L8, void>::value ? 8 :
        is_same<L9, void>::value ? 9 : 10;

      typedef
      typename SelectType<size == 1, tuple<L0>,
      typename SelectType<size == 2, tuple<L0, L1>,
      typename SelectType<size == 3, tuple<L0, L1, L2>,
      typename SelectType<size == 4, tuple<L0, L1, L2, L3>,
      typename SelectType<size == 5, tuple<L0, L1, L2, L3, L4>,
      typename SelectType<size == 6, tuple<L0, L1, L2, L3, L4, L5>,
      typename SelectType<size == 7, tuple<L0, L1, L2, L3, L4, L5, L6>,
      typename SelectType<size == 8, tuple<L0, L1, L2, L3, L4, L5, L6, L7>,
      typename SelectType<size == 9, tuple<L0, L1, L2, L3, L4, L5, L6, L7, L8>,
               tuple<L0, L1, L2, L3, L4, L5, L6, L7, L8, L9>
      >::Type>::Type>::Type>::Type>::Type>::Type>::Type>::Type>::Type
      Args;

      template<typename T> struct AddPtr { typedef T* Type; };
      typedef typename ForEachType<AddPtr, Args>::Type ArgPtrs;
      ArgPtrs lops;

    public:
      //////////////////////////////////////////////////////////////////////
      //
      //! \name Construction and modification
      //! \{
      //

      //! construct a 1-component SumLocalOperator
      SumLocalOperator(L0& l0)
        : lops(&l0)
      { }

      //! construct a 2-component SumLocalOperator
      SumLocalOperator(L0& l0, L1& l1)
        : lops(&l0, &l1)
      { }

      //! construct a 3-component SumLocalOperator
      SumLocalOperator(L0& l0, L1& l1, L2& l2)
        : lops(&l0, &l1, &l2)
      { }

      //! construct a 4-component SumLocalOperator
      SumLocalOperator(L0& l0, L1& l1, L2& l2, L3& l3)
        : lops(&l0, &l1, &l2, &l3)
      { }

      //! construct a 5-component SumLocalOperator
      SumLocalOperator(L0& l0, L1& l1, L2& l2, L3& l3, L4& l4)
        : lops(&l0, &l1, &l2, &l3, &l4)
      { }

      //! construct a 6-component SumLocalOperator
      SumLocalOperator(L0& l0, L1& l1, L2& l2, L3& l3, L4& l4, L5& l5)
        : lops(&l0, &l1, &l2, &l3, &l4, &l5)
      { }

      //! construct a 7-component SumLocalOperator
      SumLocalOperator(L0& l0, L1& l1, L2& l2, L3& l3, L4& l4, L5& l5, L6& l6)
        : lops(&l0, &l1, &l2, &l3, &l4, &l5, &l6)
      { }

      //! construct a 8-component SumLocalOperator
      SumLocalOperator(L0& l0, L1& l1, L2& l2, L3& l3, L4& l4, L5& l5, L6& l6,
                       L7& l7)
        : lops(&l0, &l1, &l2, &l3, &l4, &l5, &l6, &l7)
      { }

      //! construct a 9-component SumLocalOperator
      SumLocalOperator(L0& l0, L1& l1, L2& l2, L3& l3, L4& l4, L5& l5, L6& l6,
                       L7& l7, L8& l8)
        : lops(&l0, &l1, &l2, &l3, &l4, &l5, &l6, &l7, &l8)
      { }

      //! construct a 10-component SumLocalOperator
      SumLocalOperator(L0& l0, L1& l1, L2& l2, L3& l3, L4& l4, L5& l5, L6& l6,
                       L7& l7, L8& l8, L9& l9)
        : lops(&l0, &l1, &l2, &l3, &l4, &l5, &l6, &l7, &l8, &l9)
      { }

      //! set the i'th component of the sum
      template<std::size_t i>
      void setSummand(typename tuple_element<i,Args>::type& summand)
      { get<i>(lops) = &summand; }

      //! get the i'th component of the sum
      template<std::size_t i>
      typename tuple_element<i,Args>::type& getSummand()
      { return *get<i>(lops); }

      //! \} Construction and modification

      ////////////////////////////////////////////////////////////////////////
      //
      //! \name Control flags
      //! \{
      //

    private:
      template<typename T1, typename T2>
      struct OrOperation
        : public integral_constant<bool, T1::value || T2:: value>
      { };
      template<template<int> class Value>
      struct AccFlag : public GenericForLoop<OrOperation, Value, 0, size-1>
      { };

      template<int i>
      struct PatternVolumeValue : public integral_constant
      < bool, tuple_element<i, Args>::type::doPatternVolume>
      { };
      template<int i>
      struct PatternVolumePostSkeletonValue : public integral_constant
      < bool, tuple_element<i, Args>::type::doPatternVolumePostSkeleton>
      { };
      template<int i>
      struct PatternSkeletonValue : public integral_constant
      < bool, tuple_element<i, Args>::type::doPatternSkeleton>
      { };
      template<int i>
      struct PatternBoundaryValue : public integral_constant
      < bool, tuple_element<i, Args>::type::doPatternBoundary>
      { };

      template<int i>
      struct AlphaVolumeValue : public integral_constant
      < bool, tuple_element<i, Args>::type::doAlphaVolume>
      { };
      template<int i>
      struct AlphaVolumePostSkeletonValue : public integral_constant
      < bool, tuple_element<i, Args>::type::doAlphaVolumePostSkeleton>
      { };
      template<int i>
      struct AlphaSkeletonValue : public integral_constant
      < bool, tuple_element<i, Args>::type::doAlphaSkeleton>
      { };
      template<int i>
      struct AlphaBoundaryValue : public integral_constant
      < bool, tuple_element<i, Args>::type::doAlphaBoundary>
      { };

      template<int i>
      struct LambdaVolumeValue : public integral_constant
      < bool, tuple_element<i, Args>::type::doLambdaVolume>
      { };
      template<int i>
      struct LambdaVolumePostSkeletonValue : public integral_constant
      < bool, tuple_element<i, Args>::type::doLambdaVolumePostSkeleton>
      { };
      template<int i>
      struct LambdaSkeletonValue : public integral_constant
      < bool, tuple_element<i, Args>::type::doLambdaSkeleton>
      { };
      template<int i>
      struct LambdaBoundaryValue : public integral_constant
      < bool, tuple_element<i, Args>::type::doLambdaBoundary>
      { };

      template<int i>
      struct OneSidedSkeletonRequiredValue : public integral_constant
      < bool, ( ( tuple_element<i, Args>::type::doAlphaSkeleton ||
                  tuple_element<i, Args>::type::doLambdaSkeleton) &&
                ! tuple_element<i, Args>::type::doSkeletonTwoSided)>
      { };
      template<int i>
      struct TwoSidedSkeletonRequiredValue : public integral_constant
      < bool, ( ( tuple_element<i, Args>::type::doAlphaSkeleton ||
                  tuple_element<i, Args>::type::doLambdaSkeleton) &&
                tuple_element<i, Args>::type::doSkeletonTwoSided)>
      { };

    public:
      enum { doPatternVolume             =
             AccFlag<PatternVolumeValue>::value             };
      enum { doPatternVolumePostSkeleton =
             AccFlag<PatternVolumePostSkeletonValue>::value };
      enum { doPatternSkeleton           =
             AccFlag<PatternSkeletonValue>::value           };
      enum { doPatternBoundary           =
             AccFlag<PatternBoundaryValue>::value           };

      enum { doAlphaVolume               =
             AccFlag<AlphaVolumeValue>::value               };
      enum { doAlphaVolumePostSkeleton   =
             AccFlag<AlphaVolumePostSkeletonValue>::value   };
      enum { doAlphaSkeleton             =
             AccFlag<AlphaSkeletonValue>::value             };
      enum { doAlphaBoundary             =
             AccFlag<AlphaBoundaryValue>::value             };

      enum { doLambdaVolume              =
             AccFlag<LambdaVolumeValue>::value              };
      enum { doLambdaVolumePostSkeleton  =
             AccFlag<LambdaVolumePostSkeletonValue>::value  };
      enum { doLambdaSkeleton            =
             AccFlag<LambdaSkeletonValue>::value            };
      enum { doLambdaBoundary            =
             AccFlag<LambdaBoundaryValue>::value            };

      enum { doSkeletonTwoSided          =
             AccFlag<TwoSidedSkeletonRequiredValue>::value  };
      dune_static_assert(!(AccFlag<OneSidedSkeletonRequiredValue>::value &&
                           AccFlag<TwoSidedSkeletonRequiredValue>::value),
                         "Some summands require a one-sided skelton, others a "
                         "two-sided skeleton.  This is not supported.");

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the sparsity pattern
      //! \{
      //

    private:
      template<int i>
      struct PatternVolumeOperation {
        template<typename LFSU, typename LFSV>
        static void apply(ArgPtrs& lops,
                          const LFSU& lfsu, const LFSV& lfsv,
                          LocalSparsityPattern& pattern)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doPatternVolume>::
            pattern_volume(*get<i>(lops), lfsu, lfsv, pattern);
        }
      };

      template<int i>
      struct PatternVolumePostSkeletonOperation {
        template<typename LFSU, typename LFSV>
        static void apply(ArgPtrs& lops,
                          const LFSU& lfsu, const LFSV& lfsv,
                          LocalSparsityPattern& pattern)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doPatternVolumePostSkeleton>::
            pattern_volume_post_skeleton(*get<i>(lops), lfsu, lfsv, pattern);
        }
      };

      template<int i>
      struct PatternSkeletonOperation {
        template<typename LFSU, typename LFSV>
        static void apply(ArgPtrs& lops,
                          const LFSU& lfsu_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const LFSV& lfsv_n,
                          LocalSparsityPattern& pattern_sn,
                          LocalSparsityPattern& pattern_ns)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doPatternSkeleton>::
            pattern_skeleton(*get<i>(lops),
                             lfsu_s, lfsv_s, lfsu_n, lfsv_n,
                             pattern_sn, pattern_ns);
        }
      };

      template<int i>
      struct PatternBoundaryOperation {
        template<typename LFSU, typename LFSV>
        static void apply(ArgPtrs& lops,
                          const LFSU& lfsu_s, const LFSV& lfsv_s,
                          LocalSparsityPattern& pattern_ss)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doPatternBoundary>::
            pattern_boundary(*get<i>(lops), lfsu_s, lfsv_s, pattern_ss);
        }
      };

    public:
      //! get an element's contribution to the sparsity pattern
      /**
       * \param lfsu    LocalFunctionSpace of the trial GridFunctionSpace.
       * \param lfsv    LocalFunctionSpace of the test GridFunctionSpace.
       * \param pattern Local sparsity pattern.
       *
       * \note The method should not clear the pattern; it should just add
       *       its entries to it.
       *
       * This method is controlled by the flag \ref doPatternVolume.  For a
       * given element, it is called *before* the pattern_skeleton() and/or
       * pattern_boundary() methods are called (if they are called at all).
       */
      template<typename LFSU, typename LFSV>
      void pattern_volume
      ( const LFSU& lfsu, const LFSV& lfsv,
        LocalSparsityPattern& pattern)
      {
        ForLoop<PatternVolumeOperation, 0, size-1>::
          apply(lops, lfsu, lfsv, pattern);
      }

      //! \brief get an element's contribution to the sparsity pattern after
      //!        the intersections have been handled
      /**
       * \param lfsu    LocalFunctionSpace of the trial GridFunctionSpace.
       * \param lfsv    LocalFunctionSpace of the test GridFunctionSpace.
       * \param pattern Local sparsity pattern.
       *
       * \note The method should not clear the pattern; it should just add
       *       its entries to it.
       *
       * This method is controlled by the flag \ref doPatternVolume.  For a
       * given element, it is called *before* the pattern_skeleton() and/or
       * pattern_boundary() methods are called (if they are called at all).
       */
      template<typename LFSU, typename LFSV>
      void pattern_volume_post_skeleton
      ( const LFSU& lfsu, const LFSV& lfsv,
        LocalSparsityPattern& pattern)
      {
        ForLoop<PatternVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, lfsu, lfsv, pattern);
      }

      //! get an internal intersection's contribution to the sparsity pattern
      /**
       * \param lfsu_s     LocalFunctionSpace of the trial GridFunctionSpace
       *                   in the inside entity.
       * \param lfsv_s     LocalFunctionSpace of the test GridFunctionSpace
       *                   in the inside entity.
       * \param lfsu_n     LocalFunctionSpace of the trial GridFunctionSpace
       *                   in the outside entity.
       * \param lfsv_n     LocalFunctionSpace of the test GridFunctionSpace
       *                   in the outside entity.
       * \param pattern_sn Local sparsity pattern.
       * \param pattern_ns Local sparsity pattern.
       *
       * \note The method should not clear the patterns; it should just add
       *       its entries to them.
       *
       * This method is controlled by the flag \ref doPatternSkeleton.  For a
       * given element, it's calls happen intermingled with the calls to
       * pattern_boundary(), but after the call to pattern_volume() and before
       * the call to pattern_volume_post_skeleton().
       */
      template<typename LFSU, typename LFSV>
      void pattern_skeleton
      ( const LFSU& lfsu_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const LFSV& lfsv_n,
        LocalSparsityPattern& pattern_sn,
        LocalSparsityPattern& pattern_ns)
      {
        ForLoop<PatternSkeletonOperation, 0, size-1>::
          apply(lops, lfsu_s, lfsv_s, lfsu_n, lfsv_n,
                pattern_sn, pattern_ns);
      }

      //! get a boundary intersection's contribution to the sparsity pattern
      /**
       * \param lfsu_s     LocalFunctionSpace of the trial GridFunctionSpace
       *                   in the inside entity.
       * \param lfsv_s     LocalFunctionSpace of the test GridFunctionSpace
       *                   in the inside entity.
       * \param pattern_ss Local sparsity pattern for the inside entity.
       *
       * \note The method should not clear the pattern; it should just add
       *       its entries to it.
       *
       * This method is controlled by the flag \ref doPatternBoundary.  For a
       * given element, it's calls happen intermingled with the calls to
       * pattern_skeleton(), but after the call to pattern_volume() and before
       * the call to pattern_volume_post_skeleton().
       */
      template<typename LFSU, typename LFSV>
      void pattern_boundary
      ( const LFSU& lfsu_s, const LFSV& lfsv_s,
        LocalSparsityPattern& pattern_ss)
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
      template<int i>
      struct AlphaVolumeOperation {
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename R>
        static void apply(ArgPtrs& lops,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          R& r)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doAlphaVolume>::
            alpha_volume(*get<i>(lops), eg, lfsu, x, lfsv, r);
        }
      };

      template<int i>
      struct AlphaVolumePostSkeletonOperation {
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename R>
        static void apply(ArgPtrs& lops,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          R& r)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doAlphaVolumePostSkeleton>::
            alpha_volume_post_skeleton(*get<i>(lops), eg, lfsu, x, lfsv, r);
        }
      };

      template<int i>
      struct AlphaSkeletonOperation {
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename R>
        static void apply(ArgPtrs& lops,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          R& r_s, R& r_n)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doAlphaSkeleton>::
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
        static void apply(ArgPtrs& lops,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          R& r_s)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doAlphaBoundary>::
            alpha_boundary(*get<i>(lops), ig, lfsu_s, x_s, lfsv_s, r_s);
        }
      };

    public:
      //! get an element's contribution to alpha
      /**
       * \param eg   ElementGeometry describing the entity.
       * \param lfsu LocalFunctionSpace of the trial GridFunctionSpace.
       * \param x    Local position in the trial GridFunctionSpace.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param r    Local part of the residual.
       *
       * \note It is permissible to include contributions of the residual
       *       which are independent of \c x here (they have to be omitted
       *       from lambda_volume() in that case, of course).  This is the
       *       difference to jacobian_apply_volume().
       *
       * \note \c x and \c r are of type std::vector.
       *
       * \note The method should not clear \c r; it should just add its
       *       entries to it.
       *
       * This method is controlled by the flag \ref doAlphaVolume.  For a
       * given element, it is called *before* the alpha_skeleton() and/or
       * alpha_boundary() methods are called (if they are called at all).
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename R>
      void alpha_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        R& r)
      {
        ForLoop<AlphaVolumeOperation, 0, size-1>::
          apply(lops, eg, lfsu, x, lfsv, r);
      }

      //! \brief get an element's contribution to alpha after the
      //!        intersections have been handled
      /**
       * \param eg   ElementGeometry describing the entity.
       * \param lfsu LocalFunctionSpace of the trial GridFunctionSpace.
       * \param x    Local position in the trial GridFunctionSpace.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param r    Local part of the residual.
       *
       * \note It is permissible to include contributions of the residual
       *       which are independent of \c x here (they have to be omitted
       *       from lambda_volume_post_skeleton() in that case, of course).
       *       This is the difference to
       *       jacobian_apply_volume_post_skeleton().
       *
       * \note \c x and \c r are of type std::vector.
       *
       * \note The method should not clear \c r; it should just add its
       *       entries to it.
       *
       * This method is controlled by the flag \ref doAlphaVolumePostSkeleton.
       * For a given element, it is called *after* the alpha_skeleton() and/or
       * alpha_boundary() methods are called (if they are called at all).
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename R>
      void alpha_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        R& r)
      {
        ForLoop<AlphaVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, eg, lfsu, x, lfsv, r);
      }

      //! get an internal intersections's contribution to alpha
      /**
       * \param ig     IntersectionGeometry describing the intersection.
       * \param lfsu_s LocalFunctionSpace of the trial GridFunctionSpace in
       *               the inside entity.
       * \param x_s    Local position in the trial GridFunctionSpace in the
       *               inside entity.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param lfsu_n LocalFunctionSpace of the trial GridFunctionSpace in
       *               the outside entity.
       * \param x_n    Local position in the trial GridFunctionSpace in the
       *               outside entity.
       * \param lfsv_n LocalFunctionSpace of the test GridFunctionSpace in the
       *               outside entity.
       * \param r_s    Local part of the residual in the inside entity.
       * \param r_n    Local part of the residual in the outside entity.
       *
       * \note It is permissible to include contributions of the residual
       *       which are independent of \c x_s and \c x_n here (they have to
       *       be omitted from lambda_skeleton() in that case, of course).
       *       This is the difference to jacobian_apply_skeleton().
       *
       * \note \c x_s, \c x_n, \c r_s and \c r_n are of type std::vector.
       *
       * \note The method should not clear \c r_s and \c r_n; it should just
       *       add its entries to them.
       *
       * This method is controlled by the flag \ref doAlphaSkeleton.  For a
       * given element, it's calls happen intermingled with the calls to
       * alpha_boundary(), but after the call to alpha_volume() and before the
       * call to alpha_volume_post_skeleton().
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename R>
      void alpha_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        R& r_s, R& r_n)
      {
        ForLoop<AlphaSkeletonOperation, 0, size-1>::
          apply(lops, ig,
                lfsu_s, x_s, lfsv_s,
                lfsu_n, x_n, lfsv_n,
                r_s, r_n);
      }

      //! get a boundary intersections's contribution to alpha
      /**
       * \param ig     IntersectionGeometry describing the intersection.
       * \param lfsu_s LocalFunctionSpace of the trial GridFunctionSpace in
       *               the inside entity.
       * \param x_s    Local position in the trial GridFunctionSpace in the
       *               inside entity.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param r_s    Local part of the residual in the inside entity.
       *
       * \note It is permissible to include contributions of the residual
       *       which are independent of \c x_s here (they have to be omitted
       *       from lambda_boundary() in that case, of course).  This is the
       *       difference to jacobian_apply_boundary().
       *
       * \note \c x_s and \c r_s are of type std::vector.
       *
       * \note The method should not clear \c r_s; it should just add its
       *       entries to it.
       *
       * This method is controlled by the flag \ref doAlphaBoundary.  For a
       * given element, it's calls happen intermingled with the calls to
       * alpha_skeleton(), but after the call to alpha_volume() and before the
       * call to alpha_volume_post_skeleton().
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename R>
      void alpha_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        R& r_s)
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
      template<int i>
      struct LambdaVolumeOperation {
        template<typename EG, typename LFSV, typename R>
        static void apply(ArgPtrs& lops,
                          const EG& eg,
                          const LFSV& lfsv,
                          R& r)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doLambdaVolume>::
            lambda_volume(*get<i>(lops), eg, lfsv, r);
        }
      };

      template<int i>
      struct LambdaVolumePostSkeletonOperation {
        template<typename EG, typename LFSV, typename R>
        static void apply(ArgPtrs& lops,
                          const EG& eg,
                          const LFSV& lfsv,
                          R& r)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doLambdaVolumePostSkeleton>::
            lambda_volume_post_skeleton(*get<i>(lops), eg, lfsv, r);
        }
      };

      template<int i>
      struct LambdaSkeletonOperation {
        template<typename IG, typename LFSV, typename R>
        static void apply(ArgPtrs& lops,
                          const IG& ig,
                          const LFSV& lfsv_s,
                          const LFSV& lfsv_n,
                          R& r_s, R& r_n)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doLambdaSkeleton>::
            lambda_skeleton(*get<i>(lops), ig, lfsv_s, lfsv_n, r_s, r_n);
        }
      };

      template<int i>
      struct LambdaBoundaryOperation {
        template<typename IG, typename LFSV, typename R>
        static void apply(ArgPtrs& lops,
                          const IG& ig,
                          const LFSV& lfsv_s,
                          R& r_s)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doLambdaBoundary>::
            lambda_boundary(*get<i>(lops), ig, lfsv_s, r_s);
        }
      };

    public:
      //! get an element's contribution to lambda
      /**
       * \param eg   ElementGeometry describing the entity.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param r    Local part of the residual.
       *
       * \note \c r is of type std::vector.
       *
       * \note The method should not clear \c r; it should just add its
       *       entries to it.
       *
       * This method is controlled by the flag \ref doLambdaVolume.  For a
       * given element, it is called *before* the lambda_skeleton() and/or
       * lambda_boundary() methods are called (if they are called at all).
       */
      template<typename EG, typename LFSV, typename R>
      void lambda_volume(const EG& eg, const LFSV& lfsv, R& r)
      {
        ForLoop<LambdaVolumeOperation, 0, size-1>::
          apply(lops, eg, lfsv, r);
      }

      //! \brief get an element's contribution to lambda after the
      //!        intersections have been handled
      /**
       *
       * \param eg   ElementGeometry describing the entity.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param r    Local part of the residual.
       *
       * \note \c r is of type std::vector.
       *
       * \note The method should not clear \c r; it should just add its
       *       entries to it.
       *
       * This method is controlled by the flag \ref
       * doLambdaVolumePostSkeleton.  For a given element, it is called
       * *after* the lambda_skeleton() and/or lambda_boundary() methods are
       * called (if they are called at all).
       */
      template<typename EG, typename LFSV, typename R>
      void lambda_volume_post_skeleton(const EG& eg, const LFSV& lfsv, R& r)
      {
        ForLoop<LambdaVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, eg, lfsv, r);
      }

      //! get an internal intersections's contribution to lambda
      /**
       * \param ig     IntersectionGeometry describing the intersection.
       *               inside entity.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param lfsv_n LocalFunctionSpace of the test GridFunctionSpace in the
       *               outside entity.
       * \param r_s    Local part of the residual in the inside entity.
       * \param r_n    Local part of the residual in the outside entity.
       *
       * \note \c r_s and \c r_n are of type std::vector.
       *
       * \note The method should not clear \c r_s and \c r_n; it should just
       *       add its entries to them.
       *
       * This method is controlled by the flag \ref doLambdaSkeleton.  For a
       * given element, it's calls happen intermingled with the calls to
       * lambda_boundary(), but after the call to lambda_volume() and before
       * the call to lambda_volume_post_skeleton().
       */
      template<typename IG, typename LFSV, typename R>
      void lambda_skeleton(const IG& ig,
                           const LFSV& lfsv_s, const LFSV& lfsv_n,
                           R& r_s, R& r_n)
      {
        ForLoop<LambdaSkeletonOperation, 0, size-1>::
          apply(lops, ig, lfsv_s, lfsv_n, r_s, r_n);
      }

      //! get a boundary intersections's contribution to lambda
      /**
       * \param ig     IntersectionGeometry describing the intersection.
       *               inside entity.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param r_s    Local part of the residual in the inside entity.
       *
       * \note \c r_s is of type std::vector.
       *
       * \note The method should not clear \c r_s; it should just add its
       *       entries to it.
       *
       * This method is controlled by the flag \ref doLambdaBoundary.  For a
       * given element, it's calls happen intermingled with the calls to
       * lambda_skeleton(), but after the call to lambda_volume() and before
       * the call to lambda_volume_post_skeleton().
       */
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary(const IG& ig, const LFSV& lfsv_s, R& r_s)
      {
        ForLoop<LambdaSkeletonOperation, 0, size-1>::
          apply(lops, ig, lfsv_s, r_s);
      }

      //! \} Methods for the residual -- constant parts

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the application of the jacobian
      //! \{
      //

    private:
      template<int i>
      struct JacobianApplyVolumeOperation {
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename Y>
        static void apply(ArgPtrs& lops,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          Y& y)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doAlphaVolume>::
            jacobian_apply_volume(*get<i>(lops), eg, lfsu, x, lfsv, y);
        }
      };

      template<int i>
      struct JacobianApplyVolumePostSkeletonOperation {
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename Y>
        static void apply(ArgPtrs& lops,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          Y& y)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doAlphaVolumePostSkeleton>::
            jacobian_apply_volume_post_skeleton(*get<i>(lops), eg,
                                                lfsu, x, lfsv,
                                                y);
        }
      };

      template<int i>
      struct JacobianApplySkeletonOperation {
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename Y>
        static void apply(ArgPtrs& lops,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          Y& y_s, Y& y_n)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doAlphaSkeleton>::
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
        static void apply(ArgPtrs& lops,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          Y& y_s)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doAlphaBoundary>::
            jacobian_apply_boundary(*get<i>(lops), ig,
                                    lfsu_s, x_s, lfsv_s,
                                    y_s);
        }
      };

    public:
      //! apply an element's jacobian
      /**
       * \param eg   ElementGeometry describing the entity.
       * \param lfsu LocalFunctionSpace of the trial GridFunctionSpace.
       * \param x    Local position in the trial GridFunctionSpace.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param y    Where to store the result.
       *
       * \note This is different from alpha_volume(), since the result will be
       *       linear in \c x, whereas alpha_volume() may include
       *       contributions to the the residual which are constant in \c x.
       *
       * \note \c x and \c y are of type std::vector.
       *
       * \note The method should not clear \c y; it should just add its
       *       entries to it.
       *
       * \note \c x is both the position where the jacobian is evaluated (for
       *       non-linear problems) as well as the vector the jacobian is
       *       applied to.
       *
       * This method is controlled by the flag \ref doAlphaVolume.  For a
       * given element, it is called *before* the jacobian_apply_skeleton()
       * and/or jacobian_apply_boundary() methods are called (if they are
       * called at all).
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        Y& y)
      {
        ForLoop<JacobianApplyVolumeOperation, 0, size-1>::
          apply(lops, eg, lfsu, x, lfsv, y);
      }

      //! \brief apply an element's jacobian after the intersections have been
      //!        handled
      /**
       * \param eg   ElementGeometry describing the entity.
       * \param lfsu LocalFunctionSpace of the trial GridFunctionSpace.
       * \param x    Local position in the trial GridFunctionSpace.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param y    Where to store the result.
       *
       * \note This is different from alpha_volume_post_skeleton(), since the
       *       result will be linear in \c x, whereas
       *       alpha_volume_post_skeleton() may include contributions to the
       *       the residual which are constant in \c x.
       *
       * \note \c x and \c y are of type std::vector.
       *
       * \note The method should not clear \c y; it should just add its
       *       entries to it.
       *
       * \note \c x is both the position where the jacobian is evaluated (for
       *       non-linear problems) as well as the vector the jacobian is
       *       applied to.
       *
       * This method is controlled by the flag \ref doAlphaVolumePostSkeleton.
       * For a given element, it is called *after* the
       * jacobian_apply_skeleton() and/or jacobian_apply_boundary() methods
       * are called (if they are called at all).
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        Y& y)
      {
        ForLoop<JacobianApplyVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, eg, lfsu, x, lfsv, y);
      }

      //! apply an internal intersections's jacobians
      /**
       * \param ig     IntersectionGeometry describing the intersection.
       * \param lfsu_s LocalFunctionSpace of the trial GridFunctionSpace in
       *               the inside entity.
       * \param x_s    Local position in the trial GridFunctionSpace in the
       *               inside entity.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param lfsu_n LocalFunctionSpace of the trial GridFunctionSpace in
       *               the outside entity.
       * \param x_n    Local position in the trial GridFunctionSpace in the
       *               outside entity.
       * \param lfsv_n LocalFunctionSpace of the test GridFunctionSpace in the
       *               outside entity.
       * \param y_s    Where to store the inside entity's result.
       * \param y_n    Where to store the outside entity's result.
       *
       * \note This is different from alpha_skeleton(), since the result will
       *       be linear in \c x_s and \c x_n, whereas alpha_skeleton() may
       *       include contributions to the the residual which are constant in
       *       \c x_s and \c x_n.
       *
       * \note \c x_s, \c x_n, \c y_s and \c y_n are of type std::vector.
       *
       * \note The method should not clear \c y_s and \c y_n; it should just
       *       add its entries to them.
       *
       * \note \c x_s and \c x_n are both the positions where the jacobian is
       *       evaluated (for non-linear problems) as well as the vectors the
       *       jacobian is applied to.
       *
       * This method is controlled by the flag \ref doAlphaSkeleton.  For a
       * given element, it's calls happen intermingled with the calls to
       * jacobian_apply_boundary(), but after the call to
       * jacobian_apply_volume() and before the call to
       * jacobian_apply_volume_post_skeleton().
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n)
      {
        ForLoop<JacobianApplySkeletonOperation, 0, size-1>::
          apply(lops, ig,
                lfsu_s, x_s, lfsv_s,
                lfsu_n, x_n, lfsv_n,
                y_s, y_n);
      }

      //! apply a boundary intersections's jacobian
      /**
       * \param ig     IntersectionGeometry describing the intersection.
       * \param lfsu_s LocalFunctionSpace of the trial GridFunctionSpace in
       *               the inside entity.
       * \param x_s    Local position in the trial GridFunctionSpace in the
       *               inside entity.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param r_s    Local part of the residual in the inside entity.
       *
       * \note This is different from alpha_boundary(), since the result will
       *       be linear in \c x, whereas alpha_boundary() may include
       *       contributions to the the residual which are constant in \c x.
       *
       * \note \c x_s and \c y_s are of type std::vector.
       *
       * \note The method should not clear \c y_s; it should just add its
       *       entries to it.
       *
       * \note \c x_s is both the position where the jacobian is evaluated
       *       (for non-linear problems) as well as the vector the jacobian is
       *       applied to.
       *
       * This method is controlled by the flag \ref doAlphaBoundary.  For a
       * given element, it's calls happen intermingled with the calls to
       * jacobian_apply_skeleton(), but after the call to
       * jacobian_apply_volume() and before the call to
       * jacobian_apply_volume_post_skeleton().
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        Y& y_s)
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
      template<int i>
      struct JacobianVolumeOperation {
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename R>
        static void apply(ArgPtrs& lops,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          LocalMatrix<R>& mat)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doAlphaVolume>::
            jacobian_volume(*get<i>(lops), eg, lfsu, x, lfsv, mat);
        }
      };

      template<int i>
      struct JacobianVolumePostSkeletonOperation {
        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename R>
        static void apply(ArgPtrs& lops,
                          const EG& eg,
                          const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          LocalMatrix<R>& mat)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doAlphaVolumePostSkeleton>::
            jacobian_volume_post_skeleton(*get<i>(lops), eg,
                                          lfsu, x, lfsv,
                                          mat);
        }
      };

      template<int i>
      struct JacobianSkeletonOperation {
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename R>
        static void apply(ArgPtrs& lops,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          LocalMatrix<R>& mat_ss, LocalMatrix<R>& mat_sn,
                          LocalMatrix<R>& mat_ns, LocalMatrix<R>& mat_nn)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doAlphaSkeleton>::
            jacobian_skeleton(*get<i>(lops), ig,
                              lfsu_s, x_s, lfsv_s,
                              lfsu_n, x_n, lfsv_n,
                              mat_ss, mat_sn, mat_ns, mat_nn);
        }
      };

      template<int i>
      struct JacobianBoundaryOperation {
        template<typename IG, typename LFSU, typename X, typename LFSV,
                 typename R>
        static void apply(ArgPtrs& lops,
                          const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          LocalMatrix<R>& mat_ss)
        {
          LocalAssemblerCallSwitch<typename tuple_element<i,Args>::type,
            tuple_element<i,Args>::type::doAlphaBoundary>::
            jacobian_boundary(*get<i>(lops), ig, lfsu_s, x_s, lfsv_s, mat_ss);
        }
      };

    public:
      //! get an element's jacobian
      /**
       * \param eg   ElementGeometry describing the entity.
       * \param lfsu LocalFunctionSpace of the trial GridFunctionSpace.
       * \param x    Local position in the trial GridFunctionSpace.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param mat  Where to store the contribution to the jacobian.
       *
       * \note \c x is of type std::vector.
       *
       * \note The method should not clear \c mat; it should just add its
       *       entries to it.
       *
       * This method is controlled by the flag \ref doAlphaVolume.  For a
       * given element, it is called *before* the jacobian_skeleton() and/or
       * jacobian_boundary() methods are called (if they are called at all).
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename R>
      void jacobian_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        LocalMatrix<R>& mat)
      {
        ForLoop<JacobianVolumeOperation, 0, size-1>::
          apply(lops, eg, lfsu, x, lfsv, mat);
      }

      //! get an element's jacobian after the intersections have been handled
      /**
       * \param eg   ElementGeometry describing the entity.
       * \param lfsu LocalFunctionSpace of the trial GridFunctionSpace.
       * \param x    Local position in the trial GridFunctionSpace.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param mat  Where to store the contribution to the jacobian.
       *
       * \note \c x is of type std::vector.
       *
       * \note The method should not clear \c mat; it should just add its
       *       entries to it.
       *
       * This method is controlled by the flag \ref doAlphaVolumePostSkeleton.
       * For a given element, it is called *after* the jacobian_skeleton()
       * and/or jacobian_boundary() methods are called (if they are called at
       * all).
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename R>
      void jacobian_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        LocalMatrix<R>& mat)
      {
        ForLoop<JacobianVolumePostSkeletonOperation, 0, size-1>::
          apply(lops, eg, lfsu, x, lfsv, mat);
      }

      //! apply an internal intersections's jacobians
      /**
       * \param ig     IntersectionGeometry describing the intersection.
       * \param lfsu_s LocalFunctionSpace of the trial GridFunctionSpace in
       *               the inside entity.
       * \param x_s    Local position in the trial GridFunctionSpace in the
       *               inside entity.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param lfsu_n LocalFunctionSpace of the trial GridFunctionSpace in
       *               the outside entity.
       * \param x_n    Local position in the trial GridFunctionSpace in the
       *               outside entity.

       * \param lfsv_n LocalFunctionSpace of the test GridFunctionSpace in the
       *               outside entity.
       * \param mat_ss Where to store the contribution to the inside entity's
       *               jacobian.
       * \param mat_sn Where to store the contribution to the interaction
       *               jacobian between the inside and the outside entity.
       * \param mat_ns Where to store the contribution to the interaction
       *               jacobian between the outside and the inside entity.
       * \param mat_nn Where to store the contribution to the outside entity's
       *               jacobian.
       *
       * \note \c x_s and \c x_n are of type std::vector.
       *
       * \note The method should not clear \c mat_ss, \c mat_sn, \c mat_ns, \c
       *       mat_nn; it should just add its entries to them.
       *
       * This method is controlled by the flag \ref doAlphaSkeleton.  For a
       * given element, it's calls happen intermingled with the calls to
       * jacobian_boundary(), but after the call to jacobian_volume() and
       * before the call to jacobian_volume_post_skeleton().
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename R>
      void jacobian_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        LocalMatrix<R>& mat_ss, LocalMatrix<R>& mat_sn,
        LocalMatrix<R>& mat_ns, LocalMatrix<R>& mat_nn)
      {
        ForLoop<JacobianSkeletonOperation, 0, size-1>::
          apply(lops, ig,
                lfsu_s, x_s, lfsv_s,
                lfsu_n, x_n, lfsv_n,
                mat_ss, mat_sn, mat_ns, mat_nn);
      }

      //! get a boundary intersections's jacobian
      /**
       * \param ig     IntersectionGeometry describing the intersection.
       * \param lfsu_s LocalFunctionSpace of the trial GridFunctionSpace in
       *               the inside entity.
       * \param x_s    Local position in the trial GridFunctionSpace in the
       *               inside entity.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param mat_ss Where to store the contribution to the inside entity's
       *               jacobian.
       *
       * \note \c x_s is of type std::vector.
       *
       * \note The method should not clear \c mat_ss; it should just add its
       *       entries to it.
       *
       * This method is controlled by the flag \ref doAlphaBoundary.  For a
       * given element, it's calls happen intermingled with the calls to
       * jacobian_skeleton(), but after the call to jacobian_volume() and
       * before the call to jacobian_volume_post_skeleton().
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename R>
      void jacobian_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        LocalMatrix<R>& mat_ss)
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

      typedef typename L0::RealType RealType;

    private:
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
      /**
       * This method set the time for subsequent calls to the alpha_*(),
       * lambda_*(), jacobian_*() and jacobian_apply_*() methods.
       *
       * \note For ExplicitOneStepMethod the time given here in the
       *       first stage may be incorrect, since the time step size is
       *       only finally determined after the first stage has been
       *       assembled.
       */
      void setTime (RealType t)
      {
        ForLoop<SetTimeOperation, 0, size-1>::apply(lops, t);
      }

      //! get current time
      /**
       * \return The time previously set by setTime().
       */
      RealType getTime () const
      {
        return get<0>(lops)->getTime();
      }

      //! to be called once before each time step
      /**
       * \param time   Time at beginning of the step.
       * \param dt     Size of time step.
       * \param stages Number of stages to do in the step.  For the
       *               MultiStepMethod this is always 1.
       *
       * \note For ExplicitOneStepMethod the dt given here may be
       *       incorrect, since the time step size is only finally
       *       determined after the first stage has been assembled.
       *
       * \note For the MultiStepMethod the number of stages is given as
       *       1.  There are no since there are no times of evaluation
       *       in the middle of the step, a multi-step method is similar
       *       to a one step method with one stage.
       */
      void preStep (RealType time, RealType dt, int stages)
      {
        ForLoop<PreStepOperation, 0, size-1>::apply(lops, time, dt, stages);
      }

      //! to be called once at the end of each time step
      /**
       * \note With the OneStepMethod and the ExplicitOneStepMetod, for
       *       reasons unknown this is only called for temporal but not
       *       for spatial local operators.  With the MultiStepMethod
       *       this *is* called for all local operators.
       */
      void postStep ()
      {
        ForLoop<PostStepOperation, 0, size-1>::apply(lops);
      }

      //! to be called once before each stage
      /**
       * \param time Time of the stage
       * \param r    Number of the stage, r ∈ [1, nstages] inclusive,
       *             where nstages is the number of stage in the step
       *             given in the previous call to preStep()
       *
       * \note For ExplicitOneStepMethod the time given here for stage 1
       *       may be incorrect, since the time step size is only
       *       finally determined after the first stage has been
       *       assembled.
       *
       * \note For the MultiStepMethod, this is called once after
       *       preStep() with r=1.
       */
      void preStage (RealType time, int r)
      {
        ForLoop<PreStageOperation, 0, size-1>::apply(lops, time, r);
      }

      //! get current stage
      /**
       * \return The current stage number previously set by preStage().
       */
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
       * \note Only used by the ExplicitOneStepMethod.
       *
       * This may be called on the spatial local operator in the case of
       * an explicit one step scheme.  It is called after stage 1 has
       * been assembled (so the time given to preStep() may not apply
       * anymore in this case).  All the alpha_*() and lambda_*()
       * methods should have been called, so they are a good place to
       * generate the information returned here.
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
