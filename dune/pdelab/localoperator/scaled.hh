// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_SCALED_HH
#define DUNE_PDELAB_LOCALOPERATOR_SCALED_HH

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    //! A local operator that scales the result of another local operator
    /**
     * \nosubgrouping
     *
     * This local operator takes another local operator as a backend and a
     * scaling factor.  It forwards calls to the evaluation methods to the
     * backend local operator, but scales the result by the scaling factor
     * before returning it.  Calls to the instationary methods setTime()
     * et. al. are forwarded as well, without modification to the behaviour.
     *
     * If the scaling factor equal zero, calls to the evaluation methods of
     * the backend operator are eliminated at run time.  Note that "equals
     * zero" is determined by the expression \c factor==0, so if you have
     * values that are nearly zero and you want them to be treated as zero,
     * you should force the factor to zero yourself before passing it to this
     * class.
     *
     * \tparam Backend Type of the backend operator.
     * \tparam Factor  Type of the scaling factor.
     * \tparam Time    Type of time values.
     */
    template<typename Backend, typename Factor, typename Time = double>
    class ScaledLocalOperator
    // do *not* derive from LocalOperatorDefaultFlags -- we need to take care
    // of all the possible flags ourselves.
    //
    // If there is a new flag added to LocalOperatorDefaultFlags, it will
    // probably be false by default.  If we would derive from
    // LocalOperatorDefaultFlags, we would derive from that flag -- and claim
    // this flag to be false no matter what the underlying local operator
    // says.  It that flag controls a method that the underlying operator
    // implements, that method would never get called, and no diagnostic would
    // be given.  The only way to notice this situation would be by noticing
    // that the solution of the simulation is errornous.
    //
    // By not deriving from LocalOperatorDefaultFlags, we do not inherit any
    // new flags.  If anything tries to use this new flag, this will result in
    // a compilation error.  This way the missing implementation in this class
    // gets noticed, and the flag can be correctly forwarded from the
    // underlying local operator.
    //
    // The same argument applies for not deriving from
    // InstationaryLocalOperatorDefaultMethods.
    {
      Factor factor;
      Backend* bp;

    public:
      //////////////////////////////////////////////////////////////////////
      //
      //! \name Construction and modification
      //! \{
      //

      //! construct a ScaledLocalOperator
      /**
       * \param backend Reference to the backend local operator
       * \param factor_ Scaling factor
       */
      ScaledLocalOperator(Backend& backend, Factor factor_ = 0)
        : factor(factor_), bp(&backend) { }
      //! construct a ScaledLocalOperator
      /**
       * \note This constructor does not set a backend local operator -- you
       *       *must* call setBackend() before you can use this object
       *
       * \param factor_ Scaling factor
       */
      ScaledLocalOperator(Factor factor_ = 0) : factor(factor_), bp(0) { }

      //! set the scaling factor
      void setFactor(Factor factor_) { factor = factor_; }
      //! get the scaling factor
      Factor getFactor() const { return factor; }

      //! set the backend
      void setBackend(Backend& backend) { bp = &backend; }
      //! get a reference to the backend
      Backend& getBackend() const { return *bp; }

      //! \} Construction and modification

      ////////////////////////////////////////////////////////////////////////
      //
      //! \name Control flags
      //! \{
      //

      //! \brief Whether to assemble the pattern on the elements, i.e. whether
      //!        or not pattern_volume() should be called.
      enum { doPatternVolume = Backend::doPatternVolume };
      //! \brief Whether to assemble the pattern on the elements after the
      //!        skeleton has been handled, i.e. whether or not
      //!        pattern_volume_post_skeleton() should be called.
      enum {
        doPatternVolumePostSkeleton = Backend::doPatternVolumePostSkeleton
      };
      //! \brief Whether to assemble the pattern on the interior
      //!        intersections, i.e. whether or not pattern_skeleton()
      //!        should be called.
      enum { doPatternSkeleton = Backend::doPatternSkeleton };
      //! \brief Whether to assemble the pattern on the boundary
      //!        intersections, i.e. whether or not pattern_boundary()
      //!        should be called.
      enum { doPatternBoundary = Backend::doPatternBoundary };

      //! \brief Whether to call the local operator's alpha_volume(),
      //!        jacobian_apply_volume() and jacobian_volume().
      enum { doAlphaVolume = Backend::doAlphaVolume };
      //! \brief Whether to call the local operator's
      //!        alpha_volume_post_skeleton(),
      //!        jacobian_apply_volume_post_skeleton() and
      //!        jacobian_volume_post_skeleton().
      enum { doAlphaVolumePostSkeleton = Backend::doAlphaVolumePostSkeleton };
      //! \brief Whether to call the local operator's alpha_skeleton(),
      //!        jacobian_apply_skeleton() and jacobian_skeleton().
      enum { doAlphaSkeleton = Backend::doAlphaSkeleton };
      //! \brief Whether to call the local operator's alpha_boundary(),
      //!        jacobian_apply_boundary() and jacobian_boundary().
      enum { doAlphaBoundary = Backend::doAlphaBoundary };

      //! \brief Whether to call the local operator's lambda_volume().
      enum { doLambdaVolume = Backend::doLambdaVolume };
      //! \brief Whether to call the local operator's
      //!        lambda_volume_post_skeleton().
      enum {
        doLambdaVolumePostSkeleton = Backend::doLambdaVolumePostSkeleton
      };
      //! \brief Whether to call the local operator's lambda_skeleton().
      enum { doLambdaSkeleton = Backend::doLambdaSkeleton };
      //! \brief Whether to call the local operator's lambda_boundary().
      enum { doLambdaBoundary = Backend::doLambdaBoundary };

      //! \brief Whether to visit the skeleton methods from both sides
      enum { doSkeletonTwoSided = Backend::doSkeletonTwoSided };

      //! \} Control flags

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the sparsity pattern
      //! \{
      //

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
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_volume
      ( const LFSU& lfsu, const LFSV& lfsv,
        LocalPattern& pattern) const
      {
        if(factor != 0)
          bp->pattern_volume(lfsu, lfsv, pattern);
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
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_volume_post_skeleton
      ( const LFSU& lfsu, const LFSV& lfsv,
        LocalPattern& pattern) const
      {
        if(factor != 0)
          bp->pattern_volume_post_skeleton(lfsu, lfsv, pattern);
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
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_skeleton
      ( const LFSU& lfsu_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const LFSV& lfsv_n,
        LocalPattern& pattern_sn,
        LocalPattern& pattern_ns) const
      {
        if(factor != 0)
          bp->pattern_skeleton(lfsu_s, lfsv_s, lfsu_n, lfsv_n,
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
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_boundary
      ( const LFSU& lfsu_s, const LFSV& lfsv_s,
        LocalPattern& pattern_ss) const
      {
        if(factor != 0)
          bp->pattern_boundary(lfsu_s, lfsv_s, pattern_ss);
      }

      //! \} Methods for the sparsity pattern

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the residual -- non-constant parts
      //! \{
      //

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
        R& r) const
      {
        if(factor != 0) {
          typename R::WeightedAccumulationView
            my_r(r.weightedAccumulationView(factor));
          bp->alpha_volume(eg, lfsu, x, lfsv, my_r);
        }
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
        R& r) const
      {
        if(factor != 0) {
          typename R::WeightedAccumulationView
            my_r(r.weightedAccumulationView(factor));
          bp->alpha_volume_post_skeleton(eg, lfsu, x, lfsv, my_r);
        }
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
        R& r_s, R& r_n) const
      {
        if(factor != 0) {
          typename R::WeightedAccumulationView
            my_r_s(r_s.weightedAccumulationView(factor));
          typename R::WeightedAccumulationView
            my_r_n(r_n.weightedAccumulationView(factor));
          bp->alpha_skeleton(ig,
                             lfsu_s, x_s, lfsv_s,
                             lfsu_n, x_n, lfsv_n,
                             my_r_s, my_r_n);
        }
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
        R& r_s) const
      {
        if(factor != 0) {
          typename R::WeightedAccumulationView
            my_r_s(r_s.weightedAccumulationView(factor));
          bp->alpha_boundary(ig, lfsu_s, x_s, lfsv_s, my_r_s);
        }
      }

      //! \} Methods for the residual -- non-constant parts

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the residual -- constant parts
      //! \{
      //

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
      void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
      {
        if(factor != 0) {
          typename R::WeightedAccumulationView
            my_r(r.weightedAccumulationView(factor));
          bp->lambda_volume(eg, lfsv, my_r);
        }
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
      void lambda_volume_post_skeleton(const EG& eg,
                                       const LFSV& lfsv,
                                       R& r) const
      {
        if(factor != 0) {
          typename R::WeightedAccumulationView
            my_r(r.weightedAccumulationView(factor));
          bp->lambda_volume_post_skeleton(eg, lfsv, my_r);
        }
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
                           R& r_s, R& r_n) const
      {
        if(factor != 0) {
          typename R::WeightedAccumulationView
            my_r_s(r_s.weightedAccumulationView(factor));
          typename R::WeightedAccumulationView
            my_r_n(r_n.weightedAccumulationView(factor));
          bp->lambda_skeleton(ig, lfsv_s, lfsv_n, my_r_s, my_r_n);
        }
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
      void lambda_boundary(const IG& ig, const LFSV& lfsv_s, R& r_s) const
      {
        if(factor != 0) {
          typename R::WeightedAccumulationView
            my_r_s(r_s.weightedAccumulationView(factor));
          bp->lambda_boundary(ig, lfsv_s, my_r_s);
        }
      }

      //! \} Methods for the residual -- constant parts

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the application of the jacobian
      //! \{
      //

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
        Y& y) const
      {
        if(factor != 0) {
          typename Y::WeightedAccumulationView
            my_y(y.weightedAccumulationView(factor));
          bp->jacobian_apply_volume(eg, lfsu, x, lfsv, my_y);
        }
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
        Y& y) const
      {
        if(factor != 0) {
          typename Y::WeightedAccumulationView
            my_y(y.weightedAccumulationView(factor));
          bp->jacobian_apply_volume_post_skeleton(eg, lfsu, x, lfsv, my_y);
        }
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
        Y& y_s, Y& y_n) const
      {
        if(factor != 0) {
          typename Y::WeightedAccumulationView
            my_y_s(y_s.weightedAccumulationView(factor));
          typename Y::WeightedAccumulationView
            my_y_n(y_n.weightedAccumulationView(factor));
          bp->jacobian_apply_skeleton(ig,
                                      lfsu_s, x_s, lfsv_s,
                                      lfsu_n, x_n, lfsv_n,
                                      my_y_s, my_y_n);
        }
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
       * \param y_s    Local part of the residual in the inside entity.
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
        Y& y_s) const
      {
        if(factor != 0) {
          typename Y::WeightedAccumulationView
            my_y_s(y_s.weightedAccumulationView(factor));
          bp->jacobian_apply_boundary(ig, lfsu_s, x_s, lfsv_s, my_y_s);
        }
      }

      //! \} Methods for the application of the jacobian

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods to extract the jacobian
      //! \{
      //

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
               typename M>
      void jacobian_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        M& mat) const
      {
        if(factor != 0) {
          typename M::WeightedAccumulationView
            my_mat(mat.weightedAccumulationView(factor));
          bp->jacobian_volume(eg, lfsu, x, lfsv, my_mat);
        }
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
               typename M>
      void jacobian_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        M& mat) const
      {
        if(factor != 0) {
          typename M::WeightedAccumulationView
            my_mat(mat.weightedAccumulationView(factor));
          bp->jacobian_volume_post_skeleton(eg, lfsu, x, lfsv, my_mat);
        }
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
               typename M>
      void jacobian_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        M& mat_ss, M& mat_sn, M& mat_ns, M& mat_nn) const
      {
        if(factor != 0) {
          typename M::WeightedAccumulationView
            my_mat_ss(mat_ss.weightedAccumulationView(factor));
          typename M::WeightedAccumulationView
            my_mat_sn(mat_sn.weightedAccumulationView(factor));
          typename M::WeightedAccumulationView
            my_mat_ns(mat_ns.weightedAccumulationView(factor));
          typename M::WeightedAccumulationView
            my_mat_nn(mat_nn.weightedAccumulationView(factor));
          bp->jacobian_skeleton(ig,
                                lfsu_s, x_s, lfsv_s,
                                lfsu_n, x_n, lfsv_n,
                                my_mat_ss, my_mat_sn, my_mat_ns, my_mat_nn);
        }
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
               typename M>
      void jacobian_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        M& mat_ss) const
      {
        if(factor != 0) {
          typename M::WeightedAccumulationView
            my_mat_ss(mat_ss.weightedAccumulationView(factor));
          bp->jacobian_boundary(ig, lfsu_s, x_s, lfsv_s, my_mat_ss);
        }
      }

      //! \} Methods to extract the jacobian

      ////////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for temporal local operators
      //! \{

      typedef Time RealType;

      //! set time for subsequent evaluation
      /**
       * This method set the time for subsequent calls to the alpha_*(),
       * lambda_*(), jacobian_*() and jacobian_apply_*() methods.
       *
       * \note For ExplicitOneStepMethod the time given here in the first
       *       stage may be incorrect, since the time step size is only
       *       finally determined after the first stage has been assembled.
       */
      void setTime (Time t) { bp->setTime(t); }

      //! get current time
      /**
       * \return The time previously set by setTime().
       */
      Time getTime () const { return bp->getTime(); }

      //! to be called once before each time step
      /**
       * \param time   Time at beginning of the step.
       * \param dt     Size of time step.
       * \param stages Number of stages to do in the step.  For the
       *               MultiStepMethod this is always 1.
       *
       * \note For ExplicitOneStepMethod the dt given here may be incorrect,
       *       since the time step size is only finally determined after the
       *       first stage has been assembled.
       *
       * \note For the MultiStepMethod the number of stages is given as 1.
       *       There are no since there are no times of evaluation in the
       *       middle of the step, a multi-step method is similar to a one
       *       step method with one stage.
       */
      void preStep (Time time, Time dt, int stages)
      { bp->preStep(time, dt, stages); }

      //! to be called once at the end of each time step
      /**
       * \note With the OneStepMethod and the ExplicitOneStepMetod, for
       *       reasons unknown this is only called for temporal but not
       *       for spatial local operators.  With the MultiStepMethod
       *       this *is* called for all local operators.
       */
      void postStep () { bp->postStep(); }

      //! to be called once before each stage
      /**
       * \param time Time of the stage
       * \param r    Number of the stage, r ∈ [1, nstages] inclusive, where
       *             nstages is the number of stage in the step given in the
       *             previous call to preStep()
       *
       * \note For ExplicitOneStepMethod the time given here for stage 1 may
       *       be incorrect, since the time step size is only finally
       *       determined after the first stage has been assembled.
       *
       * \note For the MultiStepMethod, this is called once after preStep()
       *       with r=1.
       */
      void preStage (Time time, int r) { bp->preStage(time, r); }

      //! get current stage
      /**
       * \return The current stage number previously set by preStage().
       */
      int getStage () const { return bp->getStage(); }

      //! to be called once at the end of each stage
      void postStage () { bp->postStage(); }

      //! to be called after stage 1
      /**
       * \note Only used by the ExplicitOneStepMethod.
       *
       * This may be called on the spatial local operator in the case of an
       * explicit one step scheme.  It is called after stage 1 has been
       * assembled (so the time given to preStep() may not apply anymore in
       * this case).  All the alpha_*() and lambda_*() methods should have
       * been called, so they are a good place to generate the information
       * returned here.
       */
      Time suggestTimestep (Time dt) const
      { return bp->suggestTimestep(dt); }

      //! \} Methods for temporal local operators

    };

    //! \} group LocalOperator
  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_SCALED_HH
