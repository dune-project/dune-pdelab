// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_INTERFACE_HH
#define DUNE_PDELAB_LOCALOPERATOR_INTERFACE_HH

#include <dune/pdelab/localoperator/flags.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    //! Class to document the stationary local operator interface
    /**
     * \nosubgrouping
     *
     * This class is for documentation purposes only.  Each method given here
     * is controlled by a flag.  If the corresponding flag for a method is
     * false, that method is never called and it is permissible for the method
     * to be missing entirely from the local operator.  The flags are those
     * from LocalOperatorDefaultFlags.
     *
     * There are five categories of methods, denoted by the first part of
     * their name:
     * \li \c pattern_*(): Methods for the sparsity pattern, controlled by the
     *     \c doPattern* flags,
     * \li \c alpha_*(): Methods for the non-constant parts of the residual,
     *     controlled by the \c doAlpha* flags,
     * \li \c lambda_*(): Methods for the constant parts of the residual,
     *     controlled by the \c doLambda* flags,
     * \li \c jacobian_apply_*(): Methods for the application of the jacobian,
     *     controlled by the \c doAlpha* flags, and finally
     * \li \c jacobian_*(): Methods to extract the jacobian, controlled by the
     *     \c doAlpha* flags.
     *
     * There are four classes of methods, denoted by the last part of their
     * name:
     * \li \c *_volume(): methods called on the entities before iterating over
     *     the intersections, controlled by the \c do*Volume flags,
     * \li \c *_volume_post_skeleton(): methods called on the entities after
     *     iterating over the intersections, controlled by the \c
     *     do*VolumePostSkeleton flags,
     * \li \c *_skeleton(): methods called on the interior intersections,
     *     controlled by the \c do*Skeleton flags, and finally
     * \li \c *_boundary(): methods called on the boundary intersections,
     *     controlled by the \c do*Boundary flags.
     *
     * Not all combinations of categories and methods do actually exist.
     *
     * To assemble the global sparsity pattern, residual or jacobian, the
     * GridOperator iterates over the elements of the grid.  For each
     * element, it will call the appropriate \c *_volume() method.  Then it
     * will iterate through the elements intersections and call the
     * appropriate \c *_skeleton() or \c *_boundary() methods on the
     * intersection.  Finally it will call the appropriate
     * *_volume_post_skeleton() method.
     *
     * The special flag \ref doSkeletonTwoSided controls whether each interior
     * intersection is visited once or twice.  If it is true, each
     * intersection \em may be given to \c *_skeleton() twice -- the second
     * time with the meaning of \em inside and \em outside exchanged.  Note the
     * "may": In the parallel case only interior entities are visited, so
     * intersections at a processor boundary will only be visited once per
     * processor in any case.
     *
     * If \ref doSkeletonTwoSided is false (the default), each intersection
     * will only be visited once -- and the orientation in which it is visited
     * is left unspecified, except that the inside entity will always be an
     * interior entity.  Note that it will be visited once per process, so
     * intersecions at processor boundaries are still visited twice when all
     * processes are considered.
     *
     * The \c alpha and \c lambda categories are a bit special in that the
     * GridOperator uses them together -- each time a method on the local
     * operator should be called, it will first call the \c alpha_*() method
     * and then call the \c lambda_*() method.
     *
     * If the controlling flag for a method is false, the call to the method
     * is omitted in such a way that the method does not even have to be
     * present on the local operator.  If both the \c do*Skeleton and \c
     * do*Boundary flags are false, the iteration through the intersections is
     * skipped.
     */
    class LocalOperatorInterface
      : public LocalOperatorDefaultFlags
    {
    public:

      static const bool doPatternVolume = true;
      static const bool doPatternVolumePostSkeleton = true;
      static const bool doPatternSkeleton = true;
      static const bool doPatternBoundary = true;
      static const bool doAlphaVolume = true;
      static const bool doAlphaVolumePostSkeleton = true;
      static const bool doAlphaSkeleton = true;
      static const bool doAlphaBoundary = true;
      static const bool doLambdaVolume = true;
      static const bool doLambdaVolumePostSkeleton = true;
      static const bool doLambdaSkeleton = true;
      static const bool doLambdaBoundary = true;

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
      {}

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
      {}

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
      {}

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
      {}

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
      {}

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
      {}

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
      {}

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
      {}

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
       * \note The method should not clear \c r; it should just add its
       *       entries to it.
       *
       * This method is controlled by the flag \ref doLambdaVolume.  For a
       * given element, it is called *before* the lambda_skeleton() and/or
       * lambda_boundary() methods are called (if they are called at all).
       */
      template<typename EG, typename LFSV, typename R>
      void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
      {}

      //! \brief get an element's contribution to lambda after the
      //!        intersections have been handled
      /**
       *
       * \param eg   ElementGeometry describing the entity.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param r    Local part of the residual.
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
      void lambda_volume_post_skeleton(const EG& eg, const LFSV& lfsv, R& r) const
      {}

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
      {}

      //! get a boundary intersections's contribution to lambda
      /**
       * \param ig     IntersectionGeometry describing the intersection.
       *               inside entity.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param r_s    Local part of the residual in the inside entity.
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
      {}

      //! \} Methods for the residual -- constant parts

      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the application of the jacobian for linear problems
      //! \{
      //

      //! Applies an element's jacobian to a vector for a linear problem.
      /**
       * This method performs the following update:
       *
       * \f[ y = y + (\nabla_u \alpha_{\text{vol}}) z, \f]
       *
       * where \f$u\f$ is the solution and the residual is considered linear in \f$u\f$.
       *
       * \param eg   ElementGeometry describing the entity.
       * \param lfsu LocalFunctionSpace of the trial GridFunctionSpace.
       * \param z    Local position in the trial space to which to apply the Jacobian.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param y    Where to store the result.
       *
       * \note This is different from alpha_volume(), since the result will be
       *       linear in \c z, whereas alpha_volume() may include
       *       contributions to the the residual which are constant in \c z.
       *
       * \note The method should not clear \c y; it should just add its
       *       entries to it.
       *
       * \note This method assumes that the residual is linear and the Jacobian
       *       is constant; as such, the method isn't passed a linearization point
       *       at which to evaluate the Jacobian, but only the vector \c z to which
       *       the Jacobian should be applied.
       *
       * \note Your operator should export the flag isLinear=true in this case. This
       *       is done by default.
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
        const LFSU& lfsu, const X& z, const LFSV& lfsv,
        Y& y) const
      {}

      //! Applies an element's jacobian to a vector for a linear problem after the
      //! intersections have been handled.
      /**
       * This method performs the following update:
       *
       * \f[ y = y + (\nabla_u \alpha_{\text{vol\_post\_skel}}) z, \f]
       *
       * where \f$u\f$ is the solution and the residual is considered linear in \f$u\f$.
       *
       * \param eg   ElementGeometry describing the entity.
       * \param lfsu LocalFunctionSpace of the trial GridFunctionSpace.
       * \param z    Local position in the trial space to which to apply the Jacobian.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param y    Where to store the result.
       *
       * \note This is different from alpha_volume_post_skeleton(), since the result
       *       will be linear in \c z, whereas alpha_volume_post_skeleton() may include
       *       contributions to the the residual which are constant in \c z.
       *
       * \note The method should not clear \c y; it should just add its
       *       entries to it.
       *
       * \note This method assumes that the residual is linear and the Jacobian
       *       is constant; as such, the method isn't passed a linearization point
       *       at which to evaluate the Jacobian, but only the vector \c z to which
       *       the Jacobian should be applied.
       *
       * This method is controlled by the flag \ref doAlphaVolumePostSkeleton.  For a
       * given element, it is called *after* the jacobian_apply_skeleton()
       * and/or jacobian_apply_boundary() methods are called (if they are
       * called at all).
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& z, const LFSV& lfsv,
        Y& y) const
      {}

      //! Applies an internal intersections's jacobians to the appropriate vectors
      //! for a linear problem.
      /**
       * This method performs the following updates:
       * \f{align*}{
       *     y_s & = y_s + (\nabla_{u_s} \alpha_{\text{skel},s}) z_s + (\nabla_{u_n} \alpha_{\text{skel},s}) z_n, \\
       *     y_n & = y_n + (\nabla_{u_s} \alpha_{\text{skel},n}) z_s + (\nabla_{u_n} \alpha_{\text{skel},n}) z_n,
       * \f}
       * where \f$u_s\f$ and \f$u_n\f$ are the solution on the inside and outside of the intersection
       *  and the residual is considered linear in \f$u_{\{s,n\}}\f$.
       * \param ig     IntersectionGeometry describing the intersection.
       * \param lfsu_s LocalFunctionSpace of the trial GridFunctionSpace in
       *               the inside entity.
       * \param z_s    Local position in the trial space for the inside entity
       *               to which to apply the Jacobian.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param lfsu_n LocalFunctionSpace of the trial GridFunctionSpace in
       *               the outside entity.
       * \param z_n    Local position in the trial space for the outside entity
       *               to which to apply the Jacobian.
       * \param lfsv_n LocalFunctionSpace of the test GridFunctionSpace in the
       *               outside entity.
       * \param y_s    Where to store the inside entity's result.
       * \param y_n    Where to store the outside entity's result.
       *
       * \note This is different from alpha_skeleton(), since the result will
       *       be linear in \c z_s and \c z_n, whereas alpha_skeleton() may
       *       include contributions to the the residual which are constant in
       *       \c z_s and \c z_n.
       *
       * \note The method should not clear \c y_s and \c y_n; it should just
       *       add its entries to them.
       *
       * \note This method assumes that the residual is linear and the Jacobian
       *       is constant; as such, the method isn't passed linearization points
       *       at which to evaluate the Jacobian, but only the vectors \c z_s and
       *       \c z_n to which the Jacobian should be applied.
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
        const LFSU& lfsu_s, const X& z_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& z_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n) const
      {}

      //! apply a boundary intersections's jacobian for a linear problem.
      /**
       * This method performs the following update:
       *
       * \f[ y = y + (\nabla_u \alpha_{\text{bnd}}) z, \f]
       *
       * where \f$u\f$ is the solution and the residual is considered linear in \f$u\f$.
       *
       * \param ig     IntersectionGeometry describing the intersection.
       * \param lfsu_s LocalFunctionSpace of the trial GridFunctionSpace in
       *               the inside entity.
       * \param z_s    Local position in the trial space for the inside entity
       *               to which to apply the Jacobian.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param y_s    Local part of the residual in the inside entity.
       *
       * \note This is different from alpha_boundary(), since the result will
       *       be linear in \c x, whereas alpha_boundary() may include
       *       contributions to the the residual which are constant in \c x.
       *
       * \note The method should not clear \c y_s; it should just add its
       *       entries to it.
       *
       * \note This method assumes that the residual is linear and the Jacobian
       *       is constant; as such, the method isn't passed a linearization point
       *       at which to evaluate the Jacobian, but only the vector \c z_s
       *       to which the Jacobian should be applied.
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
        const LFSU& lfsu_s, const X& z_s, const LFSV& lfsv_s,
        Y& y_s) const
      {}

      //! \} Methods for the application of the jacobian for linear problems


      //////////////////////////////////////////////////////////////////////
      //
      //! \name Methods for the application of the jacobian for nonlinear problems
      //! \{
      //

      //! Applies an element's jacobian to a vector for a nonlinear problem.
      /**
       * This method performs the following update:
       *
       * \f[ y = y + (\nabla_u \alpha_{\text{vol}})\Big|_x z, \f]
       *
       * where \f$u\f$ is the solution and the Jacobian is evaluated at \f$u=x\f$.
       *
       * \param eg   ElementGeometry describing the entity.
       * \param lfsu LocalFunctionSpace of the trial GridFunctionSpace.
       * \param x    Local position in the trial space at which to evaluate the Jacobian.
       * \param z    Local position in the trial space to which to apply the Jacobian.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param y    Where to store the result.
       *
       * \note This is different from alpha_volume(), since the result will be
       *       linear in \c z, whereas alpha_volume() may include
       *       contributions to the the residual which are constant in \c z.
       *
       * \note The method should not clear \c y; it should just add its
       *       entries to it.
       *
       * \note This method assumes that the residual is not linear and the Jacobian
       *       is not constant; the Jacobian should be evaluated for the vector \c x
       *       and applied to the vector \c z.
       *
       * \note Your linear operator should export the flag isLinear=false in this
       *       case, see LocalOperatorDefaultFlags.
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
        const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv,
        Y& y) const
      {}

      //! Applies an element's jacobian to a vector for a nonlinear problem after the
      //! intersections have been handled.
      /**
       * This method performs the following update:
       *
       * \f[ y = y + (\nabla_u \alpha_{\text{vol\_post\_skel}})\Big|_x z, \f]
       *
       * where \f$u\f$ is the solution and the Jacobian is evaluated at \f$u=x\f$.
       *
       * \param eg   ElementGeometry describing the entity.
       * \param lfsu LocalFunctionSpace of the trial GridFunctionSpace.
       * \param x    Local position in the trial space at which to evaluate the Jacobian.
       * \param z    Local position in the trial space to which to apply the Jacobian.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param y    Where to store the result.
       *
       * \note This is different from alpha_volume_post_skeleton(), since the result
       *       will be linear in \c z, whereas alpha_volume_post_skeleton() may include
       *       contributions to the the residual which are constant in \c z.
       *
       * \note The method should not clear \c y; it should just add its
       *       entries to it.
       *
       * \note This method assumes that the residual is not linear and the Jacobian
       *       is not constant; the Jacobian should be evaluated for the vector \c x
       *       and applied to the vector \c z.
       *
       * This method is controlled by the flag \ref doAlphaVolumePostSkeleton.  For a
       * given element, it is called *after* the jacobian_apply_skeleton()
       * and/or jacobian_apply_boundary() methods are called (if they are
       * called at all).
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv,
        Y& y) const
      {}

      //! Applies an internal intersections's jacobians to the appropriate vectors
      //! for a nonlinear problem.
      /**
       * This method performs the following updates:
       * \f{align*}{
       *     y_s & = y_s + (\nabla_{u_s} \alpha_{\text{skel},s})\Big|_{(x_s,x_n)} z_s
       *                 + (\nabla_{u_n} \alpha_{\text{skel},s})\Big|_{(x_s,x_n)} z_n, \\
       *     y_n & = y_n + (\nabla_{u_s} \alpha_{\text{skel},n})\Big|_{(x_s,x_n)} z_s
       *                 + (\nabla_{u_n} \alpha_{\text{skel},n})\Big|_{(x_s,x_n)} z_n,
       * \f}
       * where \f$u_s\f$ and \f$u_n\f$ are the solution on the inside and the outside of
       * the intersection and the Jacobian is evaluated at \f$u = (x_s,x_n)\f$$.
       *
       * \param ig     IntersectionGeometry describing the intersection.
       * \param lfsu_s LocalFunctionSpace of the trial GridFunctionSpace in
       *               the inside entity.
       * \param x_s    Local position in the trial space for the inside entity
       *               at which to evaluate the Jacobian.
       * \param z_s    Local position in the trial space for the inside entity
       *               to which to apply the Jacobian.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param lfsu_n LocalFunctionSpace of the trial GridFunctionSpace in
       *               the outside entity.
       * \param x_n    Local position in the trial space for the outside entity
       *               at which to evaluate the Jacobian.
       * \param z_n    Local position in the trial space for the outside entity
       *               to which to apply the Jacobian.
       * \param lfsv_n LocalFunctionSpace of the test GridFunctionSpace in the
       *               outside entity.
       * \param y_s    Where to store the inside entity's result.
       * \param y_n    Where to store the outside entity's result.
       *
       * \note This is different from alpha_skeleton(), since the result will
       *       be linear in \c z_s and \c z_n, whereas alpha_skeleton() may
       *       include contributions to the the residual which are constant in
       *       \c z_s and \c z_n.
       *
       * \note The method should not clear \c y_s and \c y_n; it should just
       *       add its entries to them.
       *
       * \note This method assumes that the residual is not linear and the Jacobian
       *       is not constant; the Jacobian should be evaluated for the vector \c x
       *       and applied to the vector \c z.
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
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const X& z_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n) const
      {}

      //! apply a boundary intersections's jacobian for a nonlinear problem.
      /**
       * This method performs the following update:
       *
       * \f[ y = y + (\nabla_u \alpha_{\text{bnd}})\Big|_{x_s} z, \f]
       *
       * where \f$u\f$ is the solution and the Jacobian is evaluated at \f$u = x_s\f$.
       *
       * \param ig     IntersectionGeometry describing the intersection.
       * \param lfsu_s LocalFunctionSpace of the trial GridFunctionSpace in
       *               the inside entity.
       * \param x_s    Local position in the trial space for the inside entity
       *               at which to evaluate the Jacobian.
       * \param z_s    Local position in the trial space for the inside entity
       *               to which to apply the Jacobian.
       * \param lfsv_s LocalFunctionSpace of the test GridFunctionSpace in the
       *               inside entity.
       * \param y_s    Local part of the residual in the inside entity.
       *
       * \note This is different from alpha_boundary(), since the result will
       *       be linear in \c x, whereas alpha_boundary() may include
       *       contributions to the the residual which are constant in \c x.
       *
       * \note The method should not clear \c y_s; it should just add its
       *       entries to it.
       *
       * \note This method assumes that the residual is not linear and the Jacobian
       *       is not constant; the Jacobian should be evaluated for the vector \c x
       *       and applied to the vector \c z.
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
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        Y& y_s) const
      {}

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
       * \note The method should not clear \c mat; it should just add its
       *       entries to it.
       *
       * This method is controlled by the flag \ref doAlphaVolume.  For a
       * given element, it is called *before* the jacobian_skeleton() and/or
       * jacobian_boundary() methods are called (if they are called at all).
       */
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        LocalMatrix& mat) const
      {}

      //! get an element's jacobian after the intersections have been handled
      /**
       * \param eg   ElementGeometry describing the entity.
       * \param lfsu LocalFunctionSpace of the trial GridFunctionSpace.
       * \param x    Local position in the trial GridFunctionSpace.
       * \param lfsv LocalFunctionSpace of the test GridFunctionSpace.
       * \param mat  Where to store the contribution to the jacobian.
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
               typename LocalMatrix>
      void jacobian_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        LocalMatrix& mat) const
      {}

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
       * \note The method should not clear \c mat_ss, \c mat_sn, \c mat_ns, \c
       *       mat_nn; it should just add its entries to them.
       *
       * This method is controlled by the flag \ref doAlphaSkeleton.  For a
       * given element, it's calls happen intermingled with the calls to
       * jacobian_boundary(), but after the call to jacobian_volume() and
       * before the call to jacobian_volume_post_skeleton().
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        LocalMatrix& mat_ss, LocalMatrix& mat_sn,
        LocalMatrix& mat_ns, LocalMatrix& mat_nn) const
      {}

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
       * \note The method should not clear \c mat_ss; it should just add its
       *       entries to it.
       *
       * This method is controlled by the flag \ref doAlphaBoundary.  For a
       * given element, it's calls happen intermingled with the calls to
       * jacobian_skeleton(), but after the call to jacobian_volume() and
       * before the call to jacobian_volume_post_skeleton().
       */
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename LocalMatrix>
      void jacobian_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        LocalMatrix& mat_ss) const
      {}

      //! \} Methods to extract the jacobian
    };

    //! \} group LocalOperatorDefaultImp
  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_INTERFACE_HH
