// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_NUMERICALRESIDUAL_HH
#define DUNE_PDELAB_LOCALOPERATOR_NUMERICALRESIDUAL_HH

#include <cmath>

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup LocalOperatorDefaultImp
    //! \ingroup LocalOperator
    //! \{

    ////////////////////////////////////////////////////////////////////////
    //
    //  Implementation of alpha_*() in terms of jacobian_*()
    //

    //! Implement alpha_volume() based on jacobian_volume()
    /**
     * Derive from this class to add an alpha_volume() method.  The derived
     * class needs to implement jacobian_volume().
     *
     * \note This only works in the common case that alpha_volume() does not
     *       contain an affine shift (the affine shift can be moved to
     *       lambda_volume() instead).
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class JacobianBasedAlphaVolume
    {
    public:

      //! compute \f$\alpha_\text{vol}\f$
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename R>
      void alpha_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        R& r) const
      {
        typedef LocalMatrix<typename R::value_type> Jacobian;
        typedef typename Jacobian::WeightedAccumulationView JacobianView;

        Jacobian mat(r.size(),x.size(), 0);
        JacobianView matview = mat.weightedAccumulationView(1.0);
        asImp().jacobian_volume(eg, lfsu, x, lfsv, matview);
        // we need to include the weight here, as umv() and usmv() operate on the bare container for effiency
        mat.usmv(r.weight(),x,r);
      }

    private:
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! Implement alpha_skeleton() based on jacobian_skeleton()
    /**
     * Derive from this class to add an alpha_skeleton() method.  The derived
     * class needs to implement jacobian_skeleton().
     *
     * \note This only works in the common case that alpha_skeleton() does not
     *       contain an affine shift (the affine shift can be moved to
     *       lambda_skeleton() instead).
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class JacobianBasedAlphaSkeleton
    {
    public:

      //! compute \f$\alpha_\text{skel}\f$
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename R>
      void alpha_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        R& r_s, R& r_n) const
      {
        typedef LocalMatrix<typename R::value_type> Jacobian;
        typedef typename Jacobian::WeightedAccumulationView JacobianView;

        Jacobian mat_ss(r_s.size(),x_s.size(),0);
        Jacobian mat_sn(r_s.size(),x_n.size(),0);
        Jacobian mat_ns(r_n.size(),x_s.size(),0);
        Jacobian mat_nn(r_n.size(),x_n.size(),0);

        JacobianView view_ss = mat_ss.weightedAccumulationView(1.0);
        JacobianView view_sn = mat_sn.weightedAccumulationView(1.0);
        JacobianView view_ns = mat_ns.weightedAccumulationView(1.0);
        JacobianView view_nn = mat_nn.weightedAccumulationView(1.0);

        asImp().jacobian_skeleton(ig,
          lfsu_s, x_s, lfsv_s,
          lfsu_n, x_n, lfsv_n,
          view_ss, view_sn, view_ns, view_nn);
        // TODO: Reihenfolge der Multiplikationen!
        mat_ss.usmv(r_s.weight(),x_s,r_s);
        mat_ns.usmv(r_n.weight(),x_s,r_n);
        mat_sn.usmv(r_s.weight(),x_n,r_s);
        mat_nn.usmv(r_n.weight(),x_n,r_n);
      }

    private:
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! Implement alpha_boundary() based on jacobian_boundary()
    /**
     * Derive from this class to add an alpha_boundary() method.  The derived
     * class needs to implement jacobian_boundary().
     *
     * \note This only works in the common case that alpha_boundary() does not
     *       contain an affine shift (the affine shift can be moved to
     *       lambda_boundary() instead).
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class JacobianBasedAlphaBoundary
    {
    public:

      //! compute \f$\alpha_\text{bnd}\f$
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename R>
      void alpha_boundary
      ( const IG& ig,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        R& r) const
      {
        typedef LocalMatrix<typename R::value_type> Jacobian;
        typedef typename Jacobian::WeightedAccumulationView JacobianView;

        Jacobian mat(x.size(),r.size(), 0);
        JacobianView view = mat.weightedAccumulationView(1.0);
        asImp().jacobian_boundary(ig, lfsu, x, lfsv, view);
        // we need to include the weight here, as umv() and usmv() operate on the bare container for effiency
        mat.usmv(r.weight(),x,r);
      }

    private:
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! \} group LocalOperatorDefaultImp
  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_NUMERICALRESIDUAL_HH
