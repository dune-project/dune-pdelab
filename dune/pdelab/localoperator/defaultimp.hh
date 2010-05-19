// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_DEFAULTIMP_HH
#define DUNE_PDELAB_LOCALOPERATOR_DEFAULTIMP_HH

#include <cmath>
#include <vector>

#include <dune/pdelab/gridoperatorspace/localmatrix.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperatorDefaultImp
    //! \ingroup LocalOperator
    //! \{

    ////////////////////////////////////////////////////////////////////////
    //
    //  Numerical implementation of jacobian_*() in terms of alpha_*()
    //

    //! Implement jacobian_volume() based on alpha_volume()
    /**
     * Derive from this class to add numerical jacobian for volume.  The
     * derived class needs to implement alpha_volume().
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalJacobianVolume
    {
    public:

      //! compute local jacobian of the volume term
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename R>
      void jacobian_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        LocalMatrix<R>& mat) const
      {
        const R epsilon=1E-11; // problem: this depends on data type R!
        const int m=lfsv.size();
        const int n=lfsu.size();

        X u(x);
        std::vector<R> down(m,0.0),up(m);

        asImp().alpha_volume(eg,lfsu,u,lfsv,down);
        for (int j=0; j<n; j++) // loop over columns
        {
          for (int k=0; k<m; k++) up[k]=0.0;
          R delta = epsilon*(1.0+std::abs(u[j]));
          u[j] += delta;
          asImp().alpha_volume(eg,lfsu,u,lfsv,up);
          for (int i=0; i<m; i++)
            mat(i,j) += (up[i]-down[i])/delta;
          u[j] = x[j];
        }
      }

    private:
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! Implement jacobian_volume_post_skeleton() based on
    //! alpha_volume_post_skeleton()
    /**
     * Derive from this class to add numerical jacobian for volume
     * (post-skeleton part).  The derived class needs to implement
     * alpha_volume_post_skeleton().
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalJacobianVolumePostSkeleton
    {
    public:

      //! compute local post-skeleton jacobian of the volume term
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename R>
      void jacobian_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        LocalMatrix<R>& mat) const
      {
        const R epsilon=1E-11; // problem: this depends on data type R!
        const int m=lfsv.size();
        const int n=lfsu.size();

        X u(x);
        std::vector<R> down(m,0.0),up(m);

        asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,down);
        for (int j=0; j<n; j++) // loop over columns
        {
          for (int k=0; k<m; k++) up[k]=0.0;
          R delta = epsilon*(1.0+std::abs(u[j]));
          u[j] += delta;
          asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,up);
          for (int i=0; i<m; i++)
            mat(i,j) += (up[i]-down[i])/delta;
          u[j] = x[j];
        }
      }

    private:
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! Implement jacobian_skeleton() based on alpha_skeleton()
    /**
     * Derive from this class to add numerical jacobian for skeleton.  The
     * derived class needs to implement alpha_skeleton().
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalJacobianSkeleton
    {
    public:

      //! compute local jacobian of the skeleton term
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename R>
      void jacobian_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        LocalMatrix<R>& mat_ss, LocalMatrix<R>& mat_sn,
        LocalMatrix<R>& mat_ns, LocalMatrix<R>& mat_nn) const
      {
        const R epsilon=1E-11; // problem: this depends on data type R!
        const int m_s=lfsv_s.size();
        const int m_n=lfsv_n.size();
        const int n_s=lfsu_s.size();
        const int n_n=lfsu_n.size();

        X u_s(x_s);
        X u_n(x_n);
        std::vector<R> down_s(m_s,0.0),up_s(m_s);
        std::vector<R> down_n(m_n,0.0),up_n(m_n);

        // base line
        asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,down_s,
                               down_n);

        // jiggle in self
        for (int j=0; j<n_s; j++)
        {
          for (int k=0; k<m_s; k++) up_s[k]=0.0;
          for (int k=0; k<m_n; k++) up_n[k]=0.0;
          R delta = epsilon*(1.0+std::abs(u_s[j]));
          u_s[j] += delta;
          asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,up_s,
                                 up_n);
          for (int i=0; i<m_s; i++)
            mat_ss(i,j) += (up_s[i]-down_s[i])/delta;
          for (int i=0; i<m_n; i++)
            mat_ns(i,j) += (up_n[i]-down_n[i])/delta;
          u_s[j] = x_s[j];
        }

        // jiggle in neighbor
        for (int j=0; j<n_n; j++)
        {
          for (int k=0; k<m_s; k++) up_s[k]=0.0;
          for (int k=0; k<m_n; k++) up_n[k]=0.0;
          R delta = epsilon*(1.0+std::abs(u_s[j]));
          u_n[j] += delta;
          asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,up_s,
                                 up_n);
          for (int i=0; i<m_s; i++)
            mat_sn(i,j) += (up_s[i]-down_s[i])/delta;
          for (int i=0; i<m_n; i++)
            mat_nn(i,j) += (up_n[i]-down_n[i])/delta;
          u_n[j] = x_n[j];
        }
      }

    private:
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! Implement jacobian_boundary() based on alpha_boundary()
    /**
     * Derive from this class to add numerical jacobian for boundary.  The
     * derived class needs to implement alpha_boundary().
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalJacobianBoundary
    {
    public:

      //! compute local jacobian of the boundary term
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename R>
      void jacobian_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        LocalMatrix<R>& mat_ss) const
      {
        const R epsilon=1E-11; // problem: this depends on data type R!
        const int m_s=lfsv_s.size();
        const int n_s=lfsu_s.size();

        X u_s(x_s);
        std::vector<R> down_s(m_s,0.0),up_s(m_s);

        // base line
        asImp().alpha_boundary(ig,lfsu_s,u_s,lfsv_s,down_s);

        // jiggle in self
        for (int j=0; j<n_s; j++)
        {
          for (int k=0; k<m_s; k++) up_s[k]=0.0;
          R delta = epsilon*(1.0+std::abs(u_s[j]));
          u_s[j] += delta;
          asImp().alpha_boundary(ig,lfsu_s,u_s,lfsv_s,up_s);
          for (int i=0; i<m_s; i++)
            mat_ss(i,j) += (up_s[i]-down_s[i])/delta;
          u_s[j] = x_s[j];
        }
      }

    private:
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    ////////////////////////////////////////////////////////////////////////
    //
    //  Numerical implementation of jacobian_apply_*() in terms of alpha_*()
    //

    //! Implement jacobian_apply_volume() based on alpha_volume()
    /**
     * Derive from this class to add numerical jacobian application for
     * volume.  The derived class needs to implement alpha_volume().
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalJacobianApplyVolume
    {
    public:

      //! apply local jacobian of the volume term
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        Y& y) const
      {
        typedef typename X::value_type R;
        const R epsilon=1E-11; // problem: this depends on data type R!
        const int m=lfsv.size();
        const int n=lfsu.size();

        X u(x);
        std::vector<R> down(m,0.0),up(m);

        asImp().alpha_volume(eg,lfsu,u,lfsv,down);
        for (int j=0; j<n; j++) // loop over columns
        {
          for (int k=0; k<m; k++) up[k]=0.0;
          R delta = epsilon*(1.0+std::abs(u[j]));
          u[j] += delta;
          asImp().alpha_volume(eg,lfsu,u,lfsv,up);
          for (int i=0; i<m; i++)
            y[i] += ((up[i]-down[i])/delta)*x[j];
          u[j] = x[j];
        }
      }

    private:
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! Implement jacobian_apply_volume_post_skeleton() based on
    //! alpha_volume_post_skeleton()
    /**
     * Derive from this class to add numerical jacobian application for volume
     * (post skeleton part).  The derived class needs to implement
     * alpha_volume().
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalJacobianApplyVolumePostSkeleton
    {
    public:

      //! apply local jacobian of the volume term (post skeleton part)
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        Y& y) const
      {
        typedef typename X::value_type R;
        const R epsilon=1E-11; // problem: this depends on data type R!
        const int m=lfsv.size();
        const int n=lfsu.size();

        X u(x);
        std::vector<R> down(m,0.0),up(m);

        asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,down);
        for (int j=0; j<n; j++) // loop over columns
        {
          for (int k=0; k<m; k++) up[k]=0.0;
          R delta = epsilon*(1.0+std::abs(u[j]));
          u[j] += delta;
          asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,up);
          for (int i=0; i<m; i++)
            y[i] += ((up[i]-down[i])/delta)*x[j];
          u[j] = x[j];
        }
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };

    //! Implement jacobian_apply_skeleton() based on alpha_skeleton()
    /**
     * Derive from this class to add numerical jacobian application for
     * skeleton.  The derived class needs to implement alpha_skeleton().
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalJacobianApplySkeleton
    {
    public:

      //! apply local jacobian of the skeleton term
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n) const
      {
        typedef typename X::value_type R;
        const R epsilon=1E-11; // problem: this depends on data type R!
        const int m_s=lfsv_s.size();
        const int m_n=lfsv_n.size();
        const int n_s=lfsu_s.size();
        const int n_n=lfsu_n.size();

        X u_s(x_s);
        X u_n(x_n);
        std::vector<R> down_s(m_s,0.0),up_s(m_s);
        std::vector<R> down_n(m_n,0.0),up_n(m_n);

        // base line
        asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,down_s,
                               down_n);

        // jiggle in self
        for (int j=0; j<n_s; j++)
        {
          for (int k=0; k<m_s; k++) up_s[k]=0.0;
          for (int k=0; k<m_n; k++) up_n[k]=0.0;
          R delta = epsilon*(1.0+std::abs(u_s[j]));
          u_s[j] += delta;
          asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,up_s,
                                 up_n);
          for (int i=0; i<m_s; i++)
            y_s[i] += ((up_s[i]-down_s[i])/delta)*x_s[j];
          for (int i=0; i<m_n; i++)
            y_n[i] += ((up_n[i]-down_n[i])/delta)*x_s[j];
          u_s[j] = x_s[j];
        }

        // jiggle in neighbor
        for (int j=0; j<n_n; j++)
        {
          for (int k=0; k<m_s; k++) up_s[k]=0.0;
          for (int k=0; k<m_n; k++) up_n[k]=0.0;
          R delta = epsilon*(1.0+std::abs(u_s[j]));
          u_n[j] += delta;
          asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,up_s,
                                 up_n);
          for (int i=0; i<m_s; i++)
            y_s[i] += ((up_s[i]-down_s[i])/delta)*x_n[j];
          for (int i=0; i<m_n; i++)
            y_n[i] += ((up_n[i]-down_n[i])/delta)*x_n[j];
          u_n[j] = x_n[j];
        }
      }

    private:
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! Implement jacobian_apply_boundary() based on alpha_boundary()
    /**
     * Derive from this class to add numerical jacobian application for
     * boundary.  The derived class needs to implement alpha_boundary().
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalJacobianApplyBoundary
    {
    public:

      //! apply local jacobian of the boundaryterm
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        Y& y_s) const
      {
        typedef typename X::value_type R;
        const R epsilon=1E-11; // problem: this depends on data type R!
        const int m_s=lfsv_s.size();
        const int n_s=lfsu_s.size();

        X u_s(x_s);
        std::vector<R> down_s(m_s,0.0),up_s(m_s);

        // base line
        asImp().alpha_boundary(ig,lfsu_s,u_s,lfsv_s,down_s);

        // jiggle in self
        for (int j=0; j<n_s; j++)
        {
          for (int k=0; k<m_s; k++) up_s[k]=0.0;
          R delta = epsilon*(1.0+std::abs(u_s[j]));
          u_s[j] += delta;
          asImp().alpha_boundary(ig,lfsu_s,u_s,lfsv_s,up_s);
          for (int i=0; i<m_s; i++)
            y_s[i] += ((up_s[i]-down_s[i])/delta)*x_s[j];
          u_s[j] = x_s[j];
        }
      }

    private:
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

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
        LocalMatrix<typename R::value_type> mat(lfsu.size(),lfsu.size());
        asImp().jacobian_volume(eg, lfsu, x, lfsv, mat);
        mat.umv(x,r);
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
        LocalMatrix<typename R::value_type> mat_ss(lfsu_s.size(),
                                                   lfsu_s.size());
        LocalMatrix<typename R::value_type> mat_sn(lfsu_s.size(),
                                                   lfsu_n.size());
        LocalMatrix<typename R::value_type> mat_ns(lfsu_n.size(),
                                                   lfsu_s.size());
        LocalMatrix<typename R::value_type> mat_nn(lfsu_n.size(),
                                                   lfsu_n.size());
        asImp().jacobian_skeleton(ig,
          lfsu_s, x_s, lfsv_s,
          lfsu_n, x_n, lfsv_n,
          mat_ss, mat_sn, mat_ns, mat_nn);
        // TODO: Reihenfolge der Multiplikationen!
        mat_ss.umv(x_s,r_s);
        mat_ns.umv(x_n,r_s);
        mat_sn.umv(x_s,r_n);
        mat_nn.umv(x_n,r_n);
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
        LocalMatrix<typename R::value_type> mat(lfsu.size(),lfsu.size());
        asImp().jacobian_boundary(ig, lfsu, x, lfsv, mat);
        mat.umv(x,r);
      }

    private:
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! \} group LocalOperatorDefaultImp
  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_DEFAULTIMP_HH
