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
      NumericalJacobianVolume ()
        : epsilon(1e-7)
      {}

      NumericalJacobianVolume (double epsilon_)
        : epsilon(epsilon_)
      {}

      //! compute local jacobian of the volume term
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename Jacobian>
      void jacobian_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        Jacobian& mat) const
      {
        typedef typename X::value_type D;
        typedef typename Jacobian::value_type R;
        typedef LocalVector<R,TestSpaceTag,typename Jacobian::weight_type> ResidualVector;
        typedef typename ResidualVector::WeightedAccumulationView ResidualView;

        const int m=lfsv.size();
        const int n=lfsu.size();

        X u(x);

        ResidualVector down(m),up(m);
        ResidualView downview = down.weightedAccumulationView(mat.weight());
        ResidualView upview = up.weightedAccumulationView(mat.weight());


        asImp().alpha_volume(eg,lfsu,u,lfsv,downview);
        for (int j=0; j<n; j++) // loop over columns
        {
          for (int k=0; k<mat.nrows(); k++) up[k]=0.0;
          D delta = epsilon*(1.0+std::abs(u[lfsu.localIndex(j)]));
          u[lfsu.localIndex(j)] += delta;
          asImp().alpha_volume(eg,lfsu,u,lfsv,upview);
          for (int i=0; i<m; i++)
            mat.rawAccumulate(lfsv.localIndex(i),lfsu.localIndex(j),
                              (up[lfsv.localIndex(i)]-down[lfsv.localIndex(i)])/delta);
          u[lfsu.localIndex(j)] = x[lfsu.localIndex(j)];
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
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
      NumericalJacobianVolumePostSkeleton ()
        : epsilon(1e-7)
      {}

      NumericalJacobianVolumePostSkeleton (double epsilon_)
        : epsilon(epsilon_)
      {}

      //! compute local post-skeleton jacobian of the volume term
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename Jacobian>
      void jacobian_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        Jacobian& mat) const
      {
        typedef typename X::value_type D;
        typedef typename Jacobian::value_type R;
        typedef LocalVector<R,TestSpaceTag,typename Jacobian::weight_type> ResidualVector;
        typedef typename ResidualVector::WeightedAccumulationView ResidualView;

        const int m=lfsv.size();
        const int n=lfsu.size();

        X u(x);

        ResidualVector down(m),up(m);
        ResidualView downview = down.weightedAccumulationView(mat.weight());
        ResidualView upview = up.weightedAccumulationView(mat.weight());

        asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,downview);
        for (int j=0; j<n; j++) // loop over columns
        {
          for (int k=0; k<mat.nrows(); k++) up[k]=0.0;
          D delta = epsilon*(1.0+std::abs(u[lfsu.localIndex(j)]));
          u[lfsu.localIndex(j)] += delta;
          asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,upview);
          for (int i=0; i<m; i++)
            mat.rawAccumulate(lfsv.localIndex(i),lfsu.localIndex(j),
                              (up[lfsv.localIndex(i)]-down[lfsv.localIndex(i)])/delta);
          u[lfsu.localIndex(j)] = x[lfsu.localIndex(j)];
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
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
      NumericalJacobianSkeleton ()
        : epsilon(1e-7)
      {}

      NumericalJacobianSkeleton (double epsilon_)
        : epsilon(epsilon_)
      {}

      //! compute local jacobian of the skeleton term
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename Jacobian>
      void jacobian_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        Jacobian& mat_ss, Jacobian& mat_sn,
        Jacobian& mat_ns, Jacobian& mat_nn) const
      {
        typedef typename X::value_type D;
        typedef typename Jacobian::value_type R;
        typedef LocalVector<R,TestSpaceTag,typename Jacobian::weight_type> ResidualVector;
        typedef typename ResidualVector::WeightedAccumulationView ResidualView;

        const int m_s=lfsv_s.size();
        const int m_n=lfsv_n.size();
        const int n_s=lfsu_s.size();
        const int n_n=lfsu_n.size();

        X u_s(x_s);
        X u_n(x_n);

        ResidualVector down_s(m_s),up_s(m_s);
        ResidualView downview_s = down_s.weightedAccumulationView(1.0);
        ResidualView upview_s = up_s.weightedAccumulationView(1.0);

        ResidualVector down_n(m_n),up_n(m_n);
        ResidualView downview_n = down_n.weightedAccumulationView(1.0);
        ResidualView upview_n = up_n.weightedAccumulationView(1.0);

        // base line
        asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,downview_s,
                               downview_n);

        // jiggle in self
        for (int j=0; j<n_s; j++)
        {
          for (int k=0; k<mat_ss.nrows(); k++) up_s[k]=0.0;
          for (int k=0; k<mat_nn.nrows(); k++) up_n[k]=0.0;
          D delta = epsilon*(1.0+std::abs(u_s[lfsu_s.localIndex(j)]));
          u_s[lfsu_s.localIndex(j)] += delta;
          asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,upview_s,
                                 upview_n);
          for (int i=0; i<m_s; i++)
            mat_ss.accumulate(lfsv_s.localIndex(i),lfsu_s.localIndex(j),
                              (up_s[lfsv_s.localIndex(i)]-down_s[lfsv_s.localIndex(i)])/delta);
          for (int i=0; i<m_n; i++)
            mat_ns.accumulate(lfsv_n.localIndex(i),lfsu_s.localIndex(j),
                              (up_n[lfsv_n.localIndex(i)]-down_n[lfsv_n.localIndex(i)])/delta);
          u_s[lfsu_s.localIndex(j)] = x_s[lfsu_s.localIndex(j)];
        }

        // jiggle in neighbor
        for (int j=0; j<n_n; j++)
        {
          for (int k=0; k<mat_ss.nrows(); k++) up_s[k]=0.0;
          for (int k=0; k<mat_nn.nrows(); k++) up_n[k]=0.0;
          D delta = epsilon*(1.0+std::abs(u_n[lfsu_n.localIndex(j)]));
          u_n[lfsu_n.localIndex(j)] += delta;
          asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,upview_s,
                                 upview_n);
          for (int i=0; i<m_s; i++)
            mat_sn.accumulate(lfsv_s.localIndex(i),lfsu_n.localIndex(j),
                              (up_s[lfsv_s.localIndex(i)]-down_s[lfsv_s.localIndex(i)])/delta);
          for (int i=0; i<m_n; i++)
            mat_nn.accumulate(lfsv_n.localIndex(i),lfsu_n.localIndex(j),
                              (up_n[lfsv_n.localIndex(i)]-down_n[lfsv_n.localIndex(i)])/delta);
          u_n[lfsu_n.localIndex(j)] = x_n[lfsu_n.localIndex(j)];
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
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
      NumericalJacobianBoundary ()
        : epsilon(1e-7)
      {}

      NumericalJacobianBoundary (double epsilon_)
        : epsilon(epsilon_)
      {}

      //! compute local jacobian of the boundary term
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename Jacobian>
      void jacobian_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        Jacobian& mat_ss) const
      {
        typedef typename X::value_type D;
        typedef typename Jacobian::value_type R;
        typedef LocalVector<R,TestSpaceTag,typename Jacobian::weight_type> ResidualVector;
        typedef typename ResidualVector::WeightedAccumulationView ResidualView;

        const int m_s=lfsv_s.size();
        const int n_s=lfsu_s.size();

        X u_s(x_s);

        ResidualVector down_s(m_s),up_s(m_s);
        ResidualView downview_s = down_s.weightedAccumulationView(mat_ss.weight());
        ResidualView upview_s = up_s.weightedAccumulationView(mat_ss.weight());;


        // base line
        asImp().alpha_boundary(ig,lfsu_s,u_s,lfsv_s,downview_s);

        // jiggle in self
        for (int j=0; j<n_s; j++)
        {
          for (int k=0; k<mat_ss.nrows(); k++) up_s[k]=0.0;
          D delta = epsilon*(1.0+std::abs(u_s[lfsu_s.localIndex(j)]));
          u_s[lfsu_s.localIndex(j)] += delta;
          asImp().alpha_boundary(ig,lfsu_s,u_s,lfsv_s,upview_s);
          for (int i=0; i<m_s; i++)
            mat_ss.rawAccumulate(lfsv_s.localIndex(i),lfsu_s.localIndex(j),
                                 (up_s[lfsv_s.localIndex(i)]-down_s[lfsv_s.localIndex(i)])/delta);
          u_s[lfsu_s.localIndex(j)] = x_s[lfsu_s.localIndex(j)];
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
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
      NumericalJacobianApplyVolume ()
        : epsilon(1e-7)
      {}

      NumericalJacobianApplyVolume (double epsilon_)
        : epsilon(epsilon_)
      {}

      //! apply local jacobian of the volume term
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_volume
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        Y& y) const
      {
        typedef typename X::value_type D;
        typedef typename X::size_type size_type;
        typedef typename Y::value_type R;
        typedef LocalVector<R,TestSpaceTag,typename Y::weight_type> ResidualVector;
        typedef typename ResidualVector::WeightedAccumulationView ResidualView;

        const int m=lfsv.size();
        const int n=lfsu.size();

        X u(x);

        ResidualVector down(m),up(m);
        ResidualView downview = down.weightedAccumulationView(y.weight());
        ResidualView upview = up.weightedAccumulationView(y.weight());

        asImp().alpha_volume(eg,lfsu,u,lfsv,downview);
        for (int j=0; j<n; j++) // loop over columns
        {
          for (size_type k=0; k<y.size(); k++) up[k]=0.0;
          D delta = epsilon*(1.0+std::abs(u[lfsu.localIndex(j)]));
          u[lfsu.localIndex(j)] += delta;
          asImp().alpha_volume(eg,lfsu,u,lfsv,upview);
          for (int i=0; i<m; i++)
            y.rawAccumulate(lfsv.localIndex(i),
                            ((up[lfsv.localIndex(i)]-down[lfsv.localIndex(i)])/delta)*x[lfsu.localIndex(j)]);
          u[lfsu.localIndex(j)] = x[lfsu.localIndex(j)];
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
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
      NumericalJacobianApplyVolumePostSkeleton ()
        : epsilon(1e-7)
      {}

      NumericalJacobianApplyVolumePostSkeleton (double epsilon_)
        : epsilon(epsilon_)
      {}

      //! apply local jacobian of the volume term (post skeleton part)
      template<typename EG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_volume_post_skeleton
      ( const EG& eg,
        const LFSU& lfsu, const X& x, const LFSV& lfsv,
        Y& y) const
      {
        typedef typename X::value_type D;
        typedef typename X::size_type size_type;
        typedef typename Y::value_type R;
        typedef LocalVector<R,TestSpaceTag,typename Y::weight_type> ResidualVector;
        typedef typename ResidualVector::WeightedAccumulationView ResidualView;

        const int m=lfsv.size();
        const int n=lfsu.size();

        X u(x);

        ResidualVector down(m),up(m);
        ResidualView downview = down.weightedAccumulationView(y.weight());
        ResidualView upview = up.weightedAccumulationView(y.weight());

        asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,downview);
        for (int j=0; j<n; j++) // loop over columns
        {
          for (size_type k=0; k<y.size(); k++) up[k]=0.0;
          D delta = epsilon*(1.0+std::abs(u[lfsu.localIndex(j)]));
          u[lfsu.localIndex(j)] += delta;
          asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,upview);
          for (int i=0; i<m; i++)
            y.rawAccumulate(lfsv.localIndex(i),
                            ((up[lfsv.localIndex(i)]-down[lfsv.localIndex(i)])/delta)*x[lfsu.localIndex(j)]);
          u[lfsu.localIndex(j)] = x[lfsu.localIndex(j)];
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
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
      NumericalJacobianApplySkeleton ()
        : epsilon(1e-7)
      {}

      NumericalJacobianApplySkeleton (double epsilon_)
        : epsilon(epsilon_)
      {}

      //! apply local jacobian of the skeleton term
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_skeleton
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n) const
      {
        typedef typename X::value_type D;
        typedef typename Y::value_type R;
        typedef typename X::size_type size_type;
        typedef LocalVector<R,TestSpaceTag,typename Y::weight_type> ResidualVector;
        typedef typename ResidualVector::WeightedAccumulationView ResidualView;

        const int m_s=lfsv_s.size();
        const int m_n=lfsv_n.size();
        const int n_s=lfsu_s.size();
        const int n_n=lfsu_n.size();

        X u_s(x_s);
        X u_n(x_n);

        ResidualVector down_s(m_s),up_s(m_s);
        ResidualView downview_s = down_s.weightedAccumulationView(1.0);
        ResidualView upview_s = up_s.weightedAccumulationView(1.0);

        ResidualVector down_n(m_n),up_n(m_n);
        ResidualView downview_n = down_n.weightedAccumulationView(1.0);
        ResidualView upview_n = up_n.weightedAccumulationView(1.0);

        // base line
        asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,downview_s,
                               downview_n);

        // jiggle in self
        for (int j=0; j<n_s; j++)
        {
          for (size_type k=0; k<y_s.size(); k++) up_s[k]=0.0;
          for (size_type k=0; k<y_n.size(); k++) up_n[k]=0.0;
          D delta = epsilon*(1.0+std::abs(u_s[lfsu_s.localIndex(j)]));
          u_s[lfsu_s.localIndex(j)] += delta;
          asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,upview_s,
                                 upview_n);
          for (int i=0; i<m_s; i++)
            y_s.accumulate(lfsv_s.localIndex(i),
                           ((up_s[lfsv_s.localIndex(i)]-down_s[lfsv_s.localIndex(i)])/delta)*x_s[lfsu_s.localIndex(j)]);
          for (int i=0; i<m_n; i++)
            y_n.accumulate(lfsv_n.localIndex(i),
                           ((up_n[lfsv_n.localIndex(i)]-down_n[lfsv_n.localIndex(i)])/delta)*x_s[lfsu_s.localIndex(j)]);
          u_s[lfsu_s.localIndex(j)] = x_s[lfsu_s.localIndex(j)];
        }

        // jiggle in neighbor
        for (int j=0; j<n_n; j++)
        {
          for (size_type k=0; k<y_s.size(); k++) up_s[k]=0.0;
          for (size_type k=0; k<y_n.size(); k++) up_n[k]=0.0;
          D delta = epsilon*(1.0+std::abs(u_n[lfsu_n.localIndex(j)]));
          u_n[lfsu_n.localIndex(j)] += delta;
          asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,upview_s,
                                 upview_n);
          for (int i=0; i<m_s; i++)
            y_s.accumulate(lfsv_n.localIndex(i),
                           ((up_s[lfsv_s.localIndex(i)]-down_s[lfsv_s.localIndex(i)])/delta)*x_n[lfsu_n.localIndex(j)]);
          for (int i=0; i<m_n; i++)
            y_n.accumulate(lfsv_n.localIndex(i),
                           ((up_n[lfsv_n.localIndex(i)]-down_n[lfsv_n.localIndex(i)])/delta)*x_n[lfsu_n.localIndex(j)]);
          u_n[lfsu_n.localIndex(j)] = x_n[lfsu_n.localIndex(j)];
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
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
      NumericalJacobianApplyBoundary ()
        : epsilon(1e-7)
      {}

      NumericalJacobianApplyBoundary (double epsilon_)
        : epsilon(epsilon_)
      {}

      //! apply local jacobian of the boundaryterm
      template<typename IG, typename LFSU, typename X, typename LFSV,
               typename Y>
      void jacobian_apply_boundary
      ( const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        Y& y_s) const
      {
        typedef typename X::value_type D;
        typedef typename Y::value_type R;
        typedef typename X::size_type size_type;
        typedef LocalVector<R,TestSpaceTag,typename Y::weight_type> ResidualVector;
        typedef typename ResidualVector::WeightedAccumulationView ResidualView;

        const int m_s=lfsv_s.size();
        const int n_s=lfsu_s.size();

        X u_s(x_s);

        ResidualVector down_s(m_s),up_s(m_s);
        ResidualView downview_s = down_s.weightedAccumulationView(1.0);
        ResidualView upview_s = up_s.weightedAccumulationView(1.0);

        // base line
        asImp().alpha_boundary(ig,lfsu_s,u_s,lfsv_s,downview_s);

        // jiggle in self
        for (int j=0; j<n_s; j++)
        {
          for (size_type k=0; k<y_s.size(); k++) up_s[k]=0.0;
          D delta = epsilon*(1.0+std::abs(u_s[lfsu_s.localIndex(j)]));
          u_s[lfsu_s.localIndex(j)] += delta;
          asImp().alpha_boundary(ig,lfsu_s,u_s,lfsv_s,upview_s);
          for (int i=0; i<m_s; i++)
            y_s.rawAccumulate(lfsv_s.localIndex(i),
                              ((up_s[lfsv_s.localIndex(i)]-down_s[lfsv_s.localIndex(i)])/delta)*x_s[lfsu_s.localIndex(j)]);
          u_s[lfsu_s.localIndex(j)] = x_s[lfsu_s.localIndex(j)];
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
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
        typedef LocalMatrix<typename R::value_type> Jacobian;
        typedef typename Jacobian::WeightedAccumulationView JacobianView;

        Jacobian mat(r.size(),x.size(), 0);
        JacobianView matview = mat.weightedAccumulationView(1.0);
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
        mat_ss.umv(x_s,r_s);
        mat_ns.umv(x_s,r_n);
        mat_sn.umv(x_n,r_s);
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
        typedef LocalMatrix<typename R::value_type> Jacobian;
        typedef typename Jacobian::WeightedAccumulationView JacobianView;

        Jacobian mat(x.size(),r.size(), 0);
        JacobianView view = mat.weightedAccumulationView(1.0);
        asImp().jacobian_boundary(ig, lfsu, x, lfsv, view);
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
