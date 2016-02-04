// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_NUMERICALJACOBIAN_HH
#define DUNE_PDELAB_LOCALOPERATOR_NUMERICALJACOBIAN_HH

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

        // Notice that in general lfsv.size() != mat.nrows()
        ResidualVector down(mat.nrows(),0.),up(mat.nrows());
        ResidualView downview = down.weightedAccumulationView(mat.weight());
        ResidualView upview = up.weightedAccumulationView(mat.weight());


        asImp().alpha_volume(eg,lfsu,u,lfsv,downview);
        for (int j=0; j<n; j++) // loop over columns
        {
          up = 0.0;
          D delta = epsilon*(1.0+std::abs(u(lfsu,j)));
          u(lfsu,j) += delta;
          asImp().alpha_volume(eg,lfsu,u,lfsv,upview);
          for (int i=0; i<m; i++)
            mat.rawAccumulate(lfsv,i,lfsu,j,(up(lfsv,i)-down(lfsv,i))/delta);
          u(lfsu,j) = x(lfsu,j);
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

        // Notice that in general lfsv.size() != mat.nrows()
        ResidualVector down(mat.nrows(),0.),up(mat.nrows());
        ResidualView downview = down.weightedAccumulationView(mat.weight());
        ResidualView upview = up.weightedAccumulationView(mat.weight());

        asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,downview);
        for (int j=0; j<n; j++) // loop over columns
        {
          up = 0.0;
          D delta = epsilon*(1.0+std::abs(u(lfsu,j)));
          u(lfsu,j) += delta;
          asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,upview);
          for (int i=0; i<m; i++)
            mat.rawAccumulate(lfsv,i,lfsu,j,(up(lfsv,i)-down(lfsv,i))/delta);
          u(lfsu,j) = x(lfsu,j);
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

        // Notice that in general lfsv.size() != mat.nrows()
        ResidualVector down_s(mat_ss.nrows()),up_s(mat_ss.nrows());
        ResidualView downview_s = down_s.weightedAccumulationView(1.0);
        ResidualView upview_s = up_s.weightedAccumulationView(1.0);

        ResidualVector down_n(mat_nn.nrows()),up_n(mat_nn.nrows());
        ResidualView downview_n = down_n.weightedAccumulationView(1.0);
        ResidualView upview_n = up_n.weightedAccumulationView(1.0);

        // base line
        asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,downview_s,
                               downview_n);

        // jiggle in self
        for (int j=0; j<n_s; j++)
        {
          up_s = 0.0;
          up_n = 0.0;
          D delta = epsilon*(1.0+std::abs(u_s(lfsu_s,j)));
          u_s(lfsu_s,j) += delta;
          asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,upview_s,
                                 upview_n);
          for (int i=0; i<m_s; i++)
            mat_ss.accumulate(lfsv_s,i,lfsu_s,j,(up_s(lfsv_s,i)-down_s(lfsv_s,i))/delta);
          for (int i=0; i<m_n; i++)
            mat_ns.accumulate(lfsv_n,i,lfsu_s,j,(up_n(lfsv_n,i)-down_n(lfsv_n,i))/delta);
          u_s(lfsu_s,j) = x_s(lfsu_s,j);
        }

        // jiggle in neighbor
        for (int j=0; j<n_n; j++)
        {
          up_s = 0.0;
          up_n = 0.0;
          D delta = epsilon*(1.0+std::abs(u_n(lfsu_n,j)));
          u_n(lfsu_n,j) += delta;
          asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,upview_s,
                                 upview_n);
          for (int i=0; i<m_s; i++)
            mat_sn.accumulate(lfsv_s,i,lfsu_n,j,(up_s(lfsv_s,i)-down_s(lfsv_s,i))/delta);
          for (int i=0; i<m_n; i++)
            mat_nn.accumulate(lfsv_n,i,lfsu_n,j,(up_n(lfsv_n,i)-down_n(lfsv_n,i))/delta);
          u_n(lfsu_n,j) = x_n(lfsu_n,j);
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

        // Notice that in general lfsv.size() != mat.nrows()
        ResidualVector down_s(mat_ss.nrows()),up_s(mat_ss.nrows());
        ResidualView downview_s = down_s.weightedAccumulationView(mat_ss.weight());
        ResidualView upview_s = up_s.weightedAccumulationView(mat_ss.weight());;


        // base line
        asImp().alpha_boundary(ig,lfsu_s,u_s,lfsv_s,downview_s);

        // jiggle in self
        for (int j=0; j<n_s; j++)
        {
          up_s = 0.0;
          D delta = epsilon*(1.0+std::abs(u_s(lfsu_s,j)));
          u_s(lfsu_s,j) += delta;
          asImp().alpha_boundary(ig,lfsu_s,u_s,lfsv_s,upview_s);
          for (int i=0; i<m_s; i++)
            mat_ss.rawAccumulate(lfsv_s,i,lfsu_s,j,(up_s(lfsv_s,i)-down_s(lfsv_s,i))/delta);
          u_s(lfsu_s,j) = x_s(lfsu_s,j);
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! \} group LocalOperatorDefaultImp
  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_NUMERICALJACOBIAN_HH
