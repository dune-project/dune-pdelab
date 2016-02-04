// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_NUMERICALJACOBIANAPPLY_HH
#define DUNE_PDELAB_LOCALOPERATOR_NUMERICALJACOBIANAPPLY_HH

#include <dune/pdelab/localoperator/numericalnonlinearjacobianapply.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup LocalOperatorDefaultImp
    //! \ingroup LocalOperator
    //! \{

    ////////////////////////////////////////////////////////////////////////
    //
    //  Numerical implementation of jacobian_apply_*() in terms of alpha_*()
    //

    //! Implements linear and nonlinear versions of jacobian_apply_volume() based on alpha_volume()
    /**
     * Derive from this class to add numerical jacobian application for
     * volume.  The derived class needs to implement alpha_volume().
     *
     * \note This mixin is designed for linear problems and provides both the
     *       linear and the nonlinear methods by inheriting from the nonlinear
     *       mixin.
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalJacobianApplyVolume
      : public NumericalNonlinearJacobianApplyVolume<Imp>
    {
    public:

      // We need to reimport the name from the base class; otherwise clang plays stupid and refuses to
      // find the overload in the base class
      using NumericalNonlinearJacobianApplyVolume<Imp>::jacobian_apply_volume;

      NumericalJacobianApplyVolume ()
        : NumericalNonlinearJacobianApplyVolume<Imp>(1e-7)
        , epsilon(1e-7)
      {}

      NumericalJacobianApplyVolume (double epsilon_)
        : NumericalNonlinearJacobianApplyVolume<Imp>(epsilon_)
        , epsilon(epsilon_)
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
        typedef typename Y::value_type R;
        typedef LocalVector<R,TestSpaceTag,typename Y::weight_type> ResidualVector;
        typedef typename ResidualVector::WeightedAccumulationView ResidualView;

        const int m=lfsv.size();
        const int n=lfsu.size();

        X u(x);

        // Notice that in general lfsv.size() != y.size()
        ResidualVector down(y.size()),up(y.size());
        ResidualView downview = down.weightedAccumulationView(y.weight());
        ResidualView upview = up.weightedAccumulationView(y.weight());

        asImp().alpha_volume(eg,lfsu,u,lfsv,downview);
        for (int j=0; j<n; j++) // loop over columns
        {
          up = 0.0;
          D delta = epsilon*(1.0+std::abs(u(lfsu,j)));
          u(lfsu,j) += delta;
          asImp().alpha_volume(eg,lfsu,u,lfsv,upview);
          for (int i=0; i<m; i++)
            y.rawAccumulate(lfsv,i,((up(lfsv,i)-down(lfsv,i))/delta)*x(lfsu,j));
          u(lfsu,j) = x(lfsu,j);
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! Implements linear and nonlinear versions jacobian_apply_volume_post_skeleton() based on
    //! alpha_volume_post_skeleton()
    /**
     * Derive from this class to add numerical jacobian application for volume
     * (post skeleton part).  The derived class needs to implement
     * alpha_volume().
     *
     * \note This mixin is designed for linear problems and provides both the
     *       linear and the nonlinear methods by inheriting from the nonlinear
     *       mixin.
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalJacobianApplyVolumePostSkeleton
      : public NumericalNonlinearJacobianApplyVolumePostSkeleton<Imp>
    {
    public:

      // We need to reimport the name from the base class; otherwise clang plays stupid and refuses to
      // find the overload in the base class
      using NumericalNonlinearJacobianApplyVolumePostSkeleton<Imp>::jacobian_apply_volume_post_skeleton;

      NumericalJacobianApplyVolumePostSkeleton ()
        : NumericalNonlinearJacobianApplyVolumePostSkeleton<Imp>(1e-7)
        , epsilon(1e-7)
      {}

      NumericalJacobianApplyVolumePostSkeleton (double epsilon_)
        : NumericalNonlinearJacobianApplyVolumePostSkeleton<Imp>(epsilon_)
        , epsilon(epsilon_)
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
        typedef typename Y::value_type R;
        typedef LocalVector<R,TestSpaceTag,typename Y::weight_type> ResidualVector;
        typedef typename ResidualVector::WeightedAccumulationView ResidualView;

        const int m=lfsv.size();
        const int n=lfsu.size();

        X u(x);

        // Notice that in general lfsv.size() != y.size()
        ResidualVector down(y.size()),up(y.size());
        ResidualView downview = down.weightedAccumulationView(y.weight());
        ResidualView upview = up.weightedAccumulationView(y.weight());

        asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,downview);
        for (int j=0; j<n; j++) // loop over columns
        {
          up = 0.0;
          D delta = epsilon*(1.0+std::abs(u(lfsu,j)));
          u(lfsu,j) += delta;
          asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,upview);
          for (int i=0; i<m; i++)
            y.rawAccumulate(lfsv,i,((up(lfsv,i)-down(lfsv,i))/delta)*x(lfsu,j));
          u(lfsu,j) = x(lfsu,j);
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };

    //! Implements linear and nonlinear versions of jacobian_apply_skeleton() based on alpha_skeleton()
    /**
     * Derive from this class to add numerical jacobian application for
     * skeleton.  The derived class needs to implement alpha_skeleton().
     *
     * \note This mixin is designed for linear problems and provides both the
     *       linear and the nonlinear methods by inheriting from the nonlinear
     *       mixin.
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalJacobianApplySkeleton
      : public NumericalNonlinearJacobianApplySkeleton<Imp>
    {
    public:

      // We need to reimport the name from the base class; otherwise clang plays stupid and refuses to
      // find the overload in the base class
      using NumericalNonlinearJacobianApplySkeleton<Imp>::jacobian_apply_skeleton;

      NumericalJacobianApplySkeleton ()
        : NumericalNonlinearJacobianApplySkeleton<Imp>(1e-7)
        , epsilon(1e-7)
      {}

      NumericalJacobianApplySkeleton (double epsilon_)
        : NumericalNonlinearJacobianApplySkeleton<Imp>(epsilon_)
        , epsilon(epsilon_)
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
        typedef LocalVector<R,TestSpaceTag,typename Y::weight_type> ResidualVector;
        typedef typename ResidualVector::WeightedAccumulationView ResidualView;

        const int m_s=lfsv_s.size();
        const int m_n=lfsv_n.size();
        const int n_s=lfsu_s.size();
        const int n_n=lfsu_n.size();

        X u_s(x_s);
        X u_n(x_n);

        // Notice that in general lfsv_s.size() != y_s.size()
        ResidualVector down_s(y_s.size()),up_s(y_s.size());
        ResidualView downview_s = down_s.weightedAccumulationView(1.0);
        ResidualView upview_s = up_s.weightedAccumulationView(1.0);

        ResidualVector down_n(y_n.size()),up_n(y_n.size());
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
            y_s.accumulate(lfsv_s,i,((up_s(lfsv_s,i)-down_s(lfsv_s,i))/delta)*x_s(lfsu_s,j));
          for (int i=0; i<m_n; i++)
            y_n.accumulate(lfsv_n,i,((up_n(lfsv_n,i)-down_n(lfsv_n,i))/delta)*x_s(lfsu_s,j));
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
            y_s.accumulate(lfsv_s,i,((up_s(lfsv_s,i)-down_s(lfsv_s,i))/delta)*x_n(lfsu_n,j));
          for (int i=0; i<m_n; i++)
            y_n.accumulate(lfsv_n,i,((up_n(lfsv_n,i)-down_n(lfsv_n,i))/delta)*x_n(lfsu_n,j));
          u_n(lfsu_n,j) = x_n(lfsu_n,j);
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! Implements linear and nonlinear versions of jacobian_apply_boundary() based on alpha_boundary()
    /**
     * Derive from this class to add numerical jacobian application for
     * boundary.  The derived class needs to implement alpha_boundary().
     *
     * \note This mixin is designed for linear problems and provides both the
     *       linear and the nonlinear methods by inheriting from the nonlinear
     *       mixin.
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalJacobianApplyBoundary
      : public NumericalNonlinearJacobianApplyBoundary<Imp>
    {
    public:

      // We need to reimport the name from the base class; otherwise clang plays stupid and refuses to
      // find the overload in the base class
      using NumericalNonlinearJacobianApplyBoundary<Imp>::jacobian_apply_boundary;

      NumericalJacobianApplyBoundary ()
        : NumericalNonlinearJacobianApplyBoundary<Imp>(1e-7)
        , epsilon(1e-7)
      {}

      NumericalJacobianApplyBoundary (double epsilon_)
        : NumericalNonlinearJacobianApplyBoundary<Imp>(epsilon_)
        , epsilon(epsilon_)
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
        typedef LocalVector<R,TestSpaceTag,typename Y::weight_type> ResidualVector;
        typedef typename ResidualVector::WeightedAccumulationView ResidualView;

        const int m_s=lfsv_s.size();
        const int n_s=lfsu_s.size();

        X u_s(x_s);

        // Notice that in general lfsv_s.size() != y_s.size()
        ResidualVector down_s(y_s.size()),up_s(y_s.size());
        ResidualView downview_s = down_s.weightedAccumulationView(1.0);
        ResidualView upview_s = up_s.weightedAccumulationView(1.0);

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
            y_s.rawAccumulate(lfsv_s,i,((up_s(lfsv_s,i)-down_s(lfsv_s,i))/delta)*x_s(lfsu_s,j));
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

#endif // DUNE_PDELAB_LOCALOPERATOR_NUMERICALJACOBIANAPPLY_HH
