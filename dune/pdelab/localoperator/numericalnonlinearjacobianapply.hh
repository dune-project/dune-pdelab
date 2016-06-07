// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_NUMERICALNONLINEARJACOBIANAPPLY_HH
#define DUNE_PDELAB_LOCALOPERATOR_NUMERICALNONLINEARJACOBIANAPPLY_HH

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
    //  Numerical implementation of nonlinear jacobian_apply_*() in terms of alpha_*()
    //

    //! Implements nonlinear version of jacobian_apply_volume() based on alpha_volume()
    /**
     * Derive from this class to add numerical jacobian application for
     * volume.  The derived class needs to implement alpha_volume().
     *
     * \note This mixin is designed for nonlinear problems.
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalNonlinearJacobianApplyVolume
    {
    public:
      NumericalNonlinearJacobianApplyVolume ()
        : epsilon(1e-7)
      {}

      NumericalNonlinearJacobianApplyVolume (double epsilon_)
        : epsilon(epsilon_)
      {}

      //! apply local jacobian of the volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_volume(
        const EG& eg,
        const LFSU& lfsu, const X& x, const X& z,
        const LFSV& lfsv, Y& y) const
      {
        using R = typename Y::value_type;
        using ResidualVector = LocalVector<R,TestSpaceTag,typename Y::weight_type>;

        auto m = lfsv.size();
        auto n = lfsu.size();

        X u(x);

        // Notice that in general lfsv.size() != y.size()
        ResidualVector down(y.size()),up(y.size());
        auto downview = down.weightedAccumulationView(y.weight());
        auto upview = up.weightedAccumulationView(y.weight());

        asImp().alpha_volume(eg,lfsu,u,lfsv,downview);
        for (decltype(n) j = 0; j < n; ++j) // loop over columns
        {
          using namespace std;
          up = 0.0;
          R delta = epsilon*(1.0 + abs(u(lfsu,j)));
          u(lfsu,j) += delta;
          asImp().alpha_volume(eg,lfsu,u,lfsv,upview);
          for (decltype(m) i = 0; i < m; ++i)
            y.rawAccumulate(lfsv,i,((up(lfsv,i)-down(lfsv,i))/delta)*z(lfsu,j));
          u(lfsu,j) = x(lfsu,j);
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! Implements nonlinear version of jacobian_apply_volume_post_skeleton() based on
    //! alpha_volume_post_skeleton()
    /**
     * Derive from this class to add numerical jacobian application for volume
     * (post skeleton part).  The derived class needs to implement
     * alpha_volume().
     *
     * \note This mixin is designed for nonlinear problems.
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalNonlinearJacobianApplyVolumePostSkeleton
    {
    public:
      NumericalNonlinearJacobianApplyVolumePostSkeleton ()
        : epsilon(1e-7)
      {}

      NumericalNonlinearJacobianApplyVolumePostSkeleton (double epsilon_)
        : epsilon(epsilon_)
      {}

      //! apply local jacobian of the volume term (post skeleton part)
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_volume_post_skeleton(
        const EG& eg,
        const LFSU& lfsu, const X& x, const X& z,
        const LFSV& lfsv, Y& y) const
      {
        using R = typename Y::value_type;
        using ResidualVector = LocalVector<R,TestSpaceTag,typename Y::weight_type>;

        auto m = lfsv.size();
        auto n = lfsu.size();

        X u(x);

        // Notice that in general lfsv.size() != y.size()
        ResidualVector down(y.size()),up(y.size());
        auto downview = down.weightedAccumulationView(y.weight());
        auto upview = up.weightedAccumulationView(y.weight());

        asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,downview);
        for (decltype(n) j = 0; j < n; ++j) // loop over columns
        {
          using namespace std;
          up = 0.0;
          R delta = epsilon*(1.0 + abs(u(lfsu,j)));
          u(lfsu,j) += delta;
          asImp().alpha_volume_post_skeleton(eg,lfsu,u,lfsv,upview);
          for (decltype(m) i = 0; i < m; ++i)
            y.rawAccumulate(lfsv,i,((up(lfsv,i)-down(lfsv,i))/delta)*z(lfsu,j));
          u(lfsu,j) = x(lfsu,j);
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };

    //! Implements nonlinear version of jacobian_apply_skeleton() based on alpha_skeleton()
    /**
     * Derive from this class to add numerical jacobian application for
     * skeleton.  The derived class needs to implement alpha_skeleton().
     *
     * \note This mixin is designed for nonlinear problems.
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalNonlinearJacobianApplySkeleton
    {
    public:
      NumericalNonlinearJacobianApplySkeleton ()
        : epsilon(1e-7)
      {}

      NumericalNonlinearJacobianApplySkeleton (double epsilon_)
        : epsilon(epsilon_)
      {}

      //! apply local jacobian of the skeleton term
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_skeleton(
        const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
        const LFSU& lfsu_n, const X& x_n, const X& z_n, const LFSV& lfsv_n,
        Y& y_s, Y& y_n) const
      {
        using R = typename X::value_type;
        using ResidualVector = LocalVector<R,TestSpaceTag,typename Y::weight_type>;

        auto m_s=lfsv_s.size();
        auto m_n=lfsv_n.size();
        auto n_s=lfsu_s.size();
        auto n_n=lfsu_n.size();

        X u_s(x_s);
        X u_n(x_n);

        // Notice that in general lfsv_s.size() != y_s.size()
        ResidualVector down_s(y_s.size()),up_s(y_s.size());
        auto downview_s = down_s.weightedAccumulationView(1.0);
        auto upview_s = up_s.weightedAccumulationView(1.0);

        ResidualVector down_n(y_n.size()),up_n(y_n.size());
        auto downview_n = down_n.weightedAccumulationView(1.0);
        auto upview_n = up_n.weightedAccumulationView(1.0);

        // base line
        asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,downview_s,downview_n);

        // jiggle in self
        for (decltype(n_s) j = 0; j < n_s; ++j)
        {
          using namespace std;
          up_s = 0.0;
          up_n = 0.0;
          R delta = epsilon*(1.0 + abs(u_s(lfsu_s,j)));
          u_s(lfsu_s,j) += delta;
          asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,upview_s,upview_n);
          for (decltype(m_s) i = 0; i < m_s; ++i)
            y_s.accumulate(lfsv_s,i,((up_s(lfsv_s,i)-down_s(lfsv_s,i))/delta)*z_s(lfsu_s,j));
          for (decltype(m_n) i = 0; i < m_n; i++)
            y_n.accumulate(lfsv_n,i,((up_n(lfsv_n,i)-down_n(lfsv_n,i))/delta)*x_s(lfsu_s,j));
          u_s(lfsu_s,j) = x_s(lfsu_s,j);
        }

        // jiggle in neighbor
        for (decltype(n_n) j = 0; j < n_n; ++j)
        {
          using namespace std;
          up_s = 0.0;
          up_n = 0.0;
          R delta = epsilon*(1.0+abs(u_n(lfsu_n,j)));
          u_n(lfsu_n,j) += delta;
          asImp().alpha_skeleton(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,upview_s,upview_n);
          for (decltype(m_s) i = 0; i < m_s; ++i)
            y_s.accumulate(lfsv_s,i,((up_s(lfsv_s,i)-down_s(lfsv_s,i))/delta)*z_n(lfsu_n,j));
          for (decltype(m_n) i = 0; i < m_n; ++i)
            y_n.accumulate(lfsv_n,i,((up_n(lfsv_n,i)-down_n(lfsv_n,i))/delta)*x_n(lfsu_n,j));
          u_n(lfsu_n,j) = x_n(lfsu_n,j);
        }
      }

    private:
      const double epsilon; // problem: this depends on data type R!
      Imp& asImp () { return static_cast<Imp &> (*this); }
      const Imp& asImp () const { return static_cast<const Imp &>(*this); }
    };

    //! Implements nonlinear version of jacobian_apply_boundary() based on alpha_boundary()
    /**
     * Derive from this class to add numerical jacobian application for
     * boundary.  The derived class needs to implement alpha_boundary().
     *
     * \note This mixin is designed for nonlinear problems.
     *
     * \tparam Imp Type of the derived class (CRTP-trick).
     */
    template<typename Imp>
    class NumericalNonlinearJacobianApplyBoundary
    {
    public:
      NumericalNonlinearJacobianApplyBoundary ()
        : epsilon(1e-7)
      {}

      NumericalNonlinearJacobianApplyBoundary (double epsilon_)
        : epsilon(epsilon_)
      {}

      //! apply local jacobian of the boundaryterm
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_boundary(
        const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const X& z_s,
        const LFSV& lfsv_s, Y& y_s) const
      {
        using R = typename Y::value_type;
        using ResidualVector = LocalVector<R,TestSpaceTag,typename Y::weight_type>;

        auto m_s=lfsv_s.size();
        auto n_s=lfsu_s.size();

        X u_s(x_s);

        // Notice that in general lfsv_s.size() != y_s.size()
        ResidualVector down_s(y_s.size()),up_s(y_s.size());
        auto downview_s = down_s.weightedAccumulationView(1.0);
        auto upview_s = up_s.weightedAccumulationView(1.0);

        // base line
        asImp().alpha_boundary(ig,lfsu_s,u_s,lfsv_s,downview_s);

        // jiggle in self
        for (decltype(n_s) j = 0; j < n_s; ++j)
        {
          up_s = 0.0;
          R delta = epsilon*(1.0+std::abs(u_s(lfsu_s,j)));
          u_s(lfsu_s,j) += delta;
          asImp().alpha_boundary(ig,lfsu_s,u_s,lfsv_s,upview_s);
          for (decltype(m_s) i = 0; i < m_s; ++i)
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

#endif // DUNE_PDELAB_LOCALOPERATOR_NUMERICALNONLINEARJACOBIANAPPLY_HH
