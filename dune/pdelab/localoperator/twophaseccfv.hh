// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_TWOPHASEOP_HH
#define DUNE_PDELAB_TWOPHASEOP_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/localfunctions/raviartthomas/raviartthomas0q.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>

#include"../common/function.hh"
#include"pattern.hh"
#include"flags.hh"
#include"idefault.hh"

namespace Dune {
  namespace PDELab {

    //! traits class for two phase parameter class
    template<typename GV, typename RF>
    struct TwoPhaseParameterTraits
    {
      //! \brief the grid view
      typedef GV GridViewType;

      //! \brief Enum for domain dimension
      enum {
        //! \brief dimension of the domain
        dimDomain = GV::dimension
      };

      //! \brief Export type for domain field
      typedef typename GV::Grid::ctype DomainFieldType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType,dimDomain> DomainType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType,dimDomain-1> IntersectionDomainType;

      //! \brief Export type for range field
      typedef RF RangeFieldType;

      //! \brief range type
      typedef Dune::FieldVector<RF,GV::dimensionworld> RangeType;

      //! \brief permeability tensor type
      typedef RangeFieldType PermTensorType;

      //! grid types
      typedef typename GV::Traits::template Codim<0>::Entity ElementType;
      typedef typename GV::Intersection IntersectionType;
    };

    template<typename GV, typename RF>
    struct TwoPhaseFullTensorParameterTraits : TwoPhaseParameterTraits<GV, RF>
    {
      typedef TwoPhaseParameterTraits<GV, RF> Base;
      typedef typename Base::RangeFieldType RangeFieldType;

      //! \brief permeability tensor type
      typedef Dune::FieldMatrix<RangeFieldType,Base::dimDomain,Base::dimDomain> PermTensorType;
    };

    //! base class for parameter class
    template<class T, class Imp>
    class TwoPhaseParameterInterface
    {
    public:
      typedef T Traits;

      //! porosity
      typename Traits::RangeFieldType
      phi (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().phi(e,x);
      }

      //! capillary pressure function
      typename Traits::RangeFieldType
      pc (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
          typename Traits::RangeFieldType s_l) const
      {
        return asImp().pc(e,x,s_l);
      }

      //! inverse capillary pressure function
      typename Traits::RangeFieldType
      s_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
           typename Traits::RangeFieldType pc) const
      {
        return asImp().s_l(e,x,pc);
      }

      //! liquid phase relative permeability
      typename Traits::RangeFieldType
      kr_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType s_l) const
      {
        return asImp().kr_l(e,x,s_l);
      }

      //! gas phase relative permeability
      typename Traits::RangeFieldType
      kr_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType s_g) const
      {
        return asImp().kr_g(e,x,s_g);
      }

      //! liquid phase viscosity
      typename Traits::RangeFieldType
      mu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType p_l) const
      {
        return asImp().mu_l(e,x,p_l);
      }

      //! gas phase viscosity
      typename Traits::RangeFieldType
      mu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType p_g) const
      {
        return asImp().mu_l(e,x,p_g);
      }

      //! absolute permeability (scalar!)
      typename Traits::PermTensorType
      k_abs (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().k_abs(e,x);
      }

      //! gravity vector
      const typename Traits::RangeType& gravity () const
      {
        return asImp().gravity();
      }

      //! liquid phase molar density
      template<typename E>
      typename Traits::RangeFieldType
      nu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType p_l) const
      {
        return asImp().nu_l(e,x,p_l);
      }

      //! gas phase molar density
      typename Traits::RangeFieldType
      nu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType p_g) const
      {
        return asImp().nu_g(e,x,p_g);
      }

      //! liquid phase mass density
      typename Traits::RangeFieldType
      rho_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
             typename Traits::RangeFieldType p_l) const
      {
        return asImp().rho_l(e,x,p_l);
      }

      //! gas phase mass density
      typename Traits::RangeFieldType
      rho_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
             typename Traits::RangeFieldType p_g) const
      {
        return asImp().rho_g(e,x,p_g);
      }

      //! liquid phase boundary condition type
      int
      bc_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
      {
        return asImp().bc_l(is,x,time);
      }

      //! gas phase boundary condition type
      int
      bc_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
      {
        return asImp().bc_g(is,x,time);
      }

      //! liquid phase Dirichlet boundary condition
      typename Traits::RangeFieldType
      g_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
      {
        return asImp().g_l(is,x,time);
      }

      //! gas phase Dirichlet boundary condition
      typename Traits::RangeFieldType
      g_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
      {
        return asImp().g_g(is,x,time);
      }

      //! liquid phase Neumann boundary condition
      typename Traits::RangeFieldType
      j_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
      {
        return asImp().j_l(is,x,time);
      }

      //! gas phase Neumann boundary condition
      typename Traits::RangeFieldType
      j_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
      {
        return asImp().j_g(is,x,time);
      }

      //! liquid phase source term
      typename Traits::RangeFieldType
      q_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
           typename Traits::RangeFieldType time) const
      {
        return asImp().q_l(e,x,time);
      }

      //! gas phase source term
      typename Traits::RangeFieldType
      q_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
           typename Traits::RangeFieldType time) const
      {
        return asImp().q_g(e,x,time);
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };


    // a local operator for solving the two-phase flow in pressure-pressure formulation
    // with two-point flux approximation
    // TP : parameter class, see above
    // V  : Vector holding last time step
    template<typename TP>
    class TwoPhaseTwoPointFluxOperator
      : public NumericalJacobianSkeleton<TwoPhaseTwoPointFluxOperator<TP> >,
        public NumericalJacobianApplySkeleton<TwoPhaseTwoPointFluxOperator<TP> >,

        public NumericalJacobianBoundary<TwoPhaseTwoPointFluxOperator<TP> >,
        public NumericalJacobianApplyBoundary<TwoPhaseTwoPointFluxOperator<TP> >,

        public FullSkeletonPattern,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,

        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = 0 };
      enum { gas = 1 };

      typedef typename TP::Traits::RangeFieldType Real;
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume   = true };
      enum { doLambdaBoundary = true };

      //! constructor: pass parameter object
      TwoPhaseTwoPointFluxOperator (const TP& tp_, Real scale_l_=1.0, Real scale_g_=1.0)
        : tp(tp_), scale_l(scale_l_), scale_g(scale_g_)
      {}

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<liquid>::Type PLSpace;
        typedef typename LFSV::template Child<gas>::Type PGSpace;

        // domain and range field type
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);
        RF cell_volume = eg.geometry().volume();

        // contribution from source term
        r.accumulate(lfsv, liquid, -scale_l * tp.q_l(eg.entity(),cell_center_local,time) * cell_volume);
        r.accumulate(lfsv, gas, -scale_g * tp.q_g(eg.entity(),cell_center_local,time) * cell_volume);
      }

      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // select the two components
        typedef typename LFSV::template Child<0>::Type PLSpace;

        // domain and range field type
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        // cell geometries
        const Dune::FieldVector<DF,dim>&
          inside_cell_center_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>&
          outside_cell_center_local = Dune::ReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);
        Dune::FieldVector<DF,IG::dimension>
          inside_cell_center_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension>
          outside_cell_center_global = ig.outside()->geometry().center();

        // distance of cell centers
        Dune::FieldVector<DF,dim> d(outside_cell_center_global);
        d -= inside_cell_center_global;
        RF distance = d.two_norm();

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();

        // absolute permeability
        RF k_abs_inside = tp.k_abs(*(ig.inside()),inside_cell_center_local);
        RF k_abs_outside = tp.k_abs(*(ig.outside()),outside_cell_center_local);

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);

        // liquid phase calculation
        RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF rho_l_outside = tp.rho_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid));
        RF w_l = (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid))/distance + aavg(rho_l_inside,rho_l_outside)*gn; // determines direction
        RF pc_upwind, s_l_upwind, s_g_upwind;
        RF nu_l = aavg(tp.nu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid)),
                       tp.nu_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid)));
        if (w_l>=0) // upwind capillary pressure on face
          {
            pc_upwind = x_s(lfsu_s,gas)-x_s(lfsu_s,liquid);
            s_l_upwind = tp.s_l(*(ig.inside()),inside_cell_center_local,pc_upwind);
          }
        else
          {
            pc_upwind = x_n(lfsu_n,gas)-x_n(lfsu_n,liquid);
            s_l_upwind = tp.s_l(*(ig.outside()),outside_cell_center_local,pc_upwind);
          }
        s_g_upwind = 1-s_l_upwind;
        RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l_upwind)/
          tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF lambda_l_outside = tp.kr_l(*(ig.outside()),outside_cell_center_local,s_l_upwind)/
          tp.mu_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid));
        RF sigma_l = havg(lambda_l_inside*k_abs_inside,lambda_l_outside*k_abs_outside);

        r_s.accumulate(lfsv_s,liquid, scale_l * nu_l * sigma_l * w_l * face_volume);
        r_n.accumulate(lfsv_n,liquid, -scale_l * nu_l * sigma_l * w_l * face_volume);

        // gas phase calculation
        RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas));
        RF rho_g_outside = tp.rho_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,gas));
        RF w_g = (x_s(lfsu_s,gas)-x_n(lfsu_n,gas))/distance + aavg(rho_g_inside,rho_g_outside)*gn; // determines direction
        RF nu_g = aavg(tp.nu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas)),
                       tp.nu_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,gas)));
        if (w_l*w_g<0) // new evaluation necessary only if signs differ
          {
            if (w_g>=0) // upwind capillary pressure on face
              {
                pc_upwind = x_s(lfsu_s,gas)-x_s(lfsu_s,liquid);
                s_l_upwind = tp.s_l(*(ig.inside()),inside_cell_center_local,pc_upwind);
              }
            else
              {
                pc_upwind = x_n(lfsu_n,gas)-x_n(lfsu_n,liquid);
                s_l_upwind = tp.s_l(*(ig.outside()),outside_cell_center_local,pc_upwind);
              }
            s_g_upwind = 1-s_l_upwind;
          }
        RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g_upwind)/
          tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas));
        RF lambda_g_outside = tp.kr_g(*(ig.outside()),outside_cell_center_local,s_g_upwind)/
          tp.mu_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,gas));
        RF sigma_g = havg(lambda_g_inside*k_abs_inside,lambda_g_outside*k_abs_outside);

        r_s.accumulate(lfsv_s, gas, scale_g * nu_g * sigma_g * w_g * face_volume);
        r_n.accumulate(lfsv_n, gas, -scale_g * nu_g * sigma_g * w_g * face_volume);
      }

      // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // select the two components
        typedef typename LFSV::template Child<0>::Type PLSpace;

        // domain and range field type
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        // face geometry
        const Dune::FieldVector<DF,dim-1>&
          face_local = Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();

        // evaluate boundary condition type
        int bc_l = tp.bc_l(ig.intersection(),face_local,time);
        int bc_g = tp.bc_g(ig.intersection(),face_local,time);
        if (bc_l!=1 && bc_g!=1) return; // no Dirichlet boundary conditions

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          inside_cell_center_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        Dune::FieldVector<DF,dim>
          inside_cell_center_global = ig.inside()->geometry().global(inside_cell_center_local);

        // distance of cell center to boundary
        Dune::FieldVector<DF,dim> d = ig.geometry().global(face_local);
        d -= inside_cell_center_global;
        RF distance = d.two_norm();

        // absolute permeability
        RF k_abs_inside = tp.k_abs(*(ig.inside()),inside_cell_center_local);

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);

        // liquid phase Dirichlet boundary
        if (bc_l==1)
          {
            RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF g_l = tp.g_l(ig.intersection(),face_local,time);
            RF w_l = (x_s(lfsu_s,liquid)-g_l)/distance + rho_l_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas)-x_s(lfsu_s,liquid));
            RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l)/
              tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF sigma_l = lambda_l_inside*k_abs_inside;
            RF nu_l = tp.nu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            r_s.accumulate(lfsv_s, liquid, scale_l * nu_l * sigma_l * w_l * face_volume);
          }

        // gas phase Dirichlet boundary
        if (bc_g==1)
          {
            RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas));
            RF g_g = tp.g_g(ig.intersection(),face_local,time);
            RF w_g = (x_s(lfsu_s,gas)-g_g)/distance + rho_g_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas)-x_s(lfsu_s,liquid));
            RF s_g = 1-s_l;
            RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g)/
              tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas));
            RF sigma_g = lambda_g_inside*k_abs_inside;
            RF nu_g = tp.nu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas));
            r_s.accumulate(lfsv_s, gas, scale_l * nu_g * sigma_g * w_g * face_volume);
          }
      }

      // boundary integral independent of ansatz functions
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r_s) const
      {
        // select the two components
        typedef typename LFSV::template Child<0>::Type PLSpace;

        // domain and range field type
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        // face geometry
        const Dune::FieldVector<DF,dim-1>&
          face_local = Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().integrationElement(face_local)*
          Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).volume();

        // evaluate boundary condition type
        int bc_l = tp.bc_l(ig.intersection(),face_local,time);
        int bc_g = tp.bc_g(ig.intersection(),face_local,time);
        if (bc_l!=0 && bc_g!=0) return; // no Neumann boundary conditions

        // liquid phase Neumann boundary
        if (bc_l==0)
          {
            RF j_l = tp.j_l(ig.intersection(),face_local,time);
            r_s.accumulate(lfsv, liquid, scale_l * j_l * face_volume);
          }

        // gas phase Neumann boundary
        if (bc_g==0)
          {
            RF j_g = tp.j_g(ig.intersection(),face_local,time);
            r_s.accumulate(lfsv, gas, scale_g * j_g * face_volume);
          }
      }

      //! set time for subsequent evaluation
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
      }

    private:
      const TP& tp;  // two phase parameter class
      typename TP::Traits::RangeFieldType time;
      Real scale_l, scale_g;

      template<typename T>
      T aavg (T a, T b) const
      {
        return 0.5*(a+b);
      }

      template<typename T>
      T havg (T a, T b) const
      {
        T eps = 1E-30;
        return 2.0/(1.0/(a+eps) + 1.0/(b+eps));
      }
    };


    /** a local operator for the storage operator
     *
     * \f{align*}{
     \int_\Omega c(x) uv dx
     * \f}
     */
    template<class TP>
    class TwoPhaseOnePointTemporalOperator
      : public NumericalJacobianVolume<TwoPhaseOnePointTemporalOperator<TP> >,
        public NumericalJacobianApplyVolume<TwoPhaseOnePointTemporalOperator<TP> >,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = 0 };
      enum { gas = 1 };

      typedef typename TP::Traits::RangeFieldType Real;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      TwoPhaseOnePointTemporalOperator (TP& tp_, Real scale_l_=1.0, Real scale_g_=1.0)
        : tp(tp_), scale_l(scale_l_), scale_g(scale_g_)
      {
      }

      //! set time for subsequent evaluation
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<0>::Type PLSpace;

        // domain and range field type
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);
        RF cell_volume = eg.geometry().volume();

        RF phi = tp.phi(eg.entity(),cell_center_local);
        RF s_l = tp.s_l(eg.entity(),cell_center_local,x(lfsu,gas)-x(lfsu,liquid));

        r.accumulate(lfsv, liquid, scale_l * phi * s_l * tp.nu_l(eg.entity(),cell_center_local,x(lfsu,liquid)) * cell_volume);
        r.accumulate(lfsv, gas, scale_g * phi * (1-s_l) * tp.nu_g(eg.entity(),cell_center_local,x(lfsu,gas)) * cell_volume);
      }

    private:
      TP& tp;
      typename TP::Traits::RangeFieldType time;
      Real scale_l, scale_g;
    };


    /** \brief Provide velocity field for liquid phase

        Uses RT0 interpolation on a cell.

        - T  : provides TwoPhaseParameterInterface
        - PL : P0 function for liquid phase pressure
        - PG : P0 function for gas phase pressure
    */
    template<typename  TP, typename PL, typename PG>
    class V_l
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                               typename PL::Traits::RangeFieldType,
                                                                               PL::Traits::GridViewType::dimension,
                                                                               Dune::FieldVector<typename PL::Traits::RangeFieldType,PL::Traits::GridViewType::dimension> >,
                                              V_l<TP,PL,PG> >
    {
      // extract useful types
      typedef typename PL::Traits::GridViewType GV;
      typedef typename GV::Grid::ctype DF;
      typedef typename PL::Traits::RangeFieldType RF;
      typedef typename PL::Traits::RangeType RangeType;
      enum { dim = PL::Traits::GridViewType::dimension };
      typedef typename GV::Traits::template Codim<0>::Entity Element;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      typedef typename IntersectionIterator::Intersection Intersection;

      const TP& tp;
      const PL& pl;
      const PG& pg;
      Dune::RT0QLocalFiniteElement<DF,RF,dim> rt0fe;
      typename TP::Traits::RangeFieldType time;


      typedef typename Dune::RT0QLocalFiniteElement<DF,RF,dim>::Traits::LocalBasisType::Traits::RangeType RT0RangeType;

    public:
      typedef Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> > Traits;

      typedef Dune::PDELab::GridFunctionBase<Traits,V_l<TP,PL,PG> > BaseT;

      V_l (const TP& tp_, const PL& pl_, const PG& pg_) : tp(tp_), pl(pl_), pg(pg_), time(0) {}

      // set time where operator is to be evaluated (i.e. end of the time intervall)
      void set_time (typename TP::Traits::RangeFieldType time_)
      {
        time = time_;
      }

      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        // cell geometry
        const Dune::FieldVector<DF,dim>&
          inside_cell_center_local = Dune::ReferenceElements<DF,dim>::
          general(e.type()).position(0,0);
        Dune::FieldVector<DF,dim>
          inside_cell_center_global = e.geometry().global(inside_cell_center_local);

        // absolute permeability
        RF k_abs_inside = tp.k_abs(e,inside_cell_center_local);

        // pressure evaluation
        typename PL::Traits::RangeType pl_inside, pg_inside;
        pl.evaluate(e,inside_cell_center_local,pl_inside);
        pg.evaluate(e,inside_cell_center_local,pg_inside);

        // density
        RF nu_l_inside = tp.nu_l(e,inside_cell_center_local,pl_inside);

        // for coefficient computation
        RF vn[2*dim];    // normal velocities
        RF coeff[2*dim]; // RT0 coefficient
        typename Traits::ElementType::Geometry::JacobianInverseTransposed
          B = e.geometry().jacobianInverseTransposed(x); // the transformation. Assume it is linear
        RF determinant = B.determinant();

        // loop over cell neighbors
        IntersectionIterator endit = pl.getGridView().iend(e);
        for (IntersectionIterator iit = pl.getGridView().ibegin(e); iit!=endit; ++iit)
          {
            // set to zero for processor boundary
            vn[iit->indexInInside()] = 0.0;

            // face geometry
            const Dune::FieldVector<DF,dim-1>&
              face_local = Dune::ReferenceElements<DF,dim-1>::general(iit->geometry().type()).position(0,0);


            // interior face
            if (iit->neighbor())
              {
                const Dune::FieldVector<DF,dim>&
                  outside_cell_center_local = Dune::ReferenceElements<DF,dim>::
                  general(iit->outside()->type()).position(0,0);
                Dune::FieldVector<DF,dim>
                  outside_cell_center_global = iit->outside()->geometry().global(outside_cell_center_local);

                // distance of cell centers
                Dune::FieldVector<DF,dim> d(outside_cell_center_global);
                d -= inside_cell_center_global;
                RF distance = d.two_norm();

                // absolute permeability
                RF k_abs_outside = tp.k_abs(*(iit->outside()),outside_cell_center_local);

                // gravity times normal
                RF gn = tp.gravity()*iit->unitOuterNormal(face_local);

                // pressure evaluation
                typename PL::Traits::RangeType pl_outside, pg_outside;
                pl.evaluate(*(iit->outside()),outside_cell_center_local,pl_outside);
                pg.evaluate(*(iit->outside()),outside_cell_center_local,pg_outside);

                // density
                RF nu_l_outside = tp.nu_l(*(iit->outside()),outside_cell_center_local,pg_outside);

                // liquid phase calculation
                RF rho_l_inside = tp.rho_l(e,inside_cell_center_local,pl_inside);
                RF rho_l_outside = tp.rho_l(*(iit->outside()),outside_cell_center_local,pl_outside);
                RF w_l = (pl_inside-pl_outside)/distance + aavg(rho_l_inside,rho_l_outside)*gn; // determines direction
                RF pc_upwind, s_l_upwind, s_g_upwind;
                if (w_l>=0) // upwind capillary pressure on face
                  {
                    pc_upwind = pg_inside-pl_inside;
                    s_l_upwind = tp.s_l(e,inside_cell_center_local,pc_upwind);
                  }
                else
                  {
                    pc_upwind = pg_outside-pl_outside;
                    s_l_upwind = tp.s_l(*(iit->outside()),outside_cell_center_local,pc_upwind);
                  }
                s_g_upwind = 1-s_l_upwind;
                RF lambda_l_inside = tp.kr_l(e,inside_cell_center_local,s_l_upwind)/
                  tp.mu_l(e,inside_cell_center_local,pl_inside);
                RF lambda_l_outside = tp.kr_l(*(iit->outside()),outside_cell_center_local,s_l_upwind)/
                  tp.mu_l(*(iit->outside()),outside_cell_center_local,pl_outside);
                RF sigma_l = havg(lambda_l_inside*k_abs_inside,lambda_l_outside*k_abs_outside);
                RF nu_l = aavg(nu_l_inside,nu_l_outside);

                // set coefficient
                vn[iit->indexInInside()] = nu_l * sigma_l * w_l;
              }

            // boundary face
            if (iit->boundary())
              {
                // distance of cell center to boundary
                Dune::FieldVector<DF,dim> d = iit->geometry().global(face_local);
                d -= inside_cell_center_global;
                RF distance = d.two_norm();

                // evaluate boundary condition type
                int bc_l = tp.bc_l(*iit,face_local,time);

                // liquid phase Dirichlet boundary
                if (bc_l==1)
                  {
                    RF rho_l_inside = tp.rho_l(e,inside_cell_center_local,pl_inside);
                    RF g_l = tp.g_l(*iit,face_local,time);
                    RF gn = tp.gravity()*iit->unitOuterNormal(face_local);
                    RF w_l = (pl_inside-g_l)/distance + rho_l_inside*gn;
                    RF s_l = tp.s_l(e,inside_cell_center_local,pg_inside-pl_inside);
                    RF lambda_l_inside = tp.kr_l(e,inside_cell_center_local,s_l)/
                      tp.mu_l(e,inside_cell_center_local,pl_inside);
                    RF sigma_l = lambda_l_inside*k_abs_inside;
                    vn[iit->indexInInside()] = nu_l_inside * sigma_l * w_l;
                  }

                // liquid phase Neumann boundary
                if (bc_l==0)
                  {
                    RF j_l = tp.j_l(*iit,face_local,time);
                    vn[iit->indexInInside()] = j_l;
                  }
              }

            // compute coefficient
            Dune::FieldVector<DF,dim> vstar=iit->unitOuterNormal(face_local); // normal on tranformef element
            vstar *= vn[iit->indexInInside()];
            Dune::FieldVector<RF,dim> normalhat(0); // normal on reference element
            if (iit->indexInInside()%2==0)
              normalhat[iit->indexInInside()/2] = -1.0;
            else
              normalhat[iit->indexInInside()/2] =  1.0;
            Dune::FieldVector<DF,dim> vstarhat(0);
            B.umtv(vstar,vstarhat); // Piola backward transformation
            vstarhat *= determinant;
            coeff[iit->indexInInside()] = vstarhat*normalhat;
          }

        // compute velocity on reference element
        std::vector<RT0RangeType> rt0vectors(rt0fe.localBasis().size());
        rt0fe.localBasis().evaluateFunction(x,rt0vectors);
        typename Traits::RangeType yhat(0);
        for (unsigned int i=0; i<rt0fe.localBasis().size(); i++)
          yhat.axpy(coeff[i],rt0vectors[i]);

        // apply Piola transformation
        B.invert();
        y = 0;
        B.umtv(yhat,y);
        y /= determinant;

        //         std::cout << "vn= " ;
        //         for (int i=0; i<2*dim; i++) std::cout << vn[i] << " ";
        //         std::cout << std::endl;
        //         std::cout << "V_l: x=" << x << " y=" << y << std::endl;
      }

      inline const typename Traits::GridViewType& getGridView ()
      {
        return pl.getGridView();
      }

    private:

      template<typename T>
      T aavg (T a, T b) const
      {
        return 0.5*(a+b);
      }

      template<typename T>
      T havg (T a, T b) const
      {
        T eps = 1E-30;
        return 2.0/(1.0/(a+eps) + 1.0/(b+eps));
      }
    };

    /** \brief Provide velocity field for gas phase

        Uses RT0 interpolation on a cell.

        - T  : provides TwoPhaseParameterInterface
        - PL : P0 function for liquid phase pressure
        - PG : P0 function for gas phase pressure
    */
    template<typename  TP, typename PL, typename PG>
    class V_g
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                               typename PL::Traits::RangeFieldType,
                                                                               PL::Traits::GridViewType::dimension,
                                                                               Dune::FieldVector<typename PL::Traits::RangeFieldType,PL::Traits::GridViewType::dimension> >,
                                              V_g<TP,PL,PG> >
    {
      // extract useful types
      typedef typename PL::Traits::GridViewType GV;
      typedef typename GV::Grid::ctype DF;
      typedef typename PL::Traits::RangeFieldType RF;
      typedef typename PL::Traits::RangeType RangeType;
      enum { dim = PL::Traits::GridViewType::dimension };
      typedef typename GV::Traits::template Codim<0>::Entity Element;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      typedef typename IntersectionIterator::Intersection Intersection;

      const TP& tp;
      const PL& pl;
      const PG& pg;
      Dune::RT0QLocalFiniteElement<DF,RF,dim> rt0fe;
      typename TP::Traits::RangeFieldType time;


      typedef typename Dune::RT0QLocalFiniteElement<DF,RF,dim>::Traits::LocalBasisType::Traits::RangeType RT0RangeType;

    public:
      typedef Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> > Traits;

      typedef Dune::PDELab::GridFunctionBase<Traits,V_g<TP,PL,PG> > BaseT;

      V_g (const TP& tp_, const PL& pl_, const PG& pg_) : tp(tp_), pl(pl_), pg(pg_), time(0) {}

      // set time where operator is to be evaluated (i.e. end of the time intervall)
      void set_time (typename TP::Traits::RangeFieldType time_)
      {
        time = time_;
      }

      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        // cell geometry
        const Dune::FieldVector<DF,dim>&
          inside_cell_center_local = Dune::ReferenceElements<DF,dim>::
          general(e.type()).position(0,0);
        Dune::FieldVector<DF,dim>
          inside_cell_center_global = e.geometry().global(inside_cell_center_local);

        // absolute permeability
        RF k_abs_inside = tp.k_abs(e,inside_cell_center_local);

        // pressure evaluation
        typename PL::Traits::RangeType pl_inside, pg_inside;
        pl.evaluate(e,inside_cell_center_local,pl_inside);
        pg.evaluate(e,inside_cell_center_local,pg_inside);

        // density evaluation
        RF nu_g_inside = tp.nu_g(e,inside_cell_center_local,pg_inside);

        // for coefficient computation
        RF vn[2*dim];    // normal velocities
        RF coeff[2*dim]; // RT0 coefficient
        typename Traits::ElementType::Geometry::JacobianInverseTransposed
          B = e.geometry().jacobianInverseTransposed(x); // the transformation. Assume it is linear
        RF determinant = B.determinant();

        // loop over cell neighbors
        IntersectionIterator endit = pl.getGridView().iend(e);
        for (IntersectionIterator iit = pl.getGridView().ibegin(e); iit!=endit; ++iit)
          {
            // set to zero for processor boundary
            vn[iit->indexInInside()] = 0.0;

            // face geometry
            const Dune::FieldVector<DF,dim-1>&
              face_local = Dune::ReferenceElements<DF,dim-1>::general(iit->geometry().type()).position(0,0);

            // interior face
            if (iit->neighbor())
              {
                const Dune::FieldVector<DF,dim>&
                  outside_cell_center_local = Dune::ReferenceElements<DF,dim>::
                  general(iit->outside()->type()).position(0,0);
                Dune::FieldVector<DF,dim>
                  outside_cell_center_global = iit->outside()->geometry().global(outside_cell_center_local);

                // distance of cell centers
                Dune::FieldVector<DF,dim> d(outside_cell_center_global);
                d -= inside_cell_center_global;
                RF distance = d.two_norm();

                // absolute permeability
                RF k_abs_outside = tp.k_abs(*(iit->outside()),outside_cell_center_local);

                // gravity times normal
                RF gn = tp.gravity()*iit->unitOuterNormal(face_local);

                // pressure evaluation
                typename PL::Traits::RangeType pl_outside, pg_outside;
                pl.evaluate(*(iit->outside()),outside_cell_center_local,pl_outside);
                pg.evaluate(*(iit->outside()),outside_cell_center_local,pg_outside);

                // needed for both phases
                RF pc_upwind, s_l_upwind, s_g_upwind;

                // gas phase calculation
                RF rho_g_inside = tp.rho_g(e,inside_cell_center_local,pg_inside);
                RF rho_g_outside = tp.rho_g(e,outside_cell_center_local,pg_outside);
                RF w_g = (pg_inside-pg_outside)/distance + aavg(rho_g_inside,rho_g_outside)*gn; // determines direction
                if (w_g>=0) // upwind capillary pressure on face
                  {
                    pc_upwind = pg_inside-pl_inside;
                    s_l_upwind = tp.s_l(e,inside_cell_center_local,pc_upwind);
                  }
                else
                  {
                    pc_upwind = pg_outside-pl_outside;
                    s_l_upwind = tp.s_l(*(iit->outside()),outside_cell_center_local,pc_upwind);
                  }
                s_g_upwind = 1-s_l_upwind;
                RF lambda_g_inside = tp.kr_g(e,inside_cell_center_local,s_g_upwind)/
                  tp.mu_g(e,inside_cell_center_local,pg_inside);
                RF lambda_g_outside = tp.kr_g(*(iit->outside()),outside_cell_center_local,s_g_upwind)/
                  tp.mu_g(*(iit->outside()),outside_cell_center_local,pg_outside);
                RF sigma_g = havg(lambda_g_inside*k_abs_inside,lambda_g_outside*k_abs_outside);

                RF nu_g_outside = tp.nu_g(*(iit->outside()),outside_cell_center_local,pg_outside);

                // set coefficient
                vn[iit->indexInInside()] = aavg(nu_g_inside,nu_g_outside) * sigma_g * w_g;
              }

            // boundary face
            if (iit->boundary())
              {
                // distance of cell center to boundary
                Dune::FieldVector<DF,dim> d = iit->geometry().global(face_local);
                d -= inside_cell_center_global;
                RF distance = d.two_norm();

                // evaluate boundary condition type
                int bc_g = tp.bc_g(*iit,face_local,time);

                // gas phase Dirichlet boundary
                if (bc_g==1)
                  {
                    RF rho_g_inside = tp.rho_g(e,inside_cell_center_local,pg_inside);
                    RF g_g = tp.g_g(*iit,face_local,time);
                    RF gn = tp.gravity()*iit->unitOuterNormal(face_local);
                    RF w_g = (pg_inside-g_g)/distance + rho_g_inside*gn;
                    RF s_l = tp.s_l(e,inside_cell_center_local,pg_inside-pl_inside);
                    RF s_g = 1-s_l;
                    RF lambda_g_inside = tp.kr_g(e,inside_cell_center_local,s_g)/
                      tp.mu_g(e,inside_cell_center_local,pg_inside);
                    RF sigma_g = lambda_g_inside*k_abs_inside;

                    vn[iit->indexInInside()] = nu_g_inside * sigma_g * w_g;
                  }

                // gas phase Neumann boundary
                if (bc_g==0)
                  {
                    RF j_g = tp.j_g(*iit,face_local,time);
                    vn[iit->indexInInside()] = j_g; /* /nu_g_inside*/;
                  }
              }

            // compute coefficient
            Dune::FieldVector<DF,dim> vstar=iit->unitOuterNormal(face_local); // normal on tranformef element
            vstar *= vn[iit->indexInInside()];
            Dune::FieldVector<RF,dim> normalhat(0); // normal on reference element
            if (iit->indexInInside()%2==0)
              normalhat[iit->indexInInside()/2] = -1.0;
            else
              normalhat[iit->indexInInside()/2] =  1.0;
            Dune::FieldVector<DF,dim> vstarhat(0);
            B.umtv(vstar,vstarhat); // Piola backward transformation
            vstarhat *= determinant;
            coeff[iit->indexInInside()] = vstarhat*normalhat;
          }

        // compute velocity on reference element
        std::vector<RT0RangeType> rt0vectors(rt0fe.localBasis().size());
        rt0fe.localBasis().evaluateFunction(x,rt0vectors);
        typename Traits::RangeType yhat(0);
        for (unsigned int i=0; i<rt0fe.localBasis().size(); i++)
          yhat.axpy(coeff[i],rt0vectors[i]);

        // apply Piola transformation
        B.invert();
        y = 0;
        B.umtv(yhat,y);
        y /= determinant;
      }

      inline const typename Traits::GridViewType& getGridView ()
      {
        return pl.getGridView();
      }

    private:

      template<typename T>
      T aavg (T a, T b) const
      {
        return 0.5*(a+b);
      }

      template<typename T>
      T havg (T a, T b) const
      {
        T eps = 1E-30;
        return 2.0/(1.0/(a+eps) + 1.0/(b+eps));
      }
    };


  }
}

#endif
