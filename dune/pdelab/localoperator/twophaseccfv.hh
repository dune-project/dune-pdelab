// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_TWOPHASECCFV_HH
#define DUNE_PDELAB_LOCALOPERATOR_TWOPHASECCFV_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/localfunctions/raviartthomas/raviartthomascube.hh>

#include<dune/pdelab/common/quadraturerules.hh>
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
      using GridViewType = GV;

      //! \brief Enum for domain dimension
      enum {
        //! \brief dimension of the domain
        dimDomain = GV::dimension
      };

      //! \brief Export type for domain field
      using DomainFieldType = typename GV::Grid::ctype;

      //! \brief domain type
      using DomainType = Dune::FieldVector<DomainFieldType,dimDomain>;

      //! \brief domain type
      using IntersectionDomainType = Dune::FieldVector<DomainFieldType,dimDomain-1>;

      //! \brief Export type for range field
      using RangeFieldType = RF;

      //! \brief range type
      using RangeType = Dune::FieldVector<RF,GV::dimensionworld>;

      //! \brief permeability tensor type
      using PermTensorType = RangeFieldType;

      //! grid types
      using ElementType = typename GV::Traits::template Codim<0>::Entity;
      using IntersectionType = typename GV::Intersection;
    };

    template<typename GV, typename RF>
    struct TwoPhaseFullTensorParameterTraits : TwoPhaseParameterTraits<GV, RF>
    {
      using Base = TwoPhaseParameterTraits<GV, RF>;
      using RangeFieldType = typename Base::RangeFieldType;

      //! \brief permeability tensor type
      using PermTensorType = Dune::FieldMatrix<RangeFieldType,Base::dimDomain,Base::dimDomain>;
    };

    //! base class for parameter class
    template<class T, class Imp>
    class TwoPhaseParameterInterface
    {
    public:
      using Traits = T;

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

      using Real = typename TP::Traits::RangeFieldType;
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
        // Reference to cell
        const auto& cell = eg.entity();

        // get geometry
        auto geo = eg.geometry();

        // cell geometry
        auto ref_el = referenceElement(geo);
        auto cell_center_local = ref_el.position(0,0);
        auto cell_volume = geo.volume();

        // contribution from source term
        r.accumulate(lfsv, liquid, -scale_l * tp.q_l(cell,cell_center_local,time) * cell_volume);
        r.accumulate(lfsv, gas, -scale_g * tp.q_g(cell,cell_center_local,time) * cell_volume);
      }

      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // Define types
        using namespace Indices;
        using PLSpace = TypeTree::Child<LFSV,_0>;
        using RF = typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;

        // References to inside and outside cells
        const auto& cell_inside = ig.inside();
        const auto& cell_outside = ig.outside();

        // get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();
        auto geo_outside = cell_outside.geometry();

        // cell geometries
        auto ref_el_inside = referenceElement(geo_inside);
        auto ref_el_outside = referenceElement(geo_outside);
        auto inside_cell_center_local = ref_el_inside.position(0,0);
        auto outside_cell_center_local = ref_el_outside.position(0,0);
        auto inside_cell_center_global = geo_inside.center();
        auto outside_cell_center_global = geo_outside.center();

        // distance of cell centers
        auto d = outside_cell_center_global;
        d -= inside_cell_center_global;
        auto distance = d.two_norm();

        // face geometry
        auto ref_el = referenceElement(geo);
        auto face_local = ref_el.position(0,0);
        auto face_volume = geo.volume();

        // absolute permeability
        auto k_abs_inside = tp.k_abs(cell_inside,inside_cell_center_local);
        auto k_abs_outside = tp.k_abs(cell_outside,outside_cell_center_local);

        // gravity times normal
        auto gn = tp.gravity()*ig.unitOuterNormal(face_local);

        // liquid phase calculation
        auto rho_l_inside = tp.rho_l(cell_inside,inside_cell_center_local,x_s(lfsu_s,liquid));
        auto rho_l_outside = tp.rho_l(cell_outside,outside_cell_center_local,x_n(lfsu_n,liquid));
        auto w_l = (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid))/distance + aavg(rho_l_inside,rho_l_outside)*gn; // determines direction
        RF pc_upwind, s_l_upwind, s_g_upwind;
        auto nu_l = aavg(tp.nu_l(cell_inside,inside_cell_center_local,x_s(lfsu_s,liquid)),
                         tp.nu_l(cell_outside,outside_cell_center_local,x_n(lfsu_n,liquid)));
        if (w_l>=0) // upwind capillary pressure on face
          {
            pc_upwind = x_s(lfsu_s,gas)-x_s(lfsu_s,liquid);
            s_l_upwind = tp.s_l(cell_inside,inside_cell_center_local,pc_upwind);
          }
        else
          {
            pc_upwind = x_n(lfsu_n,gas)-x_n(lfsu_n,liquid);
            s_l_upwind = tp.s_l(cell_outside,outside_cell_center_local,pc_upwind);
          }
        s_g_upwind = 1-s_l_upwind;
        auto lambda_l_inside = tp.kr_l(cell_inside,inside_cell_center_local,s_l_upwind)/
          tp.mu_l(cell_inside,inside_cell_center_local,x_s(lfsu_s,liquid));
        auto lambda_l_outside = tp.kr_l(cell_outside,outside_cell_center_local,s_l_upwind)/
          tp.mu_l(cell_outside,outside_cell_center_local,x_n(lfsu_n,liquid));
        auto sigma_l = havg(lambda_l_inside*k_abs_inside,lambda_l_outside*k_abs_outside);

        r_s.accumulate(lfsv_s,liquid, scale_l * nu_l * sigma_l * w_l * face_volume);
        r_n.accumulate(lfsv_n,liquid, -scale_l * nu_l * sigma_l * w_l * face_volume);

        // gas phase calculation
        auto rho_g_inside = tp.rho_g(cell_inside,inside_cell_center_local,x_s(lfsu_s,gas));
        auto rho_g_outside = tp.rho_g(cell_outside,outside_cell_center_local,x_n(lfsu_n,gas));
        auto w_g = (x_s(lfsu_s,gas)-x_n(lfsu_n,gas))/distance + aavg(rho_g_inside,rho_g_outside)*gn; // determines direction
        auto nu_g = aavg(tp.nu_g(cell_inside,inside_cell_center_local,x_s(lfsu_s,gas)),
                         tp.nu_g(cell_outside,outside_cell_center_local,x_n(lfsu_n,gas)));
        if (w_l*w_g<=0) // new evaluation necessary only if signs differ
          {
            if (w_g>=0) // upwind capillary pressure on face
              {
                pc_upwind = x_s(lfsu_s,gas)-x_s(lfsu_s,liquid);
                s_l_upwind = tp.s_l(cell_inside,inside_cell_center_local,pc_upwind);
              }
            else
              {
                pc_upwind = x_n(lfsu_n,gas)-x_n(lfsu_n,liquid);
                s_l_upwind = tp.s_l(cell_outside,outside_cell_center_local,pc_upwind);
              }
            s_g_upwind = 1-s_l_upwind;
          }
        auto lambda_g_inside = tp.kr_g(cell_inside,inside_cell_center_local,s_g_upwind)/
          tp.mu_g(cell_inside,inside_cell_center_local,x_s(lfsu_s,gas));
        auto lambda_g_outside = tp.kr_g(cell_outside,outside_cell_center_local,s_g_upwind)/
          tp.mu_g(cell_outside,outside_cell_center_local,x_n(lfsu_n,gas));
        auto sigma_g = havg(lambda_g_inside*k_abs_inside,lambda_g_outside*k_abs_outside);

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
        // References to inside cell
        const auto& cell_inside = ig.inside();

        // get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // face geometry
        auto ref_el = referenceElement(geo);
        auto face_local = ref_el.position(0,0);
        auto face_volume = geo.volume();

        // evaluate boundary condition type
        auto bc_l = tp.bc_l(ig.intersection(),face_local,time);
        auto bc_g = tp.bc_g(ig.intersection(),face_local,time);
        if (bc_l!=1 && bc_g!=1) return; // no Dirichlet boundary conditions

        // cell geometry
        auto ref_el_inside = referenceElement(geo_inside);
        auto inside_cell_center_local = ref_el_inside.position(0,0);
        auto inside_cell_center_global = geo_inside.center();

        // distance of cell center to boundary
        auto d = geo.global(face_local);
        d -= inside_cell_center_global;
        auto distance = d.two_norm();

        // absolute permeability
        auto k_abs_inside = tp.k_abs(cell_inside,inside_cell_center_local);

        // gravity times normal
        auto gn = tp.gravity()*ig.unitOuterNormal(face_local);

        // liquid phase Dirichlet boundary
        if (bc_l==1)
          {
            auto rho_l_inside = tp.rho_l(cell_inside,inside_cell_center_local,x_s(lfsu_s,liquid));
            auto g_l = tp.g_l(ig.intersection(),face_local,time);
            auto w_l = (x_s(lfsu_s,liquid)-g_l)/distance + rho_l_inside*gn;
            auto s_l = tp.s_l(cell_inside,inside_cell_center_local,x_s(lfsu_s,gas)-x_s(lfsu_s,liquid));
            auto lambda_l_inside = tp.kr_l(cell_inside,inside_cell_center_local,s_l)/
              tp.mu_l(cell_inside,inside_cell_center_local,x_s(lfsu_s,liquid));
            auto sigma_l = lambda_l_inside*k_abs_inside;
            auto nu_l = tp.nu_l(cell_inside,inside_cell_center_local,x_s(lfsu_s,liquid));
            r_s.accumulate(lfsv_s, liquid, scale_l * nu_l * sigma_l * w_l * face_volume);
          }

        // gas phase Dirichlet boundary
        if (bc_g==1)
          {
            auto rho_g_inside = tp.rho_g(cell_inside,inside_cell_center_local,x_s(lfsu_s,gas));
            auto g_g = tp.g_g(ig.intersection(),face_local,time);
            auto w_g = (x_s(lfsu_s,gas)-g_g)/distance + rho_g_inside*gn;
            auto s_l = tp.s_l(cell_inside,inside_cell_center_local,x_s(lfsu_s,gas)-x_s(lfsu_s,liquid));
            auto s_g = 1-s_l;
            auto lambda_g_inside = tp.kr_g(cell_inside,inside_cell_center_local,s_g)/
              tp.mu_g(cell_inside,inside_cell_center_local,x_s(lfsu_s,gas));
            auto sigma_g = lambda_g_inside*k_abs_inside;
            auto nu_g = tp.nu_g(cell_inside,inside_cell_center_local,x_s(lfsu_s,gas));
            r_s.accumulate(lfsv_s, gas, scale_l * nu_g * sigma_g * w_g * face_volume);
          }
      }

      // boundary integral independent of ansatz functions
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r_s) const
      {
        // get geometries
        auto geo = ig.geometry();

        // face geometry
        auto ref_el = referenceElement(geo);
        auto face_local = ref_el.position(0,0);
        auto face_volume = geo.volume();

        // evaluate boundary condition type
        auto bc_l = tp.bc_l(ig.intersection(),face_local,time);
        auto bc_g = tp.bc_g(ig.intersection(),face_local,time);
        if (bc_l!=0 && bc_g!=0) return; // no Neumann boundary conditions

        // liquid phase Neumann boundary
        if (bc_l==0)
          {
            auto j_l = tp.j_l(ig.intersection(),face_local,time);
            r_s.accumulate(lfsv, liquid, scale_l * j_l * face_volume);
          }

        // gas phase Neumann boundary
        if (bc_g==0)
          {
            auto j_g = tp.j_g(ig.intersection(),face_local,time);
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

      using Real = typename TP::Traits::RangeFieldType;

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
        // Reference to cell
        const auto& cell = eg.entity();

        // get geometry
        auto geo = eg.geometry();

        // cell geometry
        auto ref_el = referenceElement(geo);
        auto cell_center_local = ref_el.position(0,0);
        auto cell_volume = geo.volume();

        auto phi = tp.phi(cell,cell_center_local);
        auto s_l = tp.s_l(cell,cell_center_local,x(lfsu,gas)-x(lfsu,liquid));

        r.accumulate(lfsv, liquid, scale_l * phi * s_l * tp.nu_l(cell,cell_center_local,x(lfsu,liquid)) * cell_volume);
        r.accumulate(lfsv, gas, scale_g * phi * (1-s_l) * tp.nu_g(cell,cell_center_local,x(lfsu,gas)) * cell_volume);
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
      using GV = typename PL::Traits::GridViewType;
      using DF = typename GV::Grid::ctype;
      using RF = typename PL::Traits::RangeFieldType;
      using RangeType = typename PL::Traits::RangeType;
      enum { dim = PL::Traits::GridViewType::dimension };
      using Element = typename GV::Traits::template Codim<0>::Entity;
      using IntersectionIterator = typename GV::IntersectionIterator;
      using Intersection = typename IntersectionIterator::Intersection;

      const TP& tp;
      const PL& pl;
      const PG& pg;
      Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0> rt0fe;
      typename TP::Traits::RangeFieldType time;


      using RT0RangeType = typename Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0>::Traits::LocalBasisType::Traits::RangeType;

    public:
      using Traits = Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> >;

      using BaseT = Dune::PDELab::GridFunctionBase<Traits,V_l<TP,PL,PG> >;

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
       // get geometries
        auto geo = e.geometry();

        // face geometry
        auto ref_el = referenceElement(geo);
        auto face_local = ref_el.position(0,0);
        auto face_volume = geo.volume();

        // cell geometry
        auto inside_cell_center_local = ref_el.position(0,0);
        auto inside_cell_center_global = geo.global(inside_cell_center_local);

        // absolute permeability
        auto k_abs_inside = tp.k_abs(e,inside_cell_center_local);

        // pressure evaluation
        typename PL::Traits::RangeType pl_inside, pg_inside;
        pl.evaluate(e,inside_cell_center_local,pl_inside);
        pg.evaluate(e,inside_cell_center_local,pg_inside);

        // density
        auto nu_l_inside = tp.nu_l(e,inside_cell_center_local,pl_inside);

        // for coefficient computation
        RF vn[2*dim];    // normal velocities
        RF coeff[2*dim]; // RT0 coefficient
        auto B = geo.jacobianInverseTransposed(x); // the transformation. Assume it is linear
        auto determinant = B.determinant();

        // loop over cell neighbors
        for (const auto& intersection : intersections(pl.getGridView(),e))
          {
            // set to zero for processor boundary
            vn[intersection.indexInInside()] = 0.0;

            // get geometry
            auto geo_intersection = intersection.geometry();

            // face geometry
            const auto& face_local = referenceElement(geo_intersection).position(0,0);

            // interior face
            if (intersection.neighbor())
              {
                // get geometry
                auto geo_outside = intersection.outside().geometry();
                auto ref_el_outside = referenceElement(geo_outside);
                auto outside_cell_center_local = ref_el_outside.position(0,0);
                auto outside_cell_center_global = geo_outside.global(outside_cell_center_local);

                // distance of cell centers
                auto d = outside_cell_center_global;
                d -= inside_cell_center_global;
                auto distance = d.two_norm();

                // absolute permeability
                auto k_abs_outside = tp.k_abs(*(intersection.outside()),outside_cell_center_local);

                // gravity times normal
                auto gn = tp.gravity()*intersection.unitOuterNormal(face_local);

                // pressure evaluation
                typename PL::Traits::RangeType pl_outside, pg_outside;
                pl.evaluate(*(intersection.outside()),outside_cell_center_local,pl_outside);
                pg.evaluate(*(intersection.outside()),outside_cell_center_local,pg_outside);

                // density
                auto nu_l_outside = tp.nu_l(*(intersection.outside()),outside_cell_center_local,pg_outside);

                // liquid phase calculation
                auto rho_l_inside = tp.rho_l(e,inside_cell_center_local,pl_inside);
                auto rho_l_outside = tp.rho_l(*(intersection.outside()),outside_cell_center_local,pl_outside);
                auto w_l = (pl_inside-pl_outside)/distance + aavg(rho_l_inside,rho_l_outside)*gn; // determines direction
                RF pc_upwind, s_l_upwind, s_g_upwind;
                if (w_l>=0) // upwind capillary pressure on face
                  {
                    pc_upwind = pg_inside-pl_inside;
                    s_l_upwind = tp.s_l(e,inside_cell_center_local,pc_upwind);
                  }
                else
                  {
                    pc_upwind = pg_outside-pl_outside;
                    s_l_upwind = tp.s_l(*(intersection.outside()),outside_cell_center_local,pc_upwind);
                  }
                s_g_upwind = 1-s_l_upwind;
                auto lambda_l_inside = tp.kr_l(e,inside_cell_center_local,s_l_upwind)/
                  tp.mu_l(e,inside_cell_center_local,pl_inside);
                auto lambda_l_outside = tp.kr_l(*(intersection.outside()),outside_cell_center_local,s_l_upwind)/
                  tp.mu_l(*(intersection.outside()),outside_cell_center_local,pl_outside);
                auto sigma_l = havg(lambda_l_inside*k_abs_inside,lambda_l_outside*k_abs_outside);
                auto nu_l = aavg(nu_l_inside,nu_l_outside);

                // set coefficient
                vn[intersection.indexInInside()] = nu_l * sigma_l * w_l;
              }

            // boundary face
            if (intersection.boundary())
              {
                // distance of cell center to boundary
                auto d = geo_intersection.global(face_local);
                d -= inside_cell_center_global;
                auto distance = d.two_norm();

                // evaluate boundary condition type
                auto bc_l = tp.bc_l(intersection,face_local,time);

                // liquid phase Dirichlet boundary
                if (bc_l==1)
                  {
                    auto rho_l_inside = tp.rho_l(e,inside_cell_center_local,pl_inside);
                    auto g_l = tp.g_l(intersection,face_local,time);
                    auto gn = tp.gravity()*intersection.unitOuterNormal(face_local);
                    auto w_l = (pl_inside-g_l)/distance + rho_l_inside*gn;
                    auto s_l = tp.s_l(e,inside_cell_center_local,pg_inside-pl_inside);
                    auto lambda_l_inside = tp.kr_l(e,inside_cell_center_local,s_l)/
                      tp.mu_l(e,inside_cell_center_local,pl_inside);
                    auto sigma_l = lambda_l_inside*k_abs_inside;
                    vn[intersection.indexInInside()] = nu_l_inside * sigma_l * w_l;
                  }

                // liquid phase Neumann boundary
                if (bc_l==0)
                  {
                    auto j_l = tp.j_l(intersection,face_local,time);
                    vn[intersection.indexInInside()] = j_l;
                  }
              }

            // compute coefficient
            auto vstar = intersection.unitOuterNormal(face_local); // normal on tranformef element
            vstar *= vn[intersection.indexInInside()];
            Dune::FieldVector<RF,dim> normalhat(0); // normal on reference element
            if (intersection.indexInInside()%2==0)
              normalhat[intersection.indexInInside()/2] = -1.0;
            else
              normalhat[intersection.indexInInside()/2] =  1.0;
            Dune::FieldVector<DF,dim> vstarhat(0);
            B.umtv(vstar,vstarhat); // Piola backward transformation
            vstarhat *= determinant;
            coeff[intersection.indexInInside()] = vstarhat*normalhat;
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
      using GV = typename PL::Traits::GridViewType;
      using DF = typename GV::Grid::ctype;
      using RF = typename PL::Traits::RangeFieldType;
      using RangeType = typename PL::Traits::RangeType;
      enum { dim = PL::Traits::GridViewType::dimension };
      using Element = typename GV::Traits::template Codim<0>::Entity;
      using IntersectionIterator = typename GV::IntersectionIterator;
      using Intersection = typename IntersectionIterator::Intersection;

      const TP& tp;
      const PL& pl;
      const PG& pg;
      Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0> rt0fe;
      typename TP::Traits::RangeFieldType time;


      using RT0RangeType = typename Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0>::Traits::LocalBasisType::Traits::RangeType;

    public:
      using Traits = Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> >;

      using BaseT = Dune::PDELab::GridFunctionBase<Traits,V_g<TP,PL,PG> >;

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
      // get geometries
        auto geo = e.geometry();

        // face geometry
        auto ref_el = referenceElement(geo);
        auto face_local = ref_el.position(0,0);
        auto face_volume = geo.volume();

        // cell geometry
        auto inside_cell_center_local = ref_el.position(0,0);
        auto inside_cell_center_global = geo.global(inside_cell_center_local);

        // absolute permeability
        auto k_abs_inside = tp.k_abs(e,inside_cell_center_local);

        // pressure evaluation
        typename PL::Traits::RangeType pl_inside, pg_inside;
        pl.evaluate(e,inside_cell_center_local,pl_inside);
        pg.evaluate(e,inside_cell_center_local,pg_inside);

        // density evaluation
        auto nu_g_inside = tp.nu_g(e,inside_cell_center_local,pg_inside);

        // for coefficient computation
        RF vn[2*dim];    // normal velocities
        RF coeff[2*dim]; // RT0 coefficient
        auto B = geo.jacobianInverseTransposed(x); // the transformation. Assume it is linear
        auto determinant = B.determinant();

        // loop over cell neighbors
        for (const auto& intersection : intersections(pl.getGridView(),e))
          {
            // set to zero for processor boundary
            vn[intersection.indexInInside()] = 0.0;

            // get geometry
            auto geo_intersection = intersection.geometry();

            // face geometry
            const auto& face_local = referenceElement(geo_intersection).position(0,0);

            // interior face
            if (intersection.neighbor())
              {
                // get geometry
                auto geo_outside = intersection.outside().geometry();
                auto ref_el_outside = referenceElement(geo_outside);
                auto outside_cell_center_local = ref_el_outside.position(0,0);
                auto outside_cell_center_global = geo_outside.global(outside_cell_center_local);

                // distance of cell centers
                auto d = outside_cell_center_global;
                d -= inside_cell_center_global;
                auto distance = d.two_norm();

                // absolute permeability
                auto k_abs_outside = tp.k_abs(*(intersection.outside()),outside_cell_center_local);

                // gravity times normal
                auto gn = tp.gravity()*intersection.unitOuterNormal(face_local);

                // pressure evaluation
                typename PL::Traits::RangeType pl_outside, pg_outside;
                pl.evaluate(*(intersection.outside()),outside_cell_center_local,pl_outside);
                pg.evaluate(*(intersection.outside()),outside_cell_center_local,pg_outside);

                // needed for both phases
                RF pc_upwind, s_l_upwind, s_g_upwind;

                // gas phase calculation
                auto rho_g_inside = tp.rho_g(e,inside_cell_center_local,pg_inside);
                auto rho_g_outside = tp.rho_g(e,outside_cell_center_local,pg_outside);
                auto w_g = (pg_inside-pg_outside)/distance + aavg(rho_g_inside,rho_g_outside)*gn; // determines direction
                if (w_g>=0) // upwind capillary pressure on face
                  {
                    pc_upwind = pg_inside-pl_inside;
                    s_l_upwind = tp.s_l(e,inside_cell_center_local,pc_upwind);
                  }
                else
                  {
                    pc_upwind = pg_outside-pl_outside;
                    s_l_upwind = tp.s_l(*(intersection.outside()),outside_cell_center_local,pc_upwind);
                  }
                s_g_upwind = 1-s_l_upwind;
                auto lambda_g_inside = tp.kr_g(e,inside_cell_center_local,s_g_upwind)/
                  tp.mu_g(e,inside_cell_center_local,pg_inside);
                auto lambda_g_outside = tp.kr_g(*(intersection.outside()),outside_cell_center_local,s_g_upwind)/
                  tp.mu_g(*(intersection.outside()),outside_cell_center_local,pg_outside);
                auto sigma_g = havg(lambda_g_inside*k_abs_inside,lambda_g_outside*k_abs_outside);

                auto nu_g_outside = tp.nu_g(*(intersection.outside()),outside_cell_center_local,pg_outside);

                // set coefficient
                vn[intersection.indexInInside()] = aavg(nu_g_inside,nu_g_outside) * sigma_g * w_g;
              }

            // boundary face
            if (intersection.boundary())
              {
                // distance of cell center to boundary
                auto d = geo_intersection.global(face_local);
                d -= inside_cell_center_global;
                auto distance = d.two_norm();

                // evaluate boundary condition type
                auto bc_g = tp.bc_g(intersection,face_local,time);

                // gas phase Dirichlet boundary
                if (bc_g==1)
                  {
                    auto rho_g_inside = tp.rho_g(e,inside_cell_center_local,pg_inside);
                    auto g_g = tp.g_g(intersection,face_local,time);
                    auto gn = tp.gravity()*intersection.unitOuterNormal(face_local);
                    auto w_g = (pg_inside-g_g)/distance + rho_g_inside*gn;
                    auto s_l = tp.s_l(e,inside_cell_center_local,pg_inside-pl_inside);
                    auto s_g = 1-s_l;
                    auto lambda_g_inside = tp.kr_g(e,inside_cell_center_local,s_g)/
                      tp.mu_g(e,inside_cell_center_local,pg_inside);
                    auto sigma_g = lambda_g_inside*k_abs_inside;

                    vn[intersection.indexInInside()] = nu_g_inside * sigma_g * w_g;
                  }

                // gas phase Neumann boundary
                if (bc_g==0)
                  {
                    auto j_g = tp.j_g(intersection,face_local,time);
                    vn[intersection.indexInInside()] = j_g; /* /nu_g_inside*/;
                  }
              }

            // compute coefficient
            auto vstar = intersection.unitOuterNormal(face_local); // normal on tranformef element
            vstar *= vn[intersection.indexInInside()];
            Dune::FieldVector<RF,dim> normalhat(0); // normal on reference element
            if (intersection.indexInInside()%2==0)
              normalhat[intersection.indexInInside()/2] = -1.0;
            else
              normalhat[intersection.indexInInside()/2] =  1.0;
            Dune::FieldVector<DF,dim> vstarhat(0);
            B.umtv(vstar,vstarhat); // Piola backward transformation
            vstarhat *= determinant;
            coeff[intersection.indexInInside()] = vstarhat*normalhat;
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

#endif // DUNE_PDELAB_LOCALOPERATOR_TWOPHASECCFV_HH
