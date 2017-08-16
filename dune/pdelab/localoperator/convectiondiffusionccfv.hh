// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONCCFV_HH
#define DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONCCFV_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/typetraits.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

namespace Dune {
  namespace PDELab {

    /** a local operator for solving the linear convection-diffusion equation with CCFV
     *
     * \f{align*}{
     *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \\
     *                                         u(t,x) &=& g(t,x) \mbox{ on } \partial\Omega_D \\
     *                (b(x) u - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
     *                        -(A(x)\nabla u) \cdot n &=& o \mbox{ on } \partial\Omega_O
     * \f}
     * Note:
     *  - This formulation is valid for velocity fields which are non-divergence free.
     *  - Assumes that the tensor is diagonal !
     *  - Outflow boundary conditions should only be set on the outflow boundary
     *
     * \tparam TP model of ConvectionDiffusionParameterInterface
     */
    template<typename TP>
    class ConvectionDiffusionCCFV :
      public Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionCCFV<TP> >,
      public FullSkeletonPattern,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      using BCType = typename ConvectionDiffusionBoundaryConditions::Type;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume    = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume   = true };
      enum { doLambdaSkeleton = false };
      enum { doLambdaBoundary = false };

      ConvectionDiffusionCCFV (TP& param_)
        : Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionCCFV<TP> >(1.0e-7)
        , param(param_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // get cell
        const auto& cell = eg.entity();

        // cell center
        auto geo = eg.geometry();
        auto ref_el = referenceElement(geo);
        auto local_inside = ref_el.position(0,0);

        // evaluate reaction term
        auto c = param.c(cell,local_inside);

        // and accumulate
        r.accumulate(lfsu,0,(c*x(lfsu,0))*geo.volume());
      }

      // apply jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
      {
        alpha_volume(eg,lfsu,x,lfsv,y);
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // get cell
        const auto& cell = eg.entity();

        // cell center
        auto geo = eg.geometry();
        auto ref_el = referenceElement(geo);
        auto local_inside = ref_el.position(0,0);

        // evaluate reaction term
        auto c = param.c(cell,local_inside);

        // and accumulate
        mat.accumulate(lfsu,0,lfsu,0,c*geo.volume());
      }

      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;

        // dimensions
        const auto dim = IG::dimension;

        // get cell entities from both sides of the intersection
        auto cell_inside = ig.inside();
        auto cell_outside = ig.outside();

        // get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();
        auto geo_outside = cell_outside.geometry();

        // get geometry of intersection in local coordinates of neighbor cells
        auto geo_in_inside = ig.geometryInInside();

        // center in face's reference element
        auto ref_el = referenceElement(geo);
        auto face_local = ref_el.position(0,0);

        // face volume for integration
        auto face_volume = geo.integrationElement(face_local) * ref_el.volume();

        // cell centers in references elements
        auto ref_el_inside = referenceElement(geo_inside);
        auto ref_el_outside = referenceElement(geo_outside);
        auto local_inside = ref_el_inside.position(0,0);
        auto local_outside = ref_el_outside.position(0,0);

        // evaluate diffusion coefficient from either side and take harmonic average
        auto tensor_inside = param.A(cell_inside,local_inside);
        auto tensor_outside = param.A(cell_outside,local_outside);
        auto n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F;
        tensor_inside.mv(n_F,An_F);
        auto k_inside = n_F*An_F;
        tensor_outside.mv(n_F,An_F);
        auto k_outside = n_F*An_F;
        auto k_avg = 2.0/(1.0/(k_inside+1E-30) + 1.0/(k_outside+1E-30));

        // evaluate convective term
        auto iplocal_s = geo_in_inside.global(face_local);
        auto b = param.b(cell_inside,iplocal_s);
        auto vn = b*n_F;
        RF u_upwind=0;
        if (vn>=0) u_upwind = x_s(lfsu_s,0); else u_upwind = x_n(lfsu_n,0);

        // cell centers in global coordinates
        auto global_inside = geo_inside.global(local_inside);
        auto global_outside = geo_outside.global(local_outside);

        // distance between the two cell centers
        global_inside -= global_outside;
        auto distance = global_inside.two_norm();

        // contribution to residual on inside element, other residual is computed by symmetric call
        r_s.accumulate(lfsu_s,0,(u_upwind*vn)*face_volume+k_avg*(x_s(lfsu_s,0)-x_n(lfsu_n,0))*face_volume/distance);
        r_n.accumulate(lfsu_n,0,-(u_upwind*vn)*face_volume-k_avg*(x_s(lfsu_s,0)-x_n(lfsu_n,0))*face_volume/distance);
      }

      // apply jacobian of skeleton term
      template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_skeleton (const IG& ig,
                                    const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                    const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                                    Y& y_s, Y& y_n) const
      {
        alpha_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,y_s,y_n);
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                              M& mat_ss, M& mat_sn,
                              M& mat_ns, M& mat_nn) const
      {
        // define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;

        // dimensions
        const auto dim = IG::dimension;

        // get cell entities from both sides of the intersection
        auto cell_inside = ig.inside();
        auto cell_outside = ig.outside();

        // get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();
        auto geo_outside = cell_outside.geometry();

        // get geometry of intersection in local coordinates of neighbor cells
        auto geo_in_inside = ig.geometryInInside();

        // center in face's reference element
        auto ref_el = referenceElement(geo);
        auto face_local = ref_el.position(0,0);

        // face volume for integration
        auto face_volume = geo.integrationElement(face_local) * ref_el.volume();

        // cell centers in references elements
        auto ref_el_inside = referenceElement(geo_inside);
        auto ref_el_outside = referenceElement(geo_outside);
        auto local_inside = ref_el_inside.position(0,0);
        auto local_outside = ref_el_outside.position(0,0);

        // evaluate diffusion coefficient from either side and take harmonic average
        auto tensor_inside = param.A(cell_inside,local_inside);
        auto tensor_outside = param.A(cell_outside,local_outside);
        auto n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F;
        tensor_inside.mv(n_F,An_F);
        auto k_inside = n_F*An_F;
        tensor_outside.mv(n_F,An_F);
        auto k_outside = n_F*An_F;
        auto k_avg = 2.0/(1.0/(k_inside+1E-30) + 1.0/(k_outside+1E-30));

        // evaluate convective term
        auto iplocal_s = geo_in_inside.global(face_local);
        auto b = param.b(cell_inside,iplocal_s);
        auto vn = b*n_F;

        // cell centers in global coordinates
        auto global_inside = geo_inside.global(local_inside);
        auto global_outside = geo_outside.global(local_outside);

        // distance between the two cell centers
        global_inside -= global_outside;
        auto distance = global_inside.two_norm();

        // contribution to jacobians on inside element and outside element for test and trial function
        mat_ss.accumulate(lfsu_s,0,lfsu_s,0,   k_avg*face_volume/distance );
        mat_ns.accumulate(lfsu_n,0,lfsu_s,0,  -k_avg*face_volume/distance );
        mat_sn.accumulate(lfsu_s,0,lfsu_n,0,  -k_avg*face_volume/distance );
        mat_nn.accumulate(lfsu_n,0,lfsu_n,0,   k_avg*face_volume/distance );
        if (vn>=0)
          {
            mat_ss.accumulate(lfsu_s,0,lfsu_s,0,   vn*face_volume );
            mat_ns.accumulate(lfsu_n,0,lfsu_s,0,  -vn*face_volume );
          }
        else
          {
            mat_sn.accumulate(lfsu_s,0,lfsu_n,0,   vn*face_volume );
            mat_nn.accumulate(lfsu_n,0,lfsu_n,0,  -vn*face_volume );
          }
      }


      // post skeleton: compute time step allowable for cell; to be done later
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                      const LFSV& lfsv, R& r) const
      {
        if (not first_stage) return; // time step calculation is only done in first stage

        // get cell
        const auto& cell = eg.entity();

        // cell center
        auto geo = eg.geometry();
        auto ref_el = referenceElement(geo);
        auto local_inside = ref_el.position(0,0);

        // compute optimal dt for this cell
        auto cellcapacity = param.d(cell,local_inside)*geo.volume();
        auto celldt = cellcapacity/(cellinflux+1E-30);
        dtmin = std::min(dtmin,celldt);
      }


      // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;

        // dimensions
        const auto dim = IG::dimension;

        // get cell entities from both sides of the intersection
        auto cell_inside = ig.inside();

        // get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // get geometry of intersection in local coordinates of neighbor cells
        auto geo_in_inside = ig.geometryInInside();

        // center in face's reference element
        auto ref_el = referenceElement(geo);
        auto face_local = ref_el.position(0,0);

        // face volume for integration
        auto face_volume = geo.integrationElement(face_local) * ref_el.volume();

        // cell centers in references elements
        auto ref_el_inside = referenceElement(geo_inside);
        auto local_inside = ref_el_inside.position(0,0);

        // evaluate boundary condition type
        auto bctype = param.bctype(ig.intersection(),face_local);

        if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet)
          {
            // Dirichlet boundary
            // distance between cell center and face center
            auto global_inside = geo_inside.global(local_inside);
            auto global_outside = geo.global(face_local);
            global_inside -= global_outside;
            auto distance = global_inside.two_norm();

            // evaluate diffusion coefficient
            auto tensor_inside = param.A(cell_inside,local_inside);
            auto n_F = ig.centerUnitOuterNormal();
            Dune::FieldVector<RF,dim> An_F;
            tensor_inside.mv(n_F,An_F);
            auto k_inside = n_F*An_F;

            // evaluate boundary condition function
            auto iplocal_s = geo_in_inside.global(face_local);
            auto g = param.g(cell_inside,iplocal_s);

            // velocity needed for convection term
            auto b = param.b(cell_inside,iplocal_s);
            auto n = ig.centerUnitOuterNormal();

            // contribution to residual on inside element, assumes that Dirichlet boundary is inflow
            r_s.accumulate(lfsu_s,0,(b*n)*g*face_volume + k_inside*(x_s(lfsu_s,0)-g)*face_volume/distance);

            return;
          }

        if (bctype==ConvectionDiffusionBoundaryConditions::Neumann)
          {
            // Neumann boundary
            // evaluate flux boundary condition

            //evaluate boundary function
            auto j = param.j(ig.intersection(),face_local);

            // contribution to residual on inside element
            r_s.accumulate(lfsu_s,0,j*face_volume);

            return;
          }

        if (bctype==ConvectionDiffusionBoundaryConditions::Outflow)
          {
            // evaluate velocity field and outer unit normal
            auto iplocal_s = geo_in_inside.global(face_local);
            auto b = param.b(cell_inside,iplocal_s);
            auto n = ig.centerUnitOuterNormal();

            // evaluate outflow boundary condition
            auto o = param.o(ig.intersection(),face_local);

            // integrate o
            r_s.accumulate(lfsu_s,0,((b*n)*x_s(lfsu_s,0) + o)*face_volume);

            return;
          }
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              M& mat_ss) const
      {
        // define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;

        // dimensions
        const auto dim = IG::dimension;

        // get cell entities from both sides of the intersection
        auto cell_inside = ig.inside();

        // get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // get geometry of intersection in local coordinates of neighbor cells
        auto geo_in_inside = ig.geometryInInside();

        // center in face's reference element
        auto ref_el = referenceElement(geo);
        auto face_local = ref_el.position(0,0);

        // face volume for integration
        auto face_volume = geo.integrationElement(face_local) * ref_el.volume();

        // cell centers in references elements
        auto ref_el_inside = referenceElement(geo_inside);
        auto local_inside = ref_el_inside.position(0,0);

        // evaluate boundary condition type
        auto bctype = param.bctype(ig.intersection(),face_local);

        if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet)
          {
            // Dirichlet boundary
            // distance between cell center and face center
            auto global_inside = geo_inside.global(local_inside);
            auto global_outside = geo.global(face_local);
            global_inside -= global_outside;
            auto distance = global_inside.two_norm();

            // evaluate diffusion coefficient
            auto tensor_inside = param.A(cell_inside,local_inside);
            auto n_F = ig.centerUnitOuterNormal();
            Dune::FieldVector<RF,dim> An_F;
            tensor_inside.mv(n_F,An_F);
            auto k_inside = n_F*An_F;

            // contribution to jacobian on inside element for test and trial function
            mat_ss.accumulate(lfsu_s,0,lfsv_s,0, k_inside*face_volume/distance );

            return;
          }

        if (bctype==ConvectionDiffusionBoundaryConditions::Outflow)
          {
            // evaluate velocity field and outer unit normal
            auto iplocal_s = geo_in_inside.global(face_local);
            auto b = param.b(cell_inside,iplocal_s);
            auto n = ig.centerUnitOuterNormal();

            // integrate o
            mat_ss.accumulate(lfsu_s,0,lfsv_s,0, (b*n)*face_volume );

            return;
          }
      }

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // get cell
        const auto& cell = eg.entity();

        // cell center
        auto geo = eg.geometry();
        auto ref_el = referenceElement(geo);
        auto local_inside = ref_el.position(0,0);

        // evaluate source and sink term
        auto f = param.f(cell,local_inside);

        r.accumulate(lfsv,0,-f*geo.volume());
      }

      //! set time in parameter class
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        param.setTime(t);
      }

      //! to be called once before each time step
      void preStep (typename TP::Traits::RangeFieldType time, typename TP::Traits::RangeFieldType dt,
                    int stages)
      {
      }

      //! to be called once before each stage
      void preStage (typename TP::Traits::RangeFieldType time, int r)
      {
        if (r==1)
          {
            first_stage = true;
            dtmin = 1E100;
          }
        else first_stage = false;
      }

      //! to be called once at the end of each stage
      void postStage ()
      {
      }

      //! to be called once before each stage
      typename TP::Traits::RangeFieldType suggestTimestep (typename TP::Traits::RangeFieldType dt) const
      {
        return dtmin;
      }

    private:
      TP& param;
      bool first_stage;
      mutable typename TP::Traits::RangeFieldType dtmin; // accumulate minimum dt here
      mutable typename TP::Traits::RangeFieldType cellinflux;
    };




    /** a local operator for the weighted mass matrix (capacity term d(x))
     *
     * \f{align*}{
         \int_\Omega d(x) uv dx
     * \f}
     */
    template<class TP>
    class ConvectionDiffusionCCFVTemporalOperator :
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      ConvectionDiffusionCCFVTemporalOperator (TP& param_)
        : param(param_)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // get cell
        const auto& cell = eg.entity();

        // cell center
        auto geo = eg.geometry();
        auto ref_el = referenceElement(geo);
        auto local_inside = ref_el.position(0,0);

        // capacity term
        auto capacity = param.d(cell,local_inside);

        // residual contribution
        r.accumulate(lfsu,0,capacity*x(lfsu,0)*geo.volume());
      }

      // apply jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
      {
        alpha_volume(eg,lfsu,x,lfsv,y);
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // get cell
        const auto& cell = eg.entity();

        // cell center
        auto geo = eg.geometry();
        auto ref_el = referenceElement(geo);
        auto local_inside = ref_el.position(0,0);

        // capacity term
        auto capacity = param.d(cell,local_inside);

        // jacobian contribution
        mat.accumulate(lfsu,0,lfsu,0,capacity*geo.volume());
      }

    private:
      TP& param;
    };


    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONCCFV_HH
