// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONCCFV_HH
#define DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONCCFV_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/typetraits.hh>
#include<dune/geometry/referenceelements.hh>

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
     *                        -(A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_O
     * \f}
     * Note:
     *  - This formulation is valid for velocity fields which are non-divergence free.
     *  - Assumes that the tensor is diagonal !
     *  - Outflow boundary conditions should only be set on the outflow boundary
     *
     * \tparam TP model of ConvectionDiffusionParameterInterface
     */
    template<typename TP>
    class ConvectionDiffusionCCFV
      :
      // public NumericalJacobianSkeleton<ConvectionDiffusionCCFV<TP> >,
      // public NumericalJacobianBoundary<ConvectionDiffusionCCFV<TP> >,
      // public NumericalJacobianVolume<ConvectionDiffusionCCFV<TP> >,
      public FullSkeletonPattern,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      typedef typename ConvectionDiffusionBoundaryConditions::Type BCType;

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

      ConvectionDiffusionCCFV (TP& param_) : param(param_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // cell center
        const Dune::FieldVector<DF,dim>
          inside_local(Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0));

        // evaluate reaction term
        typename TP::Traits::RangeFieldType c = param.c(eg.entity(),inside_local);

        // and accumulate
        r.accumulate(lfsu,0,(c*x(lfsu,0))*eg.geometry().volume());
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // cell center
        const Dune::FieldVector<DF,dim>
          inside_local(Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0));

        // evaluate reaction term
        typename TP::Traits::RangeFieldType c = param.c(eg.entity(),inside_local);

        // and accumulate
        mat.accumulate(lfsu,0,lfsu,0,c*eg.geometry().volume());
      }

      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = IG::dimension;

        // center in face's reference element
        const Dune::FieldVector<DF,dim-1>
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);

        // face volume for integration
        RF face_volume = ig.geometry().integrationElement(face_local)
          *Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).volume();

        auto cell_inside = ig.inside();
        auto cell_outside = ig.outside();

        // cell centers in references elements
        const Dune::FieldVector<DF,dim>
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(cell_inside.type()).position(0,0);
        const Dune::FieldVector<DF,dim>
          outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(cell_outside.type()).position(0,0);

        // evaluate diffusion coefficient from either side and take harmonic average
        typename TP::Traits::PermTensorType tensor_inside;
        tensor_inside = param.A(cell_inside,inside_local);
        typename TP::Traits::PermTensorType tensor_outside;
        tensor_outside = param.A(cell_outside,outside_local);
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F;
        tensor_inside.mv(n_F,An_F);
        RF k_inside = n_F*An_F;
        tensor_outside.mv(n_F,An_F);
        RF k_outside = n_F*An_F;
        RF k_avg = 2.0/(1.0/(k_inside+1E-30) + 1.0/(k_outside+1E-30));

        // evaluate convective term
        Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(face_local);
        typename TP::Traits::RangeType b = param.b(cell_inside,iplocal_s);
        RF vn = b*n_F;
        RF u_upwind=0;
        if (vn>=0) u_upwind = x_s(lfsu_s,0); else u_upwind = x_n(lfsu_n,0);

        // cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension>
          inside_global = cell_inside.geometry().global(inside_local);
        Dune::FieldVector<DF,IG::dimension>
          outside_global = cell_outside.geometry().global(outside_local);

        // distance between the two cell centers
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();

        // contribution to residual on inside element, other residual is computed by symmetric call
        r_s.accumulate(lfsu_s,0,(u_upwind*vn)*face_volume+k_avg*(x_s(lfsu_s,0)-x_n(lfsu_n,0))*face_volume/distance);
        r_n.accumulate(lfsu_n,0,-(u_upwind*vn)*face_volume-k_avg*(x_s(lfsu_s,0)-x_n(lfsu_n,0))*face_volume/distance);
      }

      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                              M& mat_ss, M& mat_sn,
                              M& mat_ns, M& mat_nn) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = IG::dimension;

        // center in face's reference element
        const Dune::FieldVector<DF,dim-1>
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);

        // face volume for integration
        RF face_volume = ig.geometry().integrationElement(face_local)
          *Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).volume();

        auto cell_inside = ig.inside();
        auto cell_outside = ig.outside();

        // cell centers in references elements
        const Dune::FieldVector<DF,dim>
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(cell_inside.type()).position(0,0);
        const Dune::FieldVector<DF,dim>
          outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(cell_outside.type()).position(0,0);

        // evaluate diffusion coefficient from either side and take harmonic average
        typename TP::Traits::PermTensorType tensor_inside;
        tensor_inside = param.A(cell_inside,inside_local);
        typename TP::Traits::PermTensorType tensor_outside;
        tensor_outside = param.A(cell_outside,outside_local);
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F;
        tensor_inside.mv(n_F,An_F);
        RF k_inside = n_F*An_F;
        tensor_outside.mv(n_F,An_F);
        RF k_outside = n_F*An_F;
        RF k_avg = 2.0/(1.0/(k_inside+1E-30) + 1.0/(k_outside+1E-30));

        // evaluate convective term
        Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(face_local);
        typename TP::Traits::RangeType b = param.b(cell_inside,iplocal_s);
        RF vn = b*n_F;

        // cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension>
          inside_global = cell_inside.geometry().global(inside_local);
        Dune::FieldVector<DF,IG::dimension>
          outside_global = cell_outside.geometry().global(outside_local);

        // distance between the two cell centers
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();

        // contribution to residual on inside element, other residual is computed by symmetric call
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
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        const int dim = EG::Geometry::mydimension;

        if (!first_stage) return; // time step calculation is only done in first stage

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // compute optimal dt for this cell
        typename TP::Traits::RangeFieldType cellcapacity = param.d(eg.entity(),inside_local)*eg.geometry().volume();
        typename TP::Traits::RangeFieldType celldt = cellcapacity/(cellinflux+1E-30);
        dtmin = std::min(dtmin,celldt);
      }



      // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = IG::dimension;

        // center in face's reference element
        const Dune::FieldVector<DF,dim-1>
          face_local = Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).position(0,0);

        // face volume for integration
        RF face_volume = ig.geometry().integrationElement(face_local)
          *Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).volume();

        auto cell_inside = ig.inside();

        // cell center in reference element
        const Dune::FieldVector<DF,dim>
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(cell_inside.type()).position(0,0);

        // evaluate boundary condition type
        BCType bctype;
        bctype = param.bctype(ig.intersection(),face_local);

        if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet)
          {
            // Dirichlet boundary
            // distance between cell center and face center
            Dune::FieldVector<DF,dim> inside_global = cell_inside.geometry().global(inside_local);
            Dune::FieldVector<DF,dim> outside_global = ig.geometry().global(face_local);
            inside_global -= outside_global;
            RF distance = inside_global.two_norm();

            // evaluate diffusion coefficient
            typename TP::Traits::PermTensorType tensor_inside;
            tensor_inside = param.A(cell_inside,inside_local);
            const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
            Dune::FieldVector<RF,dim> An_F;
            tensor_inside.mv(n_F,An_F);
            RF k_inside = n_F*An_F;

            // evaluate boundary condition function
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(face_local);
            RF g = param.g(cell_inside,iplocal_s);

            // velocity needed for convection term
            typename TP::Traits::RangeType b = param.b(cell_inside,iplocal_s);
            const Dune::FieldVector<DF,dim> n = ig.centerUnitOuterNormal();

            // contribution to residual on inside element, assumes that Dirichlet boundary is inflow
            r_s.accumulate(lfsu_s,0,(b*n)*g*face_volume + k_inside*(x_s(lfsu_s,0)-g)*face_volume/distance);

            return;
          }

        if (bctype==ConvectionDiffusionBoundaryConditions::Neumann)
          {
            // Neumann boundary
            // evaluate flux boundary condition

            //evaluate boundary function
            typename TP::Traits::RangeFieldType j = param.j(ig.intersection(),face_local);

            // contribution to residual on inside element
            r_s.accumulate(lfsu_s,0,j*face_volume);

            return;
          }

        if (bctype==ConvectionDiffusionBoundaryConditions::Outflow)
          {
            // evaluate velocity field and outer unit normal
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(face_local);
            typename TP::Traits::RangeType b = param.b(cell_inside,iplocal_s);
            const Dune::FieldVector<DF,dim> n = ig.centerUnitOuterNormal();

            // evaluate outflow boundary condition
            typename TP::Traits::RangeFieldType o = param.o(ig.intersection(),face_local);

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
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = IG::dimension;

        // center in face's reference element
        const Dune::FieldVector<DF,dim-1>
          face_local = Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).position(0,0);

        // face volume for integration
        RF face_volume = ig.geometry().integrationElement(face_local)
          *Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).volume();

        auto cell_inside = ig.inside();

        // cell center in reference element
        const Dune::FieldVector<DF,dim>
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(cell_inside.type()).position(0,0);

        // evaluate boundary condition type
        BCType bctype;
        bctype = param.bctype(ig.intersection(),face_local);


        if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet)
          {
            // Dirichlet boundary
            // distance between cell center and face center
            Dune::FieldVector<DF,dim> inside_global = cell_inside.geometry().global(inside_local);
            Dune::FieldVector<DF,dim> outside_global = ig.geometry().global(face_local);
            inside_global -= outside_global;
            RF distance = inside_global.two_norm();

            // evaluate diffusion coefficient
            typename TP::Traits::PermTensorType tensor_inside;
            tensor_inside = param.A(cell_inside,inside_local);
            const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
            Dune::FieldVector<RF,dim> An_F;
            tensor_inside.mv(n_F,An_F);
            RF k_inside = n_F*An_F;

            // contribution to residual on inside element
            mat_ss.accumulate(lfsu_s,0,lfsv_s,0, k_inside*face_volume/distance );

            return;
          }

        if (bctype==ConvectionDiffusionBoundaryConditions::Outflow)
          {
            // evaluate velocity field and outer unit normal
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(face_local);
            typename TP::Traits::RangeType b = param.b(cell_inside,iplocal_s);
            const Dune::FieldVector<DF,dim> n = ig.centerUnitOuterNormal();

            // integrate o
            mat_ss.accumulate(lfsu_s,0,lfsv_s,0, (b*n)*face_volume );

            return;
          }
      }

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        const int dim = EG::Geometry::mydimension;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate source and sink term
        typename TP::Traits::RangeFieldType f = param.f(eg.entity(),inside_local);

        r.accumulate(lfsv,0,-f*eg.geometry().volume());
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
    class ConvectionDiffusionCCFVTemporalOperator
      :
      // public NumericalJacobianApplyVolume<ConvectionDiffusionCCFVTemporalOperator<TP> >,
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
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);
        // capacity term
        typename TP::Traits::RangeFieldType capacity = param.d(eg.entity(),inside_local);

        // residual contribution
        r.accumulate(lfsu,0,capacity*x(lfsu,0)*eg.geometry().volume());
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // capacity term
        typename TP::Traits::RangeFieldType capacity = param.d(eg.entity(),inside_local);

        // residual contribution
        mat.accumulate(lfsu,0,lfsu,0,capacity*eg.geometry().volume());
      }

    private:
      TP& param;
    };


    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
