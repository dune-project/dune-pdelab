// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CONVECTIONDIFFUSIONCCFV_HH
#define DUNE_PDELAB_CONVECTIONDIFFUSIONCCFV_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

namespace Dune {
  namespace PDELab {

    /** a local operator for solving convection-diffusion equation with standard FEM
     *
     * \f{align*}{
     *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \\
     *                                              u &=& g \mbox{ on } \partial\Omega_D \\
     *                (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
     *                        -(A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_O
     * \f}
     * Note:
     *  - This formulation is valid for velocity fields which are non-divergence free.
     *  - Assumes that the tensor is diagonal !
     *  - Outflow boundary conditions should only be set on the outflow boundary
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T>
    class ConvectionDiffusionCCFV : public NumericalJacobianApplySkeleton<ConvectionDiffusionCCFV<T> >,
                                    public NumericalJacobianApplyBoundary<ConvectionDiffusionCCFV<T> >,
                                    public NumericalJacobianApplyVolume<ConvectionDiffusionCCFV<T> >,
    // public NumericalJacobianSkeleton<ConvectionDiffusionCCFV<T> >,
     // public NumericalJacobianBoundary<ConvectionDiffusionCCFV<T> >,
     // public NumericalJacobianVolume<ConvectionDiffusionCCFV<T> >,
                                    public FullSkeletonPattern,
                                    public FullVolumePattern,
                                    public LocalOperatorDefaultFlags,
                                    public InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>

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
      enum { doLambdaVolume   = false };
      enum { doLambdaSkeleton = false };
      enum { doLambdaBoundary = false };

      ConvectionDiffusionCCFV (T& param_) : param(param_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>
          inside_local(Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0));

        // evaluate source and sink term
        typename T::Traits::RangeFieldType f = param.f(eg.entity(),inside_local);
        typename T::Traits::RangeFieldType c = param.c(eg.entity(),inside_local);

        // and accumulate
        r.accumulate(lfsu,0,(c*x(lfsu,0)-f)*eg.geometry().volume());
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>
          inside_local(Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0));

        // evaluate source and sink term
        typename T::Traits::RangeFieldType f = param.f(eg.entity(),inside_local);
        typename T::Traits::RangeFieldType c = param.c(eg.entity(),inside_local);

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

        // cell centers in references elements
        const Dune::FieldVector<DF,dim>
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>
          outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);

        // evaluate diffusion coefficient from either side and take harmonic average
        typename T::Traits::PermTensorType tensor_inside;
        tensor_inside = param.A(*(ig.inside()),inside_local);
        typename T::Traits::PermTensorType tensor_outside;
        tensor_outside = param.A(*(ig.outside()),outside_local);
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F;
        tensor_inside.mv(n_F,An_F);
        RF k_inside = n_F*An_F;
        tensor_outside.mv(n_F,An_F);
        RF k_outside = n_F*An_F;
        RF k_avg = 2.0/(1.0/(k_inside+1E-30) + 1.0/(k_outside+1E-30));

        // evaluate convective term
        Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(face_local);
        typename T::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
        RF vn = b*n_F;
        RF u_upwind=0;
        if (vn>=0) u_upwind = x_s(lfsu_s,0); else u_upwind = x_n(lfsu_n,0);

        // cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension>
          inside_global = ig.inside()->geometry().global(inside_local);
        Dune::FieldVector<DF,IG::dimension>
          outside_global = ig.outside()->geometry().global(outside_local);

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

        // cell centers in references elements
        const Dune::FieldVector<DF,dim>
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>
          outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);

        // evaluate diffusion coefficient from either side and take harmonic average
        typename T::Traits::PermTensorType tensor_inside;
        tensor_inside = param.A(*(ig.inside()),inside_local);
        typename T::Traits::PermTensorType tensor_outside;
        tensor_outside = param.A(*(ig.outside()),outside_local);
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F;
        tensor_inside.mv(n_F,An_F);
        RF k_inside = n_F*An_F;
        tensor_outside.mv(n_F,An_F);
        RF k_outside = n_F*An_F;
        RF k_avg = 2.0/(1.0/(k_inside+1E-30) + 1.0/(k_outside+1E-30));

        // evaluate convective term
        Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(face_local);
        typename T::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
        RF vn = b*n_F;

        // cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension>
          inside_global = ig.inside()->geometry().global(inside_local);
        Dune::FieldVector<DF,IG::dimension>
          outside_global = ig.outside()->geometry().global(outside_local);

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

        // cell center in reference element
        const Dune::FieldVector<DF,dim>
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);

        // evaluate boundary condition type
        BCType bctype;
        bctype = param.bctype(ig.intersection(),face_local);

        if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet)
          {
            // Dirichlet boundary
            // distance between cell center and face center
            Dune::FieldVector<DF,dim> inside_global = ig.inside()->geometry().global(inside_local);
            Dune::FieldVector<DF,dim> outside_global = ig.geometry().global(face_local);
            inside_global -= outside_global;
            RF distance = inside_global.two_norm();

            // evaluate diffusion coefficient
            typename T::Traits::PermTensorType tensor_inside;
            tensor_inside = param.A(*(ig.inside()),inside_local);
            const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
            Dune::FieldVector<RF,dim> An_F;
            tensor_inside.mv(n_F,An_F);
            RF k_inside = n_F*An_F;

            // evaluate boundary condition function
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(face_local);
            RF g = param.g(*(ig.inside()),iplocal_s);

            // velocity needed for convection term
            typename T::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
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
            typename T::Traits::RangeFieldType j = param.j(ig.intersection(),face_local);

            // contribution to residual on inside element
            r_s.accumulate(lfsu_s,0,j*face_volume);

            return;
          }

        if (bctype==ConvectionDiffusionBoundaryConditions::Outflow)
          {
            // evaluate velocity field and outer unit normal
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(face_local);
            typename T::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
            const Dune::FieldVector<DF,dim> n = ig.centerUnitOuterNormal();

            // evaluate outflow boundary condition
            typename T::Traits::RangeFieldType o = param.o(ig.intersection(),face_local);

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

        // cell center in reference element
        const Dune::FieldVector<DF,dim>
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);

        // evaluate boundary condition type
        BCType bctype;
        bctype = param.bctype(ig.intersection(),face_local);


        if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet)
          {
            // Dirichlet boundary
            // distance between cell center and face center
            Dune::FieldVector<DF,dim> inside_global = ig.inside()->geometry().global(inside_local);
            Dune::FieldVector<DF,dim> outside_global = ig.geometry().global(face_local);
            inside_global -= outside_global;
            RF distance = inside_global.two_norm();

            // evaluate diffusion coefficient
            typename T::Traits::PermTensorType tensor_inside;
            tensor_inside = param.A(*(ig.inside()),inside_local);
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
            typename T::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
            const Dune::FieldVector<DF,dim> n = ig.centerUnitOuterNormal();

            // integrate o
            mat_ss.accumulate(lfsu_s,0,lfsv_s,0, (b*n)*face_volume );

            return;
          }
      }

      //! set time in parameter class
      void setTime (double t)
      {
        param.setTime(t);
      }

    private:
      T& param;
    };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
