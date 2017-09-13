// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONFEM_HH
#define DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONFEM_HH

#include<vector>

#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

#include"convectiondiffusionparameter.hh"

namespace Dune {
  namespace PDELab {

    /** a local operator for solving the linear convection-diffusion equation with standard FEM
     *
     * \f{align*}{
     *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \\
     *                                              u &=& g \mbox{ on } \partial\Omega_D \\
     *                (b(x) u - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
     *                        -(A(x)\nabla u) \cdot n &=& o \mbox{ on } \partial\Omega_O
     * \f}
     * Note:
     *  - This formulation is valid for velocity fields which are non-divergence free.
     *  - Outflow boundary conditions should only be set on the outflow boundary
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T, typename FiniteElementMap>
    class ConvectionDiffusionFEM :
      public Dune::PDELab::NumericalJacobianApplyVolume<ConvectionDiffusionFEM<T,FiniteElementMap> >,
      public Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionFEM<T,FiniteElementMap> >,
      public Dune::PDELab::NumericalJacobianVolume<ConvectionDiffusionFEM<T,FiniteElementMap> >,
      public Dune::PDELab::NumericalJacobianBoundary<ConvectionDiffusionFEM<T,FiniteElementMap> >,
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags,
      public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doAlphaBoundary = true };

      ConvectionDiffusionFEM (T& param_, int intorderadd_=0)
        : param(param_), intorderadd(intorderadd_)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // Define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSU::Traits::SizeType;

        // dimensions
        const int dim = EG::Entity::dimension;

        // Reference to cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        auto ref_el = referenceElement(geo);
        auto localcenter = ref_el.position(0,0);
        auto tensor = param.A(cell,localcenter);

        // Initialize vectors outside for loop
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
        Dune::FieldVector<RF,dim> gradu(0.0);
        Dune::FieldVector<RF,dim> Agradu(0.0);

        // Transformation matrix
        typename EG::Geometry::JacobianInverseTransposed jac;

        // loop over quadrature points
        auto intorder = intorderadd+2*lfsu.finiteElement().localBasis().order();
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            auto& phi = cache.evaluateFunction(ip.position(),lfsu.finiteElement().localBasis());

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            auto& js = cache.evaluateJacobian(ip.position(),lfsu.finiteElement().localBasis());

            // transform gradients of shape functions to real element
            jac = geo.jacobianInverseTransposed(ip.position());
            for (size_type i=0; i<lfsu.size(); i++)
              jac.mv(js[i][0],gradphi[i]);

            // compute gradient of u
            gradu = 0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              gradu.axpy(x(lfsu,i),gradphi[i]);

            // compute A * gradient of u
            tensor.mv(gradu,Agradu);

            // evaluate velocity field, sink term and source term
            auto b = param.b(cell,ip.position());
            auto c = param.c(cell,ip.position());
            auto f = param.f(cell,ip.position());

            // integrate (A grad u)*grad phi_i - u b*grad phi_i + c*u*phi_i
            RF factor = ip.weight() * geo.integrationElement(ip.position());
            for (size_type i=0; i<lfsu.size(); i++)
              r.accumulate(lfsu,i,( Agradu*gradphi[i] - u*(b*gradphi[i]) + (c*u-f)*phi[i] )*factor);
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // Define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSU::Traits::SizeType;

        // dimensions
        const int dim = EG::Entity::dimension;

        // Reference to cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        auto ref_el = referenceElement(geo);
        auto localcenter = ref_el.position(0,0);
        auto tensor = param.A(cell,localcenter);

        // Initialize vectors outside for loop
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
        std::vector<Dune::FieldVector<RF,dim> > Agradphi(lfsu.size());

        // Transformation matrix
        typename EG::Geometry::JacobianInverseTransposed jac;

        // loop over quadrature points
        auto intorder = intorderadd+2*lfsu.finiteElement().localBasis().order();
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            auto& js = cache.evaluateJacobian(ip.position(),lfsu.finiteElement().localBasis());

            // transform gradient to real element
            jac = geo.jacobianInverseTransposed(ip.position());
            for (size_type i=0; i<lfsu.size(); i++)
              {
                jac.mv(js[i][0],gradphi[i]);
                tensor.mv(gradphi[i],Agradphi[i]);
              }

            // evaluate basis functions
            auto& phi = cache.evaluateFunction(ip.position(),lfsu.finiteElement().localBasis());

            // evaluate velocity field, sink term and source term
            auto b = param.b(cell,ip.position());
            auto c = param.c(cell,ip.position());

            // integrate (A grad phi_j)*grad phi_i - phi_j b*grad phi_i + c*phi_j*phi_i
            RF factor = ip.weight() * geo.integrationElement(ip.position());
            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat.accumulate(lfsu,i,lfsu,j,( Agradphi[j]*gradphi[i]-phi[j]*(b*gradphi[i])+c*phi[j]*phi[i] )*factor);
          }
      }

      // boundary integral
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // Define types
        using RF = typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSV::Traits::SizeType;

        // Reference to the inside cell
        const auto& cell_inside = ig.inside();

        // Get geometry
        auto geo = ig.geometry();

        // Get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // evaluate boundary condition type
        auto ref_el = referenceElement(geo_in_inside);
        auto local_face_center = ref_el.position(0,0);
        auto intersection = ig.intersection();
        auto bctype = param.bctype(intersection,local_face_center);

        // skip rest if we are on Dirichlet boundary
        if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet) return;

        // loop over quadrature points and integrate normal flux
        auto intorder = intorderadd+2*lfsu_s.finiteElement().localBasis().order();
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // position of quadrature point in local coordinates of element
            auto local = geo_in_inside.global(ip.position());

            // evaluate shape functions (assume Galerkin method)
            auto& phi = cache.evaluateFunction(local,lfsu_s.finiteElement().localBasis());

            if (bctype==ConvectionDiffusionBoundaryConditions::Neumann)
              {
                // evaluate flux boundary condition
                auto j = param.j(intersection,ip.position());

                // integrate j
                auto factor = ip.weight()*geo.integrationElement(ip.position());
                for (size_type i=0; i<lfsu_s.size(); i++)
                  r_s.accumulate(lfsu_s,i,j*phi[i]*factor);
              }

            if (bctype==ConvectionDiffusionBoundaryConditions::Outflow)
              {
                // evaluate u
                RF u=0.0;
                for (size_type i=0; i<lfsu_s.size(); i++)
                  u += x_s(lfsu_s,i)*phi[i];

                // evaluate velocity field and outer unit normal
                auto b = param.b(cell_inside,local);
                auto n = ig.unitOuterNormal(ip.position());

                // evaluate outflow boundary condition
                auto o = param.o(intersection,ip.position());

                // integrate o
                auto factor = ip.weight()*geo.integrationElement(ip.position());
                for (size_type i=0; i<lfsu_s.size(); i++)
                  r_s.accumulate(lfsu_s,i,( (b*n)*u + o)*phi[i]*factor);
              }
          }
      }

      // jacobian contribution from boundary
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              M& mat_s) const
      {
        // Define types
        using size_type = typename LFSV::Traits::SizeType;

        // Reference to the inside cell
        const auto& cell_inside = ig.inside();

        // Get geometry
        auto geo = ig.geometry();

        // Get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // evaluate boundary condition type
        auto ref_el = referenceElement(geo_in_inside);
        auto local_face_center = ref_el.position(0,0);
        auto intersection = ig.intersection();
        auto bctype = param.bctype(intersection,local_face_center);

        // skip rest if we are on Dirichlet or Neumann boundary
        if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet) return;
        if (bctype==ConvectionDiffusionBoundaryConditions::Neumann) return;

        // loop over quadrature points and integrate normal flux
        auto intorder = intorderadd+2*lfsu_s.finiteElement().localBasis().order();
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // position of quadrature point in local coordinates of element
            auto local = geo_in_inside.global(ip.position());

            // evaluate shape functions (assume Galerkin method)
            auto& phi = cache.evaluateFunction(local,lfsu_s.finiteElement().localBasis());

            // evaluate velocity field and outer unit normal
            auto b = param.b(cell_inside,local);
            auto n = ig.unitOuterNormal(ip.position());

            // integrate
            auto factor = ip.weight()*geo.integrationElement(ip.position());
            for (size_type j=0; j<lfsu_s.size(); j++)
              for (size_type i=0; i<lfsu_s.size(); i++)
                mat_s.accumulate(lfsu_s,i,lfsu_s,j,(b*n)*phi[j]*phi[i]*factor);
          }
      }


      //! set time in parameter class
      void setTime (double t)
      {
        param.setTime(t);
      }

    private:
      T& param;
      int intorderadd;
      using LocalBasisType = typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType;
      Dune::PDELab::LocalBasisCache<LocalBasisType> cache;
    };



    /** a local operator for residual-based error estimation
     *
     * A call to residual() of a grid operator space will assemble
     * the quantity \f$\eta_T^2\f$ for each cell. Note that the squares
     * of the cell indicator \f$\eta_T\f$ is stored. To compute the global
     * error estimate sum up all values and take the square root.
     *
     * Assumptions and limitations:
     * - Assumes that LFSU is \f$P_k\f$/\f$Q_k\f$ finite element space
     *   and LFSV is a \f$P_0\f$ finite element space (one value per cell).
     *   However, the second order derivatives are ignored!
     * - Convection term is ignored (but reaction term is included)
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T>
    class ConvectionDiffusionFEMResidualEstimator
      : public Dune::PDELab::LocalOperatorDefaultFlags
    {
      enum { dim = T::Traits::GridViewType::dimension };

      using Real = typename T::Traits::RangeFieldType;
      using BCType = typename ConvectionDiffusionBoundaryConditions::Type;

    public:
      // pattern assembly flags
      enum { doPatternVolume = false };
      enum { doPatternSkeleton = false };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };

      //! constructor: pass parameter object
      ConvectionDiffusionFEMResidualEstimator (T& param_)
        : param(param_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // Define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using RangeType = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using size_type = typename LFSU::Traits::SizeType;

        // Reference to cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // Initialize vectors outside for loop
        std::vector<RangeType> phi(lfsu.size());

        // loop over quadrature points
        RF sum(0.0);
        auto intorder = 2*lfsu.finiteElement().localBasis().order();
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            lfsu.finiteElement().localBasis().evaluateFunction(ip.position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // evaluate reaction term
            auto c = param.c(cell,ip.position());

            // evaluate right hand side parameter function
            auto f = param.f(cell,ip.position());

            // integrate f^2
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            sum += (f*f-c*c*u*u)*factor;
          }

        // accumulate cell indicator
        auto h_T = diameter(geo);
        r.accumulate(lfsv,0,h_T*h_T*sum);
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
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using JacobianType = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType;
        using size_type = typename LFSU::Traits::SizeType;

        // dimensions
        const int dim = IG::dimension;

        // References to inside and outside cells
        const auto& cell_inside = ig.inside();
        const auto& cell_outside = ig.outside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();
        auto geo_outside = cell_outside.geometry();

        // Get geometry of intersection in local coordinates of cell_inside and cell_outside
        auto geo_in_inside = ig.geometryInInside();
        auto geo_in_outside = ig.geometryInOutside();

        // evaluate permeability tensors
        auto ref_el_inside = referenceElement(geo_inside);
        auto ref_el_outside = referenceElement(geo_outside);
        auto local_inside = ref_el_inside.position(0,0);
        auto local_outside = ref_el_outside.position(0,0);
        auto A_s = param.A(cell_inside,local_inside);
        auto A_n = param.A(cell_outside,local_outside);

        // tensor times normal
        auto n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        Dune::FieldVector<RF,dim> An_F_n;
        A_n.mv(n_F,An_F_n);

        // Initialize vectors outside for loop
        std::vector<JacobianType> gradphi_s(lfsu_s.size());
        std::vector<JacobianType> gradphi_n(lfsu_n.size());
        std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
        std::vector<Dune::FieldVector<RF,dim> > tgradphi_n(lfsu_n.size());
        Dune::FieldVector<RF,dim> gradu_s(0.0);
        Dune::FieldVector<RF,dim> gradu_n(0.0);

        // Transformation matrix
        typename IG::Entity::Geometry::JacobianInverseTransposed jac;

        // loop over quadrature points and integrate normal flux
        RF sum(0.0);
        auto intorder = 2*lfsu_s.finiteElement().localBasis().order();
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // position of quadrature point in local coordinates of elements
            auto iplocal_s = geo_in_inside.global(ip.position());
            auto iplocal_n = geo_in_outside.global(ip.position());

            // evaluate gradient of basis functions
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
            lfsu_n.finiteElement().localBasis().evaluateJacobian(iplocal_n,gradphi_n);

            // transform gradients of shape functions to real element
            jac = geo_inside.jacobianInverseTransposed(iplocal_s);
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
            jac = geo_outside.jacobianInverseTransposed(iplocal_n);
            for (size_type i=0; i<lfsu_n.size(); i++) jac.mv(gradphi_n[i][0],tgradphi_n[i]);

            // compute gradient of u
            gradu_s = 0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);
            gradu_n = 0.0;
            for (size_type i=0; i<lfsu_n.size(); i++)
              gradu_n.axpy(x_n(lfsu_n,i),tgradphi_n[i]);

            // integrate
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            auto jump = (An_F_s*gradu_s)-(An_F_n*gradu_n);
            sum += 0.25*jump*jump*factor;
          }

        // accumulate indicator
        // DF h_T = diameter(ig.geometry());
        auto h_T = std::max(diameter(geo_inside),diameter(geo_outside));
        r_s.accumulate(lfsv_s,0,h_T*sum);
        r_n.accumulate(lfsv_n,0,h_T*sum);
      }


      // boundary integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // Define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using JacobianType = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType;
        using size_type = typename LFSU::Traits::SizeType;

        // dimensions
        const int dim = IG::dimension;

        // References to the inside cell
        const auto& cell_inside = ig.inside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // evaluate permeability tensors
        auto ref_el_inside = referenceElement(geo_inside);
        auto local_inside = ref_el_inside.position(0,0);
        auto A_s = param.A(cell_inside,local_inside);

        // tensor times normal
        auto n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);

        // get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // evaluate boundary condition
        auto ref_el_in_inside = referenceElement(geo_in_inside);
        auto face_local = ref_el_in_inside.position(0,0);
        auto bctype = param.bctype(ig.intersection(),face_local);
        if (bctype != ConvectionDiffusionBoundaryConditions::Neumann)
          return;

        // Initialize vectors outside for loop
        std::vector<JacobianType> gradphi_s(lfsu_s.size());
        std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
        Dune::FieldVector<RF,dim> gradu_s(0.0);

        // Transformation matrix
        typename IG::Entity::Geometry::JacobianInverseTransposed jac;

        // loop over quadrature points and integrate normal flux
        RF sum(0.0);
        auto intorder = 2*lfsu_s.finiteElement().localBasis().order();
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // position of quadrature point in local coordinates of elements
            auto iplocal_s = geo_in_inside.global(ip.position());

            // evaluate gradient of basis functions
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);

            // transform gradients of shape functions to real element
            jac = geo_inside.jacobianInverseTransposed(iplocal_s);
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);

            // compute gradient of u
            gradu_s = 0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);

            // evaluate flux boundary condition
            auto j = param.j(ig.intersection(),ip.position());

            // integrate
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            auto jump = j+(An_F_s*gradu_s);
            sum += jump*jump*factor;
          }

        // accumulate indicator
        //DF h_T = diameter(ig.geometry());
        auto h_T = diameter(geo_inside);
        r_s.accumulate(lfsv_s,0,h_T*sum);
      }

    private:
      T& param;  // two phase parameter class

      template<class GEO>
      typename GEO::ctype diameter (const GEO& geo) const
      {
        using DF = typename GEO::ctype;
        DF hmax = -1.0E00;
        const int dim = GEO::coorddimension;
        for (int i=0; i<geo.corners(); i++)
          {
            Dune::FieldVector<DF,dim> xi = geo.corner(i);
            for (int j=i+1; j<geo.corners(); j++)
              {
                Dune::FieldVector<DF,dim> xj = geo.corner(j);
                xj -= xi;
                hmax = std::max(hmax,xj.two_norm());
              }
          }
        return hmax;
      }

    };


    /** a local operator for evaluating the temporal part of error estimator
     *
     * A call to residual() of a grid operator space will assemble
     * the quantity \f$\eta_T^2\f$ for each cell. Note that the squares
     * of the cell indicator \f$\eta_T\f$ is stored. To compute the global
     * error estimate sum up all values and take the square root.
     *
     * Assumptions and limitations:
     * - Assumes that LFSU is \f$P_k\f$/\f$Q_k\f$ finite element space
     *   and LFSV is a \f$P_0\f$ finite element space (one value per cell).
     * - Assumes that x is the jump from one time interval to the next, i.e. x=xnew-xold.
     * - Data oscillation part is currently not yet implemented
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T>
    class ConvectionDiffusionTemporalResidualEstimator1
      : public Dune::PDELab::LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<double>
    {
      enum { dim = T::Traits::GridViewType::dimension };

      using Real = typename T::Traits::RangeFieldType;
      using BCType = typename ConvectionDiffusionBoundaryConditions::Type;

    public:
      // pattern assembly flags
      enum { doPatternVolume = false };
      enum { doPatternSkeleton = false };

      // residual assembly flags
      enum { doAlphaVolume  = true };

      //! constructor: pass parameter object
      // supply time step from implicit Euler scheme
      ConvectionDiffusionTemporalResidualEstimator1 (T& param_, double time_, double dt_)
        : param(param_), time(time_), dt(dt_), cmax(0)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // Define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using RangeType = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using size_type = typename LFSU::Traits::SizeType;

        // Reference to the cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // Initialize vectors outside for loop
        std::vector<RangeType> phi(lfsu.size());

        // loop over quadrature points
        RF sum(0.0);
        RF fsum_up(0.0);
        RF fsum_mid(0.0);
        RF fsum_down(0.0);
        auto intorder = 2*lfsu.finiteElement().localBasis().order();
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            lfsu.finiteElement().localBasis().evaluateFunction(ip.position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // integrate f^2
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            sum += u*u*factor;

            // evaluate right hand side parameter function
            param.setTime(time);
            auto f_down = param.f(cell,ip.position());
            param.setTime(time+0.5*dt);
            auto f_mid = param.f(cell,ip.position());
            param.setTime(time+dt);
            auto f_up = param.f(cell,ip.position());
            auto f_average = (1.0/6.0)*f_down + (2.0/3.0)*f_mid + (1.0/6.0)*f_up;

            // integrate f-f_average
            fsum_down += (f_down-f_average)*(f_down-f_average)*factor;
            fsum_mid += (f_mid-f_average)*(f_mid-f_average)*factor;
            fsum_up += (f_up-f_average)*(f_up-f_average)*factor;
          }

        // accumulate cell indicator
        auto h_T = diameter(geo);
        r.accumulate(lfsv,0,(h_T*h_T/dt)*sum); // h^2*k_n||jump/k_n||^2
        r.accumulate(lfsv,0,h_T*h_T * dt*((1.0/6.0)*fsum_down+(2.0/3.0)*fsum_mid+(1.0/6.0)*fsum_up) ); // h^2*||f-time_average(f)||^2_0_s_t
      }

      void clearCmax ()
      {
        cmax = 0;
      }

      double getCmax () const
      {
        return cmax;
      }

    private:
      T& param;  // two phase parameter class
      double time;
      double dt;
      mutable double cmax;

      template<class GEO>
      typename GEO::ctype diameter (const GEO& geo) const
      {
        using DF = typename GEO::ctype;
        DF hmax = -1.0E00;
        const int dim = GEO::coorddimension;
        for (int i=0; i<geo.corners(); i++)
          {
            Dune::FieldVector<DF,dim> xi = geo.corner(i);
            for (int j=i+1; j<geo.corners(); j++)
              {
                Dune::FieldVector<DF,dim> xj = geo.corner(j);
                xj -= xi;
                hmax = std::max(hmax,xj.two_norm());
              }
          }
        return hmax;
      }

    };

    // a functor that can be used to evaluate rhs parameter function in interpolate
    template<typename T, typename EG>
    class CD_RHS_LocalAdapter
    {
    public:
      CD_RHS_LocalAdapter (const T& t_, const EG& eg_) : t(t_), eg(eg_)
      {}

      template<typename X, typename Y>
      inline void evaluate (const X& x, Y& y) const
      {
        y[0] = t.f(eg.entity(),x);
      }

    private:
      const T& t;
      const EG& eg;
    };

    /** a local operator for evaluating the temporal part of error estimator
     *
     * A call to residual() of a grid operator space will assemble
     * the quantity \f$\eta_T^2\f$ for each cell. Note that the squares
     * of the cell indicator \f$\eta_T\f$ is stored. To compute the global
     * error estimate sum up all values and take the square root.
     *
     * Assumptions and limitations:
     * - Assumes that LFSU is \f$P_k\f$/\f$Q_k\f$ finite element space
     *   and LFSV is a \f$P_0\f$ finite element space (one value per cell).
     * - Assumes that x is the jump from one time interval to the next, i.e. x=xnew-xold.
     * - Neumann boundary part not completely implemented
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T>
    class ConvectionDiffusionTemporalResidualEstimator
      : public Dune::PDELab::LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<double>
    {
      enum { dim = T::Traits::GridViewType::dimension };

      using Real = typename T::Traits::RangeFieldType;
      using BCType = typename ConvectionDiffusionBoundaryConditions::Type;

    public:
      // pattern assembly flags
      enum { doPatternVolume = false };
      enum { doPatternSkeleton = false };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaBoundary  = true };

      //! constructor: pass parameter object
      // supply time step from implicit Euler scheme
      ConvectionDiffusionTemporalResidualEstimator (T& param_, double time_, double dt_)
        : param(param_), time(time_), dt(dt_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // Define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using RangeType = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using size_type = typename LFSU::Traits::SizeType;
        using JacobianType = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType;

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // Get geometry
        auto geo = eg.geometry();

        // interpolate f as finite element function to compute the gradient
        CD_RHS_LocalAdapter<T,EG> f_adapter(param,eg);
        std::vector<RF> f_up, f_down, f_mid;
        param.setTime(time);
        lfsu.finiteElement().localInterpolation().interpolate(f_adapter,f_down);
        param.setTime(time+0.5*dt);
        lfsu.finiteElement().localInterpolation().interpolate(f_adapter,f_mid);
        param.setTime(time+dt);
        lfsu.finiteElement().localInterpolation().interpolate(f_adapter,f_up);

        // Initialize vectors outside for loop
        std::vector<JacobianType> js(lfsu.size());
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
        Dune::FieldVector<RF,dim> gradu(0.0);
        Dune::FieldVector<RF,dim> gradf_down(0.0);
        Dune::FieldVector<RF,dim> gradf_mid(0.0);
        Dune::FieldVector<RF,dim> gradf_up(0.0);
        Dune::FieldVector<RF,dim> gradf_average(0.0);

        // Transformation matrix
        typename EG::Geometry::JacobianInverseTransposed jac;

        // loop over quadrature points
        RF sum(0.0);
        RF sum_grad(0.0);
        RF fsum_grad_up(0.0);
        RF fsum_grad_mid(0.0);
        RF fsum_grad_down(0.0);
        auto intorder = 2*lfsu.finiteElement().localBasis().order();
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(ip.position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // integrate jump
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            sum += u*u*factor;

            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            lfsu.finiteElement().localBasis().evaluateJacobian(ip.position(),js);

            // transform gradients of shape functions to real element
            jac = geo.jacobianInverseTransposed(ip.position());
            for (size_type i=0; i<lfsu.size(); i++)
              jac.mv(js[i][0],gradphi[i]);

            // compute gradient of u
            gradu = 0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              gradu.axpy(x(lfsu,i),gradphi[i]);

            // integrate jump of gradient
            sum_grad += (gradu*gradu)*factor;

            // compute gradients of f
            gradf_down = 0.0;
            for (size_type i=0; i<lfsu.size(); i++) gradf_down.axpy(f_down[i],gradphi[i]);
            gradf_mid = 0.0;
            for (size_type i=0; i<lfsu.size(); i++) gradf_mid.axpy(f_mid[i],gradphi[i]);
            gradf_up = 0.0;
            for (size_type i=0; i<lfsu.size(); i++) gradf_up.axpy(f_up[i],gradphi[i]);
            gradf_average = 0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              gradf_average.axpy((1.0/6.0)*f_down[i]+(2.0/3.0)*f_mid[i]+(1.0/6.0)*f_up[i],gradphi[i]);

            // integrate grad(f-f_average)
            gradf_down -= gradf_average;
            fsum_grad_down += (gradf_down*gradf_down)*factor;
            gradf_mid -= gradf_average;
            fsum_grad_mid += (gradf_mid*gradf_mid)*factor;
            gradf_up -= gradf_average;
            fsum_grad_up += (gradf_up*gradf_up)*factor;
          }

        // accumulate cell indicator
        auto h_T = diameter(geo);
        r.accumulate(lfsv,0,dt    * sum_grad);  // k_n*||grad(jump)||^2
        r.accumulate(lfsv,0,dt*dt * dt*((1.0/6.0)*fsum_grad_down+(2.0/3.0)*fsum_grad_mid+(1.0/6.0)*fsum_grad_up)); // k_n^2*||grad(f-time_average(f))||^2_s_t
      }

      // boundary integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // Define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;

        // Reference to the inside cell
        const auto& cell_inside = ig.inside();

        // Get geometry
        auto geo = ig.geometry();
        auto geo_inside = cell_inside.geometry();

        // Get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // evaluate boundary condition
        auto ref_el_in_inside = referenceElement(geo_in_inside);
        auto face_local = ref_el_in_inside.position(0,0);
        auto bctype = param.bctype(ig.intersection(),face_local);
        if (bctype != ConvectionDiffusionBoundaryConditions::Neumann)
          return;

        // loop over quadrature points and integrate
        RF sum_up(0.0);
        RF sum_down(0.0);
        const int intorder = 2*lfsu_s.finiteElement().localBasis().order();
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate flux boundary condition
            param.setTime(time);
            auto j_down = param.j(ig.intersection(),ip.position());
            param.setTime(time+0.5*dt);
            auto j_mid = param.j(ig.intersection(),ip.position());
            param.setTime(time+dt);
            auto j_up = param.j(ig.intersection(),ip.position());

            // integrate
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            sum_down += (j_down-j_mid)*(j_down-j_mid)*factor;
            sum_up += (j_up-j_mid)*(j_up-j_mid)*factor;
          }

        // accumulate indicator
        //DF h_T = diameter(ig.geometry());
        auto h_T = diameter(geo_inside);
        r_s.accumulate(lfsv_s,0,(h_T+dt*dt/h_T)*dt*0.5*(sum_down+sum_up));
      }

    private:
      T& param;  // two phase parameter class
      double time;
      double dt;

      template<class GEO>
      typename GEO::ctype diameter (const GEO& geo) const
      {
        using DF = typename GEO::ctype;
        DF hmax = -1.0E00;
        const int dim = GEO::coorddimension;
        for (int i=0; i<geo.corners(); i++)
          {
            Dune::FieldVector<DF,dim> xi = geo.corner(i);
            for (int j=i+1; j<geo.corners(); j++)
              {
                Dune::FieldVector<DF,dim> xj = geo.corner(j);
                xj -= xi;
                hmax = std::max(hmax,xj.two_norm());
              }
          }
        return hmax;
      }

    };

  } // end namespace PDELab
} // end namespace Dune
#endif // DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONFEM_HH
