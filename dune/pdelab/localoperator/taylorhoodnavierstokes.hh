// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_TAYLORHOODNAVIERSTOKES_HH
#define DUNE_PDELAB_LOCALOPERATOR_TAYLORHOODNAVIERSTOKES_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

#include<dune/geometry/type.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/gridfunctionspace/localvector.hh>

#include"defaultimp.hh"
#include"pattern.hh"
#include"idefault.hh"
#include"flags.hh"
#include"l2.hh"
#include"stokesparameter.hh"
#include"navierstokesmass.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** \brief A local operator for the Navier-Stokes equations.

        This class implements a local operator for conforming finite
        element discretizations of the Navier-Stokes equations with
        TaylorHood basis.

        \f{align*}{
        u \cdot \nabla u \cdot v - \nabla \cdot ( \nabla u + (\nabla u)^T + p I) &=& 0 \mbox{ in } \Omega, \\
        \nabla \cdot u &=& 0 \mbox{ in } \Omega \\
        u &=& g \mbox{ on } \partial\Omega_D \\
        -\nu (\nabla u + p I ) \nu &=& j \mbox{ on } \partial\Omega_N \\
        \f}

        As indicated in the equation above, this implementation
        utilizes only scalar Neumann conditions.

      \tparam P A suitable parameter class with the interface of
      TaylorHoodNavierStokesDefaultParameters

        \tparam navier May be set to false, to avoid assembling of
        navier term in case rho=0.

        \tparam q Quadrature order.
    */

    template<typename P>
    class TaylorHoodNavierStokes :
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename P::Traits::RangeField>
    {
    public:
      //! Boundary condition indicator type
      using BC = StokesBoundaryCondition;

      using InstatBase = InstationaryLocalOperatorDefaultMethods<typename P::Traits::RangeField>;
      using Real = typename InstatBase::RealType;

      static const bool navier = P::assemble_navier;
      static const bool full_tensor = P::assemble_full_tensor;

      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doLambdaVolume = true };
      enum { doLambdaBoundary = true };

      using PhysicalParameters = P;

      TaylorHoodNavierStokes (PhysicalParameters& p, int superintegration_order_ = 0)

        : _p(p)
        , superintegration_order(superintegration_order_)
      {}

      // set time in parameter class
      void setTime(Real t)
      {
        InstatBase::setTime(t);
        _p.setTime(t);
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // define types
        using namespace Indices;
        using LFSU_V_PFS = TypeTree::Child<LFSU,_0>;
        using LFSU_V = TypeTree::Child<LFSU_V_PFS,_0>;
        using LFSU_P = TypeTree::Child<LFSU,_1>;
        using RF = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using RT_V = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using JacobianType_V = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType;
        using RT_P = typename LFSU_P::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;

        // extract local function spaces
        const auto& lfsu_v_pfs = child(lfsu,_0);
        const unsigned int vsize = lfsu_v_pfs.child(0).size();
        const auto& lfsu_p = child(lfsu,_1);
        const unsigned int psize = lfsu_p.size();

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // get geometry
        auto geo = eg.geometry();

        // determine quadrature order
        const int v_order = lfsu_v_pfs.child(0).finiteElement().localBasis().order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int jac_order = geo.type().isSimplex() ? 0 : 1;
        const int qorder = 3*v_order - 1 + jac_order + det_jac_order + superintegration_order;

        // Initialize vectors outside for loop
        typename EG::Geometry::JacobianInverseTransposed jac;
        std::vector<Dune::FieldVector<RF,dim> > gradphi(vsize);
        std::vector<RT_P> psi(psize);
        Dune::FieldVector<RF,dim> vu(0.0);
        std::vector<RT_V> phi(vsize);
        Dune::FieldMatrix<RF,dim,dim> jacu(0.0);

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType_V> js(vsize);
            lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateJacobian(ip.position(),js);

            // transform gradient to real element
            jac = geo.jacobianInverseTransposed(ip.position());
            for (size_t i=0; i<vsize; i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // evaluate basis functions
            lfsu_p.finiteElement().localBasis().evaluateFunction(ip.position(),psi);

            // compute u (if Navier term enabled)
            if(navier)
              {
                lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateFunction(ip.position(),phi);

                for(int d=0; d<dim; ++d)
                  {
                    vu[d] = 0.0;
                    const auto& lfsu_v = lfsu_v_pfs.child(d);
                    for (size_t i=0; i<lfsu_v.size(); i++)
                      vu[d] += x(lfsu_v,i) * phi[i];
                  }
              }

            // Compute velocity jacobian
            for(int d=0; d<dim; ++d){
              jacu[d] = 0.0;
              const auto& lfsu_v = lfsu_v_pfs.child(d);
              for (size_t i=0; i<lfsu_v.size(); i++)
                jacu[d].axpy(x(lfsu_v,i),gradphi[i]);
            }

            // compute pressure
            RT_P func_p(0.0);
            for (size_t i=0; i<lfsu_p.size(); i++)
              func_p += x(lfsu_p,i) * psi[i];

            // Viscosity and density
            const auto mu = _p.mu(eg,ip.position());
            const auto rho = _p.rho(eg,ip.position());

            // geometric weight
            const auto factor = ip.weight() * geo.integrationElement(ip.position());

            for(int d=0; d<dim; ++d){

              const auto& lfsu_v = lfsu_v_pfs.child(d);

              //compute u * grad u_d
              const auto u_nabla_u = vu * jacu[d];

              for (size_t i=0; i<vsize; i++){

                // integrate grad u * grad phi_i
                r.accumulate(lfsu_v,i, mu * (jacu[d] * gradphi[i]) * factor);

                // integrate (grad u)^T * grad phi_i
                if (full_tensor)
                  for(int dd=0; dd<dim; ++dd)
                    r.accumulate(lfsu_v,i, mu * (jacu[dd][d] * gradphi[i][dd]) * factor);

                // integrate div phi_i * p
                r.accumulate(lfsu_v,i,- (func_p * gradphi[i][d]) * factor);

                // integrate u * grad u * phi_i
                if(navier)
                  r.accumulate(lfsu_v,i, rho * u_nabla_u * phi[i] * factor);
              }

            }

            // integrate div u * psi_i
            for (size_t i=0; i<psize; i++)
                for(int d=0; d<dim; ++d)
                    // divergence of u is the trace of the velocity jacobian
                    r.accumulate(lfsu_p,i, -1.0 * jacu[d][d] * psi[i] * factor);

          }
      }


      // volume integral depending on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // define types
        using namespace Indices;
        using LFSV_V_PFS = TypeTree::Child<LFSV,_0>;
        using LFSV_V = TypeTree::Child<LFSV_V_PFS,_0>;
        using LFSV_P = TypeTree::Child<LFSV,_1>;
        using RT_V = typename LFSV_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using RT_P = typename LFSV_P::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;

        // extract local function spaces
        const auto& lfsv_v_pfs = child(lfsv,_0);
        const unsigned int vsize = lfsv_v_pfs.child(0).size();
        const auto& lfsv_p = child(lfsv,_1);
        const unsigned int psize = lfsv_p.size();

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // get cell
        const auto& cell = eg.entity();

        // get geometry
        auto geo = eg.geometry();

        // determine quadrature order
        const int v_order = lfsv_v_pfs.child(0).finiteElement().localBasis().order();
        const int det_jac_order = geo.type().isSimplex() ?  0 : (dim-1);
        const int qorder = 2*v_order + det_jac_order + superintegration_order;

        // Initialize vectors outside for loop
        std::vector<RT_V> phi(vsize);
        std::vector<RT_P> psi(psize);

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            lfsv_v_pfs.child(0).finiteElement().localBasis().evaluateFunction(ip.position(),phi);

            lfsv_p.finiteElement().localBasis().evaluateFunction(ip.position(),psi);

            // forcing term
            const auto f1 = _p.f(cell,ip.position());

            // geometric weight
            const auto factor = ip.weight() * geo.integrationElement(ip.position());

            for(int d=0; d<dim; ++d){

              const auto& lfsv_v = lfsv_v_pfs.child(d);

              for (size_t i=0; i<vsize; i++)
                {
                  // integrate f1 * phi_i
                  r.accumulate(lfsv_v,i, -f1[d]*phi[i] * factor);
                }

            }

            const auto g2 = _p.g2(eg,ip.position());

            // integrate div u * psi_i
            for (size_t i=0; i<psize; i++)
              {
                r.accumulate(lfsv_p,i, g2 * psi[i] * factor);
              }

          }
      }


      // residual of boundary term
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
        // define types
        using namespace Indices;
        using LFSV_V_PFS = TypeTree::Child<LFSV,_0>;
        using LFSV_V = TypeTree::Child<LFSV_V_PFS,_0>;
        using RT_V = typename LFSV_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;

        // extract local velocity function spaces
        const auto& lfsv_v_pfs = child(lfsv,_0);
        const unsigned int vsize = lfsv_v_pfs.child(0).size();

        // dimensions
        static const unsigned int dim = IG::dimension;

        // get geometry
        auto geo = ig.geometry();

        // Get geometry of intersection in local coordinates of inside cell
        auto geo_in_inside = ig.geometryInInside();

        // determine quadrature order
        const int v_order = lfsv_v_pfs.child(0).finiteElement().localBasis().order();
        const int det_jac_order = geo_in_inside.type().isSimplex() ? 0 : (dim-2);
        const int jac_order = geo_in_inside.type().isSimplex() ? 0 : 1;
        const int qorder = 2*v_order + det_jac_order + jac_order + superintegration_order;

        // Initialize vectors outside for loop
        std::vector<RT_V> phi(vsize);

        // loop over quadrature points and integrate normal flux
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            // evaluate boundary condition type
            auto bctype = _p.bctype(ig,ip.position());

            // skip rest if we are on Dirichlet boundary
            if (bctype != BC::StressNeumann)
              continue;

            // position of quadrature point in local coordinates of element
            auto local = geo_in_inside.global(ip.position());

            // evaluate basis functions
            lfsv_v_pfs.child(0).finiteElement().localBasis().evaluateFunction(local,phi);

            const auto factor = ip.weight() * geo.integrationElement(ip.position());
            const auto normal = ig.unitOuterNormal(ip.position());

            // evaluate flux boundary condition
            const auto neumann_stress = _p.j(ig,ip.position(),normal);

            for(unsigned int d=0; d<dim; ++d)
              {

                const auto& lfsv_v = lfsv_v_pfs.child(d);

                for (size_t i=0; i<vsize; i++)
                  {
                    r.accumulate(lfsv_v,i, neumann_stress[d] * phi[i] * factor);
                  }

              }
          }
      }


      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // define types
        using namespace Indices;
        using LFSU_V_PFS = TypeTree::Child<LFSU,_0>;
        using LFSU_V = TypeTree::Child<LFSU_V_PFS,_0>;
        using LFSU_P = TypeTree::Child<LFSU,_1>;
        using RF = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using RT_V = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using JacobianType_V = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType;
        using RT_P = typename LFSU_P::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;

        // extract local function spaces
        const auto& lfsu_v_pfs = child(lfsu,_0);
        const unsigned int vsize = lfsu_v_pfs.child(0).size();
        const auto& lfsu_p = child(lfsu,_1);
        const unsigned int psize = lfsu_p.size();

        // dimensions
        const int dim = EG::Geometry::mydimension;

        // get geometry
        auto geo = eg.geometry();

        // determine quadrature order
        const int v_order = lfsu_v_pfs.child(0).finiteElement().localBasis().order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int jac_order = geo.type().isSimplex() ? 0 : 1;
        const int qorder = 3*v_order - 1 + jac_order + det_jac_order + superintegration_order;

        // Initialize vectors outside for loop
        typename EG::Geometry::JacobianInverseTransposed jac;
        std::vector<JacobianType_V> js(vsize);
        std::vector<Dune::FieldVector<RF,dim> > gradphi(vsize);
        std::vector<RT_P> psi(psize);
        std::vector<RT_V> phi(vsize);
        Dune::FieldVector<RF,dim> vu(0.0);
        Dune::FieldVector<RF,dim> gradu_d(0.0);

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateJacobian(ip.position(),js);

            // transform gradient to real element
            jac = geo.jacobianInverseTransposed(ip.position());
            for (size_t i=0; i<vsize; i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // evaluate basis functions
            lfsu_p.finiteElement().localBasis().evaluateFunction(ip.position(),psi);

            // compute u (if Navier term enabled)
            if(navier){
              lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateFunction(ip.position(),phi);

              for(int d = 0; d < dim; ++d){
                vu[d] = 0.0;
                const auto& lfsv_v = lfsu_v_pfs.child(d);
                for(size_t l = 0; l < vsize; ++l)
                  vu[d] += x(lfsv_v,l) * phi[l];
              }
            }

            // Viscosity and density
            const auto mu = _p.mu(eg,ip.position());
            const auto rho = _p.rho(eg,ip.position());

            const auto factor = ip.weight() * geo.integrationElement(ip.position());

            for(int d=0; d<dim; ++d){

              const auto& lfsv_v = lfsu_v_pfs.child(d);
              const auto& lfsu_v = lfsv_v;

              // Derivatives of d-th velocity component
              gradu_d = 0.0;
              if(navier)
                for(size_t l =0; l < vsize; ++l)
                  gradu_d.axpy(x(lfsv_v,l), gradphi[l]);

              for (size_t i=0; i<lfsv_v.size(); i++){

                // integrate grad phi_u_i * grad phi_v_i (viscous force)
                for (size_t j=0; j<lfsv_v.size(); j++){
                  mat.accumulate(lfsv_v,i,lfsu_v,j, mu * (gradphi[i] * gradphi[j]) * factor);

                  // integrate (grad phi_u_i)^T * grad phi_v_i (viscous force)
                  if(full_tensor)
                    for(int dd=0; dd<dim; ++dd){
                      const auto& lfsu_v = lfsu_v_pfs.child(dd);
                      mat.accumulate(lfsv_v,i,lfsu_v,j, mu * (gradphi[j][d] * gradphi[i][dd]) * factor);
                    }

                }

                // integrate grad_d phi_v_d * p_u (pressure force)
                for (size_t j=0; j<psize; j++)
                  mat.accumulate(lfsv_v,i,lfsu_p,j, - (gradphi[i][d] * psi[j]) * factor);

                if(navier){
                  for(int k =0; k < dim; ++k){
                    const auto& lfsu_v = lfsu_v_pfs.child(k);

                    const auto pre_factor = factor * rho * gradu_d[k] * phi[i];

                    for(size_t j=0; j< lfsu_v.size(); ++j)
                      mat.accumulate(lfsv_v,i,lfsu_v,j, pre_factor * phi[j]);
                  } // k
                }

                if(navier){
                  const auto pre_factor = factor * rho *  phi[i];
                  for(size_t j=0; j< lfsu_v.size(); ++j)
                    mat.accumulate(lfsv_v,i,lfsu_v,j,  pre_factor * (vu * gradphi[j]));
                }

              }

              for (size_t i=0; i<psize; i++){
                for (size_t j=0; j<lfsu_v.size(); j++)
                  mat.accumulate(lfsu_p,i,lfsu_v,j, - (gradphi[j][d] * psi[i]) * factor);
              }
            } // d
          } // it
      }

    private:
      P& _p;
      const int superintegration_order;
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_TAYLORHOODNAVIERSTOKES_HH
