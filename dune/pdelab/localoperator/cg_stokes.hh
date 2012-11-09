// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CGSTOKES_HH
#define DUNE_PDELAB_CGSTOKES_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>

#include<dune/geometry/type.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/gridfunctionspace/localvector.hh>
#include<dune/pdelab/gridoperator/common/localmatrix.hh>

#include"defaultimp.hh"
#include"pattern.hh"
#include"idefault.hh"
#include"flags.hh"
#include"l2.hh"
#include"stokesparameter.hh"

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
      typedef StokesBoundaryCondition BC;

      static const bool navier = P::assemble_navier;
      static const bool full_tensor = P::assemble_full_tensor;

      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doLambdaVolume = true };
      enum { doLambdaBoundary = true };

      typedef P PhysicalParameters;

      TaylorHoodNavierStokes (const PhysicalParameters & p, std::size_t quadrature_order = 4)

        : _p(p)
        , _quadrature_order(quadrature_order)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // extract local function spaces
        typedef typename LFSU::template Child<0>::Type LFSU_V_PFS;
        const LFSU_V_PFS& lfsu_v_pfs = lfsu.template child<0>();

        typedef typename LFSU_V_PFS::template Child<0>::Type LFSU_V;
        const unsigned int vsize = lfsu_v_pfs.child(0).size();

        // domain and range field type
        typedef typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RT_V;
        typedef typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType_V;
        typedef typename LFSU_V::Traits::SizeType size_type;


        typedef typename LFSU::template Child<1>::Type LFSU_P;
        const LFSU_P& lfsu_p = lfsu.template child<1>();
        const unsigned int psize = lfsu_p.size();

        typedef typename LFSU_P::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU_P::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RT_P;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,_quadrature_order);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it = rule.begin(),
               endit = rule.end();
             it != endit;
             ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType_V> js(vsize);
            lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradient to real element
            const typename EG::Geometry::JacobianInverseTransposed jac =
              eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(vsize);
            for (size_t i=0; i<vsize; i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // evaluate basis functions
            std::vector<RT_P> psi(psize);
            lfsu_p.finiteElement().localBasis().evaluateFunction(it->position(),psi);

            // compute u (if Navier term enabled)
            Dune::FieldVector<RF,dim> vu(0.0);

            std::vector<RT_V> phi(vsize);
            if(navier)
              {
                lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateFunction(it->position(),phi);

                for(int d=0; d<dim; ++d)
                  {
                    const LFSU_V & lfsu_v = lfsu_v_pfs.child(d);
                    for (size_t i=0; i<lfsu_v.size(); i++)
                      vu[d] += x(lfsu_v,i) * phi[i];
                  }
              }

            // Compute velocity jacobian
            Dune::FieldMatrix<RF,dim,dim> jacu(0.0);
            for(int d=0; d<dim; ++d){
              const LFSU_V & lfsu_v = lfsu_v_pfs.child(d);
              for (size_t i=0; i<lfsu_v.size(); i++)
                jacu[d].axpy(x(lfsu_v,i),gradphi[i]);
            }

            // compute pressure
            RT_P func_p(0.0);
            for (size_t i=0; i<lfsu_p.size(); i++)
              func_p += x(lfsu_p,i) * psi[i];

            // Viscosity and density
            const RF mu = _p.mu(eg,it->position());
            const RF rho = _p.rho(eg,it->position());

            // geometric weight
            const RF factor = it->weight() * eg.geometry().integrationElement(it->position());

            for(int d=0; d<dim; ++d){

              const LFSU_V & lfsu_v = lfsu_v_pfs.child(d);

              //compute u * grad u_d
              const RF u_nabla_u = vu * jacu[d];

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

            // compute divergence of u
            RF divu(0.0);
            for (size_t i=0; i<vsize; i++)
              for(int d=0; d<dim; ++d)
                divu += x(lfsu_v_pfs.child(d),i) * gradphi[i][d];

            // integrate div u * psi_i
            for (size_t i=0; i<lfsu_p.size(); i++)
              {
                r.accumulate(lfsu_p,i, -1.0 * divu * psi[i] * factor);
              }

          }
      }


      // volume integral depending on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // dimensions
        const int dim = EG::Geometry::dimension;

        // extract local function spaces
        typedef typename LFSV::template Child<0>::Type LFSV_V_PFS;
        const LFSV_V_PFS& lfsv_v_pfs = lfsv.template child<0>();

        typedef typename LFSV_V_PFS::template Child<0>::Type LFSV_V;
        const unsigned int vsize = lfsv_v_pfs.child(0).size();

        // domain and range field type
        typedef typename LFSV_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RT_V;
        typedef typename LFSV_V::Traits::SizeType size_type;

        typedef typename LFSV::template Child<1>::Type LFSV_P;
        const LFSV_P& lfsv_p = lfsv.template child<1>();
        const unsigned int psize = lfsv_p.size();

        typedef typename LFSV_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV_P::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RT_P;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,_quadrature_order);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it = rule.begin(),
               endit = rule.end();
             it != endit;
             ++it)
          {
            std::vector<RT_V> phi(vsize);
            lfsv_v_pfs.child(0).finiteElement().localBasis().evaluateFunction(it->position(),phi);

            std::vector<RT_P> psi(psize);
            lfsv_p.finiteElement().localBasis().evaluateFunction(it->position(),psi);

            // forcing term
            const Dune::FieldVector<RF,dim> f1 = _p.f(eg,it->position());

            // geometric weight
            const RF factor = it->weight() * eg.geometry().integrationElement(it->position());

            for(int d=0; d<dim; ++d){

              const LFSV_V & lfsv_v = lfsv_v_pfs.child(d);

              for (size_t i=0; i<vsize; i++)
                {
                  // integrate f1 * phi_i
                  r.accumulate(lfsv_v,i, -f1[d]*phi[i] * factor);
                }

            }

            const RF g2 = _p.g2(eg,it->position());

            // integrate div u * psi_i
            for (size_t i=0; i<lfsv_p.size(); i++)
              {
                r.accumulate(lfsv_p,i, g2 * psi[i] * factor);
              }

          }
      }


      // residual of boundary term
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
        // dimensions
        static const unsigned int dim = IG::Geometry::dimension;
        static const unsigned int dimw = IG::Geometry::dimensionworld;

        // extract local velocity function spaces
        typedef typename LFSV::template Child<0>::Type LFSV_V_PFS;
        const LFSV_V_PFS& lfsv_v_pfs = lfsv.template child<0>();

        typedef typename LFSV_V_PFS::template Child<0>::Type LFSV_V;
        const unsigned int vsize = lfsv_v_pfs.child(0).size();

        typedef typename LFSV_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType_V;
        typedef typename LFSV_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RT_V;

        // the range field type (equal for velocity and pressure)
        typedef typename LFSV_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // the size type (equal for velocity and pressure)
        typedef typename LFSV_V::Traits::SizeType size_type;

        // the domain field type (equal for velocity and pressure)
        typedef typename LFSV_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,_quadrature_order);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it = rule.begin(),
               endit = rule.end();
             it != endit;
             ++it)
          {
            // evaluate boundary condition type
            typename BC::Type bctype = _p.bctype(ig,it->position());

            // skip rest if we are on Dirichlet boundary
            if (bctype != BC::StressNeumann)
              continue;

            // position of quadrature point in local coordinates of element
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

            // evaluate basis functions
            std::vector<RT_V> phi(vsize);
            lfsv_v_pfs.child(0).finiteElement().localBasis().evaluateFunction(local,phi);

            const RF factor = it->weight() * ig.geometry().integrationElement(it->position());
            const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(it->position());

            // evaluate flux boundary condition
            const Dune::FieldVector<DF,dimw> neumann_stress = _p.j(ig,it->position(),normal);

            for(unsigned int d=0; d<dim; ++d)
              {

                const LFSV_V & lfsv_v = lfsv_v_pfs.child(d);

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
        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;


        // extract local function spaces
        typedef typename LFSU::template Child<0>::Type LFSU_V_PFS;
        const LFSU_V_PFS& lfsu_v_pfs = lfsu.template child<0>();
        const unsigned int vsize = lfsu_v_pfs.child(0).size();

        typedef typename LFSU_V_PFS::template Child<0>::Type LFSU_V;

        // domain and range field type
        typedef typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF_V;
        typedef typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RT_V;
        typedef typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType_V;
        typedef typename LFSU_V::Traits::SizeType size_type;


        typedef typename LFSU::template Child<1>::Type LFSU_P;
        const LFSU_P& lfsu_p = lfsu.template child<1>();
        const unsigned int psize = lfsu_p.size();

        typedef typename LFSU_P::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        typedef typename LFSU_P::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RT_P;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,_quadrature_order);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it = rule.begin(),
               endit = rule.end();
             it != endit;
             ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType_V> js(vsize);
            lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradient to real element
            const typename EG::Geometry::JacobianInverseTransposed jac =
              eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(vsize);
            for (size_t i=0; i<vsize; i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // evaluate basis functions
            std::vector<RT_P> psi(psize);
            lfsu_p.finiteElement().localBasis().evaluateFunction(it->position(),psi);

            // compute u (if Navier term enabled)
            std::vector<RT_V> phi(vsize);
            Dune::FieldVector<RF,dim> vu(0.0);
            if(navier){
              lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateFunction(it->position(),phi);

              for(int d = 0; d < dim; ++d){
                const LFSU_V & lfsv_v = lfsu_v_pfs.child(d);
                for(size_t l = 0; l < vsize; ++l)
                  vu[d] += x(lfsv_v,l) * phi[l];
              }
            }

            // Viscosity and density
            const RF mu = _p.mu(eg,it->position());
            const RF rho = _p.rho(eg,it->position());

            const RF factor = it->weight() * eg.geometry().integrationElement(it->position());

            for(int d=0; d<dim; ++d){

              const LFSU_V & lfsv_v = lfsu_v_pfs.child(d);
              const LFSU_V & lfsu_v = lfsv_v;

              // Derivatives of d-th velocity component
              Dune::FieldVector<RF,dim> gradu_d(0.0);
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
                      const LFSU_V & lfsu_v = lfsu_v_pfs.child(dd);
                      mat.accumulate(lfsv_v,i,lfsu_v,j, mu * (gradphi[j][d] * gradphi[i][dd]) * factor);
                    }

                }

                // integrate grad_d phi_v_d * p_u (pressure force)
                for (size_t j=0; j<lfsu_p.size(); j++)
                  mat.accumulate(lfsv_v,i,lfsu_p,j, - (gradphi[i][d] * psi[j]) * factor);

                if(navier){
                  for(int k =0; k < dim; ++k){
                    const LFSU_V & lfsu_v = lfsu_v_pfs.child(k);

                    const RF pre_factor = factor * rho * gradu_d[k] * phi[i];

                    for(size_t j=0; j< lfsu_v.size(); ++j)
                      mat.accumulate(lfsv_v,i,lfsu_v,j, pre_factor * phi[j]);
                  } // k
                }

                if(navier){
                  const RF pre_factor = factor * rho *  phi[i];
                  for(size_t j=0; j< lfsu_v.size(); ++j)
                    mat.accumulate(lfsv_v,i,lfsu_v,j,  pre_factor * (vu * gradphi[j]));
                }

              }

              for (size_t i=0; i<lfsu_p.size(); i++){
                for (size_t j=0; j<lfsu_v.size(); j++)
                  mat.accumulate(lfsu_p,i,lfsu_v,j, - (gradphi[j][d] * psi[i]) * factor);
              }

            } // d

          } // it
      }

    private:
      const P& _p;
      const std::size_t _quadrature_order;
    };


    /** a local operator for the mass term corresponding to the
     * instationary local operator TaylorHoodNavierStokes(Jacobian)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     */
    template< typename P >
    class NavierStokesMass
      : public NumericalJacobianApplyVolume<NavierStokesMass<P> >
      , public FullVolumePattern
      , public LocalOperatorDefaultFlags
      , public InstationaryLocalOperatorDefaultMethods<double>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      NavierStokesMass (const P & p_, int intorder_=4)
        : p(p_), intorder(intorder_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        typedef typename LFSV::template Child<0>::Type LFSV_PFS_V;
        const LFSV_PFS_V& lfsv_pfs_v = lfsv.template child<0>();

        for(unsigned int i=0; i<LFSV_PFS_V::CHILDREN; ++i)
          {
            scalar_alpha_volume(eg,lfsv_pfs_v.child(i),x,lfsv_pfs_v.child(i),r);
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        typedef typename LFSV::template Child<0>::Type LFSV_PFS_V;
        const LFSV_PFS_V& lfsv_pfs_v = lfsv.template child<0>();

        for(unsigned int i=0; i<LFSV_PFS_V::CHILDREN; ++i)
          {
            scalar_jacobian_volume(eg,lfsv_pfs_v.child(i),x,lfsv_pfs_v.child(i),mat);
          }
      }

    private:
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void scalar_alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            R& r) const
      {

        // Switches between local and global interface
        typedef FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType
          > FESwitch;
        typedef BasisInterfaceSwitch<
          typename FESwitch::Basis
          > BasisSwitch;

        // domain and range field type
        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::RangeField RF;
        typedef typename BasisSwitch::Range RangeType;

        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            FESwitch::basis(lfsu.finiteElement()).evaluateFunction(it->position(),phi);

            RF rho = p.rho(eg,it->position());
            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // u*phi_i
            RF factor = it->weight() * rho * eg.geometry().integrationElement(it->position());

            for (size_type i=0; i<lfsu.size(); i++)
              r.accumulate(lfsv,i, u*phi[i]*factor);
          }
      }

      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void scalar_jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {

        // Switches between local and global interface
        typedef FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType
          > FESwitch;
        typedef BasisInterfaceSwitch<
          typename FESwitch::Basis
          > BasisSwitch;

        // domain and range field type
        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::RangeField RF;
        typedef typename BasisSwitch::Range RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            FESwitch::basis(lfsu.finiteElement()).evaluateFunction(it->position(),phi);

            // integrate phi_j*phi_i
            RF rho = p.rho(eg,it->position());
            RF factor = it->weight() * rho * eg.geometry().integrationElement(it->position());
            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat.accumulate(lfsv,i,lfsu,j, phi[j]*phi[i]*factor);
          }
      }

      const P & p;
      int intorder;
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
