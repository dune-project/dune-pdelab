// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_CGSTOKES_HH
#define DUNE_PDELAB_CGSTOKES_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>

#include"../common/geometrywrapper.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"pattern.hh"
#include"idefault.hh"
#include"flags.hh"
#include"l2.hh"

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
        u \cdot \nabla u \cdot v - \Delta u + \nabla p &=& 0 \mbox{ in } \Omega, \\
        \nabla \cdot u &=& 0 \mbox{ in } \Omega \\
        u &=& g \mbox{ on } \partial\Omega_D \\
        -\nu (\nabla u + p I ) \nu &=& j \mbox{ on } \partial\Omega_N \\
        \f}
      
        As indicated in the equation above, this implementation
        utilizes only scalar Neumann conditions.

      \tparam B Grid function type selecting boundary condition

      \tparam J Scalar grid function type giving j 

      \tparam P A suitable parameter class with the interface of
      TaylorHoodNavierStokesParameters

      \tparam navier May be set to false, to avoid assembling of
      navier term in case rho=0.

      \tparam q Quadrature order.
     */

    template<typename B, typename J, typename P, bool navier = true, int qorder=3>
	class TaylorHoodNavierStokes :
      public NumericalJacobianApplyVolume<TaylorHoodNavierStokes<B,J,P,navier,qorder> >,
      public NumericalJacobianVolume<TaylorHoodNavierStokes<B,J,P,navier,qorder> >,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<double>
	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };

	  // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doLambdaBoundary = true };

      typedef B BoundaryFunction;
      typedef J FluxFunction;
      typedef P PhysicalParameters;

      TaylorHoodNavierStokes (const BoundaryFunction& b_, 
                              const FluxFunction& j_,
                              const PhysicalParameters & p_
                              )

        : NumericalJacobianVolume< TaylorHoodNavierStokes<B,J,P,navier,qorder> >(1e-7), 
          b(b_), j(j_), p(p_)
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
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType_V> js(vsize);
            lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradient to real element
            const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
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
            if(navier){
              lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateFunction(it->position(),phi);

              for(int d=0; d<dim; ++d){
                const LFSU_V & lfsu_v = lfsu_v_pfs.child(d);
                for (size_t i=0; i<lfsu_v.size(); i++)
                  vu[d] += x(lfsu_v,i) * phi[i];
              }
            }

            for(int d=0; d<dim; ++d){

              const LFSU_V & lfsu_v = lfsu_v_pfs.child(d);

              // compute gradient of u
              Dune::FieldVector<RF,dim> gradu(0.0);
              for (size_t i=0; i<lfsu_v.size(); i++)
                gradu.axpy(x(lfsu_v,i),gradphi[i]);

              // compute pressure
              RT_P func_p(0.0);
              for (size_t i=0; i<lfsu_p.size(); i++)
                func_p += x(lfsu_p,i) * psi[i];

              //compute u * grad u_d
              const RF u_nabla_u = vu * gradu;

              // geometric weight 
              RF factor = it->weight() * eg.geometry().integrationElement(it->position());

              for (size_t i=0; i<vsize; i++){

                // integrate grad u * grad phi_i
                r.accumulate(lfsu_v,i, p.mu() * (gradu * gradphi[i]) * factor);

                // integrate div phi_i * p
                r.accumulate(lfsu_v,i,- (func_p * gradphi[i][d]) * factor);

                // integrate u * grad u * phi_i
                if(navier)
                  r.accumulate(lfsu_v,i, p.rho() * u_nabla_u * phi[i] * factor);
              }

            }

            // compute divergence of u
            RF divu(0.0);
            for (size_t i=0; i<vsize; i++)
              for(int d=0; d<dim; ++d)
                divu += x(lfsu_v_pfs.child(d),i) * gradphi[i][d];

            // integrate div u * psi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_t i=0; i<lfsu_p.size(); i++){
              r.accumulate(lfsu_p,i, - (divu * psi[i]) * factor);
            }

          }
	  }

      // jacobian of boundary term
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

        // extract velocity boundary function (equal for all
        // dimensions)
        typedef typename B::template Child<0>::Type B_V_PFS;
        const B_V_PFS & b_v_pfs = b.template child<0>();
        typedef typename B_V_PFS::template Child<0>::Type B_V;
        const B_V & boundary_function = b_v_pfs.template child<0>();

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate boundary condition type
            typename B_V::Traits::RangeType bctype;
            boundary_function.evaluate(ig,it->position(),bctype);
 
            // skip rest if we are on Dirichlet boundary
            if (bctype>0) continue;

            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

            // evaluate basis functions
            std::vector<RT_V> phi(vsize);
            lfsv_v_pfs.child(0).finiteElement().localBasis().evaluateFunction(local,phi);

            // evaluate flux boundary condition. the scalar flux is
            // assumed to be in normal direction
            typename J::Traits::RangeType neumann_flux(0);
            j.evaluate(*(ig.inside()), local, neumann_flux);
            
            const RF factor = it->weight() * ig.geometry().integrationElement(it->position());
            const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(it->position());

            for(unsigned int d=0; d<dim; ++d){

              const LFSV_V & lfsv_v = lfsv_v_pfs.child(d);

              for (size_t i=0; i<vsize; i++){
                r.accumulate(lfsv_v,i, p.mu() * neumann_flux * normal[d] * phi[i] * factor);
              }

            }
          }
      }

    protected:
      const B& b;
      const J& j;
      const P& p;
    };

    /**
       \brief A local operator for the Navier-Stokes equations with
       direct assembling of jacobian.

       This class is derived from TaylorHoodNavierStokes and provides
       the same interface and functionality.
     */

    template<typename B, typename J, typename P, bool navier, int qorder=2>
	class TaylorHoodNavierStokesJacobian :
      public JacobianBasedAlphaVolume< TaylorHoodNavierStokesJacobian<B,J,P,navier,qorder> >,
      public TaylorHoodNavierStokes<B,J,P,navier,qorder>
	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };

	  // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doLambdaBoundary = true };

      typedef B BoundaryFunction;
      typedef J FluxFunction;
      typedef P PhysicalParameters;

      typedef TaylorHoodNavierStokes<B,J,P,navier,qorder> Base;

      TaylorHoodNavierStokesJacobian (const BoundaryFunction& b_, 
                                      const FluxFunction& j_, 
                                      const PhysicalParameters & p_)

        : Base(b_,j_,p_), b(b_), j(j_), p(p_)
      {}


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
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType_V> js(vsize);
            lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradient to real element
            const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
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

            for(int d=0; d<dim; ++d){

              const LFSU_V & lfsv_v = lfsu_v_pfs.child(d);
              const LFSU_V & lfsu_v = lfsv_v;

              // Derivatives of d-th velocity component
              Dune::FieldVector<RF,dim> gradu_d(0.0);
              if(navier)
                for(size_t l =0; l < vsize; ++l)
                  gradu_d.axpy(x(lfsv_v,l), gradphi[l]);

              RF factor = it->weight() * eg.geometry().integrationElement(it->position());
              for (size_t i=0; i<lfsv_v.size(); i++){

                // integrate grad phi_u_i * grad phi_v_i (viscous force)
                for (size_t j=0; j<lfsv_v.size(); j++){
                  mat.accumulate(lfsv_v,i,lfsu_v,j, p.mu() * (gradphi[i] * gradphi[j]) * factor);
                }

                // integrate grad_d phi_v_d * p_u (pressure force)
                for (size_t j=0; j<lfsu_p.size(); j++)
                  mat.accumulate(lfsv_v,i,lfsu_p,j, - (gradphi[i][d] * psi[j]) * factor);

                if(navier){
                  for(int k =0; k < dim; ++k){
                    const LFSU_V & lfsu_v = lfsu_v_pfs.child(k);

                    const RF pre_factor = factor * p.rho() * gradu_d[k] * phi[i];

                    for(size_t j=0; j< lfsu_v.size(); ++j)
                      mat.accumulate(lfsv_v,i,lfsu_v,j, pre_factor * phi[j]);
                  } // k
                }

                if(navier){
                  const RF pre_factor = factor * p.rho() *  phi[i];
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

	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
        Base::alpha_volume(eg,lfsu,x,lfsv,r);
      }

      template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
        Base::lambda_boundary(ig,lfsv,r);
      }

    protected:
      const B& b;
      const J& j;
      const P& p;
    };

    //! Interface for the parameter class required by the classes
    //! TaylorHoodNavierStokes and TaylorHoodNavierStokesJacobian.
    template <class CT>
    class TaylorHoodNavierStokesParameters{
    public:
      CT rho() const{ return CT(1.0); }
      CT mu()  const{ return CT(1.0); }
    };

    //! Interface for the parameter class required by the classes
    //! TaylorHoodNavierStokes and TaylorHoodNavierStokesJacobian.
    template <class CT>
    class InstationaryTaylorHoodNavierStokesParameters{
    public:
      CT rho() const{ return CT(1.0); }
      CT mu()  const{ return CT(1.0); }
      CT tau()  const{ return CT(0.1); }
      CT time()  const{ return CT(1.0); }
    };

    /** a local operator for the mass term corresponding to the
     * instationary local operator TaylorHoodNavierStokes(Jacobian)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     */
	class NavierStokesMass : public NumericalJacobianApplyVolume<NavierStokesMass>,
                             public FullVolumePattern,
                             public LocalOperatorDefaultFlags,
                             public InstationaryLocalOperatorDefaultMethods<double>
	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };

	  // residual assembly flags
      enum { doAlphaVolume = true };

      NavierStokesMass (int intorder_=4)
        : intorder(intorder_),
          scalar_operator(intorder_)
      {}

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
        scalar_operator.alpha_volume(eg,lfsu.template child<0>(),x,lfsu.template child<0>(),r);
	  }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
	  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
                            M& mat) const
      {
        scalar_operator.jacobian_volume(eg,lfsu.template child<0>(),x,lfsu.template child<0>(),mat);
      }

    private:
      int intorder;
      PowerL2 scalar_operator;
	};

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
