// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_DIFFUSIONMIXED_HH
#define DUNE_PDELAB_DIFFUSIONMIXED_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>

#include"../common/geometrywrapper.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"../gridoperatorspace/gridoperatorspaceutilities.hh"
#include"pattern.hh"
#include"flags.hh"

namespace Dune {
  namespace PDELab {

	// a local operator for solving the Poisson equation
	//     div sigma +a_0 u = f         in \Omega, 
	//                sigma = -K grad u in \Omega, 
    //                    u = g         on \partial\Omega_D
    //      sigma \cdot \nu = j         on \partial\Omega_N
	// with H(div) conforming (mixed) finite elements
    // K : diffusion tensor dependent on position
    // A0: Helmholtz term
    // F : grid function type giving f
    // B : grid function type selecting boundary condition
    // G : grid function type giving g
    template<typename K, typename A0, typename F, typename B, typename G>
	class DiffusionMixed : public NumericalJacobianApplyVolume<DiffusionMixed<K,A0,F,B,G> >,
                           public NumericalJacobianVolume<DiffusionMixed<K,A0,F,B,G> >,
                           public FullVolumePattern,
                           public LocalOperatorDefaultFlags
	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };

	  // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doLambdaVolume = true };
      enum { doLambdaBoundary = true };

      DiffusionMixed (const K& k_, const A0& a0_, const F& f_, const B& b_, const G& g_, int qorder_v_=2, int qorder_p_=1)
        : k(k_), a0(a0_), f(f_), b(b_), g(g_), qorder_v(qorder_v_), qorder_p(qorder_p_)
      {}

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
        // select the two components
        typedef typename LFSU::template Child<0>::Type VelocitySpace;
        const VelocitySpace& velocityspace = lfsu.template getChild<0>();
        typedef typename LFSU::template Child<1>::Type PressureSpace;
        const PressureSpace& pressurespace = lfsu.template getChild<1>();

		// domain and range field type
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType VelocityJacobianType;
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType VelocityRangeType;
		typedef typename PressureSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType PressureRangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // evaluate transformation which must be linear
        Dune::FieldVector<DF,dim> pos; pos = 0.0;
        Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(pos);
        jac.invert();
        RF det = eg.geometry().integrationElement(pos);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        typename K::Traits::RangeType tensor;
        Dune::GeometryType gt = eg.geometry().type();
        Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
        k.evaluate(eg.entity(),localcenter,tensor);
        tensor.invert(); // need iverse for mixed method

        // \sigma\cdot v term
        const Dune::QuadratureRule<DF,dim>& vrule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder_v);
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=vrule.begin(); it!=vrule.end(); ++it)
          {
            // evaluate shape functions at ip (this is a Galerkin method)
            std::vector<VelocityRangeType> vbasis(velocityspace.size());
            velocityspace.localFiniteElement().localBasis().evaluateFunction(it->position(),vbasis);

            // transform basis vectors
            std::vector<VelocityRangeType> vtransformedbasis(velocityspace.size());
            for (int i=0; i<velocityspace.size(); i++)
              {
                vtransformedbasis[i] = 0.0;
                jac.umtv(vbasis[i],vtransformedbasis[i]);
              }

            // compute sigma
            VelocityRangeType sigma; sigma = 0.0;
            for (int i=0; i<velocityspace.size(); i++)
              sigma.axpy(x[velocityspace.localIndex(i)],vtransformedbasis[i]);

            // K^{-1} * sigma
            VelocityRangeType Kinvsigma;
            tensor.mv(sigma,Kinvsigma);

            // integrate  (K^{-1}*sigma) * phi_i
            RF factor = it->weight() / det;
            for (int i=0; i<velocityspace.size(); i++)
              r[velocityspace.localIndex(i)] += (Kinvsigma*vtransformedbasis[i])*factor;
          }

        // u div v term, div sigma q term, a0*u term
        const Dune::QuadratureRule<DF,dim>& prule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder_p);
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=prule.begin(); it!=prule.end(); ++it)
          {
            // evaluate shape functions at ip (this is a Galerkin method)
            std::vector<VelocityJacobianType> vbasis(velocityspace.size());
            velocityspace.localFiniteElement().localBasis().evaluateJacobian(it->position(),vbasis);
            std::vector<PressureRangeType> pbasis(pressurespace.size());
            pressurespace.localFiniteElement().localBasis().evaluateFunction(it->position(),pbasis);

            // compute u
            PressureRangeType u; u = 0.0;
            for (int i=0; i<pressurespace.size(); i++)
              u.axpy(x[pressurespace.localIndex(i)],pbasis[i]);

            // evaluate Helmholtz term
            typename A0::Traits::RangeType a0value;
            a0.evaluate(eg.entity(),it->position(),a0value);

            // integrate a0 * u * q
            RF factor = it->weight();
            for (int i=0; i<pressurespace.size(); i++)
              r[pressurespace.localIndex(i)] -= a0value*u*pbasis[i]*factor;

            // compute divergence of velocity basis functions on reference element
            std::vector<RF> divergence(velocityspace.size(),0.0);
            for (int i=0; i<velocityspace.size(); i++)
              for (int j=0; j<dim; j++) 
                divergence[i] += vbasis[i][j][j];

            // integrate sigma * phi_i
            for (int i=0; i<velocityspace.size(); i++)
              r[velocityspace.localIndex(i)] -= u*divergence[i]*factor;

            // compute divergence of sigma
            RF divergencesigma = 0.0;
            for (int i=0; i<velocityspace.size(); i++)
              divergencesigma += x[velocityspace.localIndex(i)]*divergence[i];

            // integrate div sigma * q
            for (int i=0; i<pressurespace.size(); i++)
              r[pressurespace.localIndex(i)] -= divergencesigma*pbasis[i]*factor;
          }
	  }

 	  // volume integral depending only on test functions
	  template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<1>::Type PressureSpace;
        const PressureSpace& pressurespace = lfsv.template getChild<1>();

		// domain and range field type
		typedef typename PressureSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename PressureSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename PressureSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType PressureRangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder_p);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate shape functions 
            std::vector<PressureRangeType> pbasis(pressurespace.size());
            pressurespace.localFiniteElement().localBasis().evaluateFunction(it->position(),pbasis);

            // evaluate right hand side parameter function
            typename F::Traits::RangeType y;
            f.evaluate(eg.entity(),it->position(),y);

            // integrate f
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (int i=0; i<pressurespace.size(); i++)
              r[pressurespace.localIndex(i)] += y*pbasis[i]*factor;
          }
      }

      // boundary integral independen of ansatz functions
 	  template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<0>::Type VelocitySpace;
        const VelocitySpace& velocityspace = lfsv.template getChild<0>();

		// domain and range field type
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType VelocityRangeType;

        // dimensions
        const int dim = IG::dimension;
        const int dimw = IG::dimensionworld;

        // evaluate transformation which must be linear
        Dune::FieldVector<DF,dim> pos; pos = 0.0;
        Dune::FieldMatrix<DF,dimw,dim> jac = ig.inside()->geometry().jacobianInverseTransposed(pos);
        jac.invert();
        RF det = ig.inside()->geometry().integrationElement(pos);

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder_v);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate boundary condition type
            typename B::Traits::RangeType bctype;
            b.evaluate(ig,it->position(),bctype);
 
            // skip rest if we are on Neumann boundary
            if (DiffusionBoundaryCondition::isNeumann(bctype)) continue;

            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

            // evaluate test shape functions 
            std::vector<VelocityRangeType> vbasis(velocityspace.size());
            velocityspace.localFiniteElement().localBasis().evaluateFunction(local,vbasis);
            
            // transform basis vectors
            std::vector<VelocityRangeType> vtransformedbasis(velocityspace.size());
            for (int i=0; i<velocityspace.size(); i++)
              {
                vtransformedbasis[i] = 0.0;
                jac.umtv(vbasis[i],vtransformedbasis[i]);
              }

            // evaluate Dirichlet boundary condition
            typename G::Traits::RangeType y;
            g.evaluate(*(ig.inside()),local,y);
            
            // integrate g v*normal
            RF factor = it->weight()*ig.geometry().integrationElement(it->position())/det;
            for (int i=0; i<velocityspace.size(); i++)
              r[velocityspace.localIndex(i)] += y*(vtransformedbasis[i]*ig.unitOuterNormal(it->position()))*factor;
          }
      }

    private:
      const K& k;
      const A0& a0;
      const F& f;
      const B& b;
      const G& g;
      int qorder_v; 
      int qorder_p;
	};

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
