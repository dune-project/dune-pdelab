// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_POISSON_HH
#define DUNE_PDELAB_POISSON_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/grid/common/quadraturerules.hh>

#include"../common/geometrywrapper.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"../gridoperatorspace/gridoperatorspaceutilities.hh"
#include"pattern.hh"


namespace Dune {
  namespace PDELab {

	// a local operator for solving the Poisson equation
	//           - \Delta u = f in \Omega, 
    //                    u = g on \partial\Omega_D
    //  -\nabla u \cdot \nu = j on \partial\Omega_N
	// with conforming finite elements on all types of grids in any dimension
    // F : grid function type giving f
    // B : grid function type selecting boundary condition
    // J : grid function type giving j
    template<typename F, typename B, typename J, int qorder=1>
	class Poisson : public NumericalJacobianApplyVolume<Poisson<F,B,J,qorder> >,
                    public NumericalJacobianVolume<Poisson<F,B,J,qorder> >,
                    public FullVolumePattern
	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = false };

	  // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doAlphaSkeleton = false };
      enum { doAlphaBoundary = false };
      enum { doLambdaVolume = true };
      enum { doLambdaSkeleton = false };
      enum { doLambdaBoundary = true };

      Poisson (const F& f_ const B& b_ const J& j)
        : f(f_), b(b_), j(j_)
      {}

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;

		// define integration point (hard coded quadrature)
		Dune::FieldVector<DF,2> integrationpoint(1.0/3.0);

		// gradient of shape functions at integration point
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JT;
		std::vector<JT> gradients;
		lfsu.localFiniteElement().localBasis().evaluateJacobian(integrationpoint,gradients);

		// transformation of gradients to real element
		const Dune::FieldMatrix<DF,2,2> 
		  jac = eg.geometry().jacobianInverseTransposed(integrationpoint);
		Dune::FieldVector<RF,2> gradphi[3];
		for (int i=0; i<3; i++)
		  {
			gradphi[i] = 0.0;
			jac.umv(gradients[i][0],gradphi[i]);
		  }

		// compute gradient of solution at integration point
		Dune::FieldVector<RF,2> gradu(0.0);
		for (int i=0; i<3; i++)
		  gradu.axpy(x[lfsu.localIndex(i)],gradphi[i]);

		// integrate grad u * grad phi_i (0.5 is the area of the reference element)
		RF area = 0.5*eg.geometry().integrationElement(integrationpoint);
		for (int i=0; i<3; i++)
		  r[lfsv.localIndex(i)] += (gradu*gradphi[i])*area;
	  }

 	  // volume integral depending only on test functions
	  template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate shape functions 
            std::vector<RangeType> phi;
            lfsv.evaluateFunction(it->position(),phi);

            // evaluate right hand side parameter function
            typename F::Traits::RangeType y;
            f.evaluate(eg.entity(),it->position(),y);

            // integrate f
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (int i=0; i<lfsv.size(); i++)
              r[i] -= y*phi[i]*factor;
          }
      }

      // boundary integral
 	  template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
       {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = IG::dimension;
        const int dimw = IG::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gtface = ig.intersectionSelfLocal().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate boundary condition type
            typename B::Traits::RangeType bctype;
            b.evaluate(ig,it->position(),bctype);
 
            // skip rest if we are on Dirichlet boundary
            if (bctype>0) continue;

            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.intersectionSelfLocal().global(it->position());

            // evaluate test shape functions 
            std::vector<RangeType> phi;
            lfsv.evaluateFunction(local,phi);
            
            // evaluate flux boundary condition
            typename J::Traits::RangeType y;
            j.evaluate(*(ig.inside()),local,y);
            
            // integrate J
            RF factor = it->weight()*ig.intersectionGlobal().integrationElement(it->position());
            for (int i=0; i<lfsv.size(); i++)
              r[i] += y*phi[i]*factor;
          }
      }

    private:
      const F& f;
      const B& b;
      const J& j;
	};

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
