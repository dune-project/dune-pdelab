// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_P12DLAPLACE_HH
#define DUNE_PDELAB_P12DLAPLACE_HH

#include<dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>

#include"../common/geometrywrapper.hh"


namespace Dune {
  namespace PDELab {

	// a local operator for solving the Poisson equation
	// - \Delta u = f in \Omega, u = g on \partial\Omega
	// with P1 conforming finite elements on triangles
	class LaplaceDirichletP12D
	{
	  typedef typename B::Traits::JacobianType JacobianType;
	  typedef typename B::Traits::DomainFieldType DomainFieldType;
	  typedef typename B::Traits::RangeType RangeType;

	public:

	  // export assembler type
      enum { doAlphaVolume = true };
      enum { doAlphaSkeleton = false };
      enum { doAlphaBoundary = false };
      enum { doLambdaVolume = false };
      enum { doLambdaSkeleton = false };
      enum { doLambdaBoundary = false };

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
		  gradu.axpy(x[lfs.localIndex(i)],gradphi[i]);

		// integrate grad u * grad phi_i (0.5 is the area of the reference element)
		RF area = 0.5*eg.geometry().integrationElement(integrationpoint);
		for (int i=0; i<3; i++)
		  r[lfsv.localIndex(i)] += (gradu*gradphi[i])*area;
	  }
	};


    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune
