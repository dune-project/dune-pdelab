// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_DIFFUSIONCCFV_HH
#define DUNE_PDELAB_DIFFUSIONCCFV_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/geometry/referenceelements.hh>

#include"defaultimp.hh"
#include"pattern.hh"
#include"flags.hh"
#include"diffusionparam.hh"


namespace Dune {
  namespace PDELab {

	// a local operator for solving the Laplace equation with Dirichlet boundary conditions
	//     - div (k(x) grad u) + a0*u = f in \Omega, 
    //                              u = g on \partial\Omega_D
    //            - k(x) grad u * \nu = j on \partial\Omega_N
	// with cell centered finite volumes on axiparallel cube grids
    // K : scalar permeability field
    // A0: Helmholtz Term
    // F : source term
    // B : boundary condition function
    // J : flux boundary condition
    // G : grid function for Dirichlet boundary conditions
    template<typename K, typename A0, typename F, typename B, typename J, typename G>
	class DiffusionCCFV : public NumericalJacobianApplySkeleton<DiffusionCCFV<K,A0,F,B,J,G> >,
                          public NumericalJacobianApplyBoundary<DiffusionCCFV<K,A0,F,B,J,G> >,
                          public NumericalJacobianApplyVolume<DiffusionCCFV<K,A0,F,B,J,G> >,
                          public NumericalJacobianSkeleton<DiffusionCCFV<K,A0,F,B,J,G> >,
                          public NumericalJacobianBoundary<DiffusionCCFV<K,A0,F,B,J,G> >,
                          public NumericalJacobianVolume<DiffusionCCFV<K,A0,F,B,J,G> >,
                          public FullSkeletonPattern, 
                          public FullVolumePattern,
                          public LocalOperatorDefaultFlags

	{
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

      DiffusionCCFV (const K& k_, const A0& a0_, const F& f_, const B& b_, const J& j_, const G& g_)
        : k(k_), a0(a0_), f(f_), b(b_), j(j_), g(g_)
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
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate Helmholtz term
        typename A0::Traits::RangeType a0value;
        a0.evaluate(eg.entity(),inside_local,a0value);

        r.accumulate(lfsu,0,a0value*x(lfsu,0)*eg.geometry().volume());

        // evaluate source term
        typename F::Traits::RangeType fvalue;
        f.evaluate(eg.entity(),inside_local,fvalue);

        r.accumulate(lfsu,0,-fvalue*eg.geometry().volume());
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

        // center in face's reference element
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);

        // face volume for integration
        RF face_volume = ig.geometry().integrationElement(face_local)
          *Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).volume();

        // cell centers in references elements
        const Dune::FieldVector<DF,IG::dimension>& 
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>& 
          outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);

        // evaluate diffusion coefficient
        typename K::Traits::RangeType k_inside, k_outside;
        k.evaluate(*(ig.inside()),inside_local,k_inside);
        k.evaluate(*(ig.outside()),outside_local,k_outside);
        typename K::Traits::RangeType k_avg = 2.0/(1.0/(k_inside+1E-30) + 1.0/(k_outside+1E-30));

        // cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension> 
          inside_global = ig.inside()->geometry().global(inside_local);
        Dune::FieldVector<DF,IG::dimension> 
          outside_global = ig.outside()->geometry().global(outside_local);

        // distance between the two cell centers
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();
        
        // contribution to residual on inside element, other residual is computed by symmetric call
        r_s.accumulate(lfsu_s,0,k_avg*(x_s(lfsu_s,0)-x_n(lfsu_n,0))*face_volume/distance);
        r_n.accumulate(lfsu_n,0,-k_avg*(x_s(lfsu_s,0)-x_n(lfsu_n,0))*face_volume/distance);
	  }

	  // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save som e geometry evaluations
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

        // center in face's reference element
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);

        // face volume for integration
        RF face_volume = ig.geometry().integrationElement(face_local)
          *Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).volume();
        
        // cell center in reference element
        const Dune::FieldVector<DF,IG::dimension>& 
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);

        // evaluate boundary condition type
        typename B::Traits::RangeType bctype;
        b.evaluate(ig,face_local,bctype);

        if (DiffusionBoundaryCondition::isDirichlet(bctype))
          {
            // Dirichlet boundary
            // distance between cell center and face center
            Dune::FieldVector<DF,IG::dimension> 
              inside_global = ig.inside()->geometry().global(inside_local);
            Dune::FieldVector<DF,IG::dimension> 
              outside_global = ig.geometry().global(face_local);
            inside_global -= outside_global;
            RF distance = inside_global.two_norm();
            
            // evaluate diffusion coefficient
            typename K::Traits::RangeType k_inside;
            k.evaluate(*(ig.inside()),inside_local,k_inside);
            
            // evaluate boundary condition function
            typename G::Traits::DomainType x = ig.geometryInInside().global(face_local);
            typename G::Traits::RangeType y;
            g.evaluate(*(ig.inside()),x,y);
            
            // contribution to residual on inside element
            r_s.accumulate(lfsu_s,0,k_inside*(x_s(lfsu_s,0)-y[0])*face_volume/distance);
          }
        else // if (DiffusionBoundaryCondition::isNeumann(bctype))
          {
            // Neumann boundary
            // evaluate flux boundary condition

            //evaluate boundary function
            typename J::Traits::DomainType x = ig.geometryInInside().global(face_local);
            typename J::Traits::RangeType jvalue;
            j.evaluate(*(ig.inside()),x,jvalue);

            // contribution to residual on inside element
            r_s.accumulate(lfsu_s,0,jvalue*face_volume);
          }
	  }

    private:
      const K& k;
      const A0& a0;
      const F& f;
      const B& b;
      const J& j;
      const G& g;
	};

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
