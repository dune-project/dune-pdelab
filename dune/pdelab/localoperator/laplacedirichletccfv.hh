// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LAPLACEDIRICHLETCCFV_HH
#define DUNE_PDELAB_LAPLACEDIRICHLETCCFV_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/common/referenceelements.hh>

#include"../common/geometrywrapper.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"../gridoperatorspace/gridoperatorspaceutilities.hh"
#include"pattern.hh"


namespace Dune {
  namespace PDELab {

	// a local operator for solving the Laplace equation with Dirichlet boundary conditions
	//     - \Delta u = 0 in \Omega, 
    //              u = g on \partial\Omega
	// with cell centered finite volumes on axiparallel cube grids
    // G : grid function for Dirichlet boundary conditions
    template<typename G>
	class LaplaceDirichletCCFV : public NumericalJacobianApplySkeleton<LaplaceDirichletCCFV<G> >,
                                 public NumericalJacobianApplyBoundary<LaplaceDirichletCCFV<G> >,
                                 public NumericalJacobianSkeleton<LaplaceDirichletCCFV<G> >,
                                 public NumericalJacobianBoundary<LaplaceDirichletCCFV<G> >,
                                 public FullSkeletonPattern, 
                                 public FullVolumePattern
	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

	  // residual assembly flags
      enum { doAlphaVolume    = false };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume   = false };
      enum { doLambdaSkeleton = false };
      enum { doLambdaBoundary = false };

      LaplaceDirichletCCFV (const G& g_) : g(g_) {}

	  // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_skeleton (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           R& r_s, R& r_n) const
	  {
 		// domain and range field type
 		typedef typename LFSU::Traits::LocalFiniteElementType::
 		  Traits::LocalBasisType::Traits::DomainFieldType DF;
 		typedef typename LFSU::Traits::LocalFiniteElementType::
 		  Traits::LocalBasisType::Traits::RangeFieldType RF;

        // center in face's reference element
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.intersectionGlobal().type()).position(0,0);

        // face volume for integration
        RF face_volume = ig.intersectionGlobal().integrationElement(face_local)
          *Dune::ReferenceElements<DF,IG::dimension>::general(ig.intersectionGlobal().type()).volume();

        // cell centers in references elements
        const Dune::FieldVector<DF,IG::dimension>& 
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>& 
          outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);

        // distance between the two cell centers
        Dune::FieldVector<DF,IG::dimension> 
          inside_global = ig.inside()->geometry().global(inside_local);
        Dune::FieldVector<DF,IG::dimension> 
          outside_global = ig.outside()->geometry().global(outside_local);
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();
        
        // contribution to residual on inside element, other residual is computed by symmetric call
        r_s[0] += (x_s[0]-x_n[0])*face_volume/distance;
        r_n[0] -= (x_s[0]-x_n[0])*face_volume/distance;
	  }

	  // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save som e geometry evaluations
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_boundary (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
	  {
 		// domain and range field type
 		typedef typename LFSU::Traits::LocalFiniteElementType::
 		  Traits::LocalBasisType::Traits::DomainFieldType DF;
 		typedef typename LFSU::Traits::LocalFiniteElementType::
 		  Traits::LocalBasisType::Traits::RangeFieldType RF;

        // center in face's reference element
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.intersectionGlobal().type()).position(0,0);

        // face volume for integration
        RF face_volume = ig.intersectionGlobal().integrationElement(face_local)
          *Dune::ReferenceElements<DF,IG::dimension>::general(ig.intersectionGlobal().type()).volume();
        
        // cell center in reference element
        const Dune::FieldVector<DF,IG::dimension>& 
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);

        // distance between cell center and face center
        Dune::FieldVector<DF,IG::dimension> 
          inside_global = ig.inside()->geometry().global(inside_local);
        Dune::FieldVector<DF,IG::dimension> 
          outside_global = ig.intersectionGlobal().global(face_local);
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();

        // evaluate boundary condition function
        typename G::Traits::DomainType x = ig.intersectionSelfLocal().global(face_local);
        typename G::Traits::RangeType y;
        g.evaluate(*(ig.inside()),x,y);

        // contribution to residual on inside element
        r_s[0] += (x_s[0]-y[0])*face_volume/distance;
	  }

    private:
      const G& g;
	};

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
