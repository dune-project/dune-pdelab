// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_DIFFUSION_HH
#define DUNE_PDELAB_DIFFUSION_HH

#include<vector>

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
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** a local operator for solving the diffusion equation
     *
     * \f{eqnarray*}
     * - \nabla\cdot\{K(x) \nabla u\} + a_0 u &=& f \mbox{ in } \Omega,          \ \
     *                                      u &=& g \mbox{ on } \partial\Omega_D \\
     *              -(K(x)\nabla u) \cdot \nu &=& j \mbox{ on } \partial\Omega_N \\
     * \f}
     * with conforming finite elements on all types of grids in any dimension
     * \tparam F grid function type giving f
     * \tparam B grid function type selecting boundary condition
     * \tparam J grid function type giving j
     */
    template<typename K, typename A0, typename F, typename B, typename J, int qorder=1>
	class Diffusion : public NumericalJacobianApplyVolume<Diffusion<K,A0,F,B,J,qorder> >,
                      public FullVolumePattern
      //,public NumericalJacobianVolume<Diffusion<K,A0,F,B,J,qorder> >
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

      Diffusion (const K& k_, const A0& a0_, const F& f_, const B& b_, const J& j_)
        : k(k_), a0(a0_), f(f_), b(b_), j(j_)
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
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        typename K::Traits::RangeType tensor;
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
        k.evaluate(eg.entity(),localcenter,tensor);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> js(lfsu.size());
            lfsu.localFiniteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradient to real element
            const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            for (int i=0; i<lfsu.size(); i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (int i=0; i<lfsu.size(); i++)
              gradu.axpy(x[i],gradphi[i]);

            // compute K * gradient of u
            Dune::FieldVector<RF,dim> Kgradu(0.0);
            tensor.umv(gradu,Kgradu);

            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.localFiniteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate u
            RF u=0.0;
            for (int i=0; i<lfsu.size(); i++)
              u += x[i]*phi[i];

            // evaluate Helmholtz term
            typename A0::Traits::RangeType y;
            a0.evaluate(eg.entity(),it->position(),y);

            // integrate (K grad u)*grad phi_i + a_0*u*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (int i=0; i<lfsu.size(); i++)
              r[i] += ( Kgradu*gradphi[i] + y*u*phi[i] )*factor;
          }
	  }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
                            LocalMatrix<R>& mat) const
      {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        typename K::Traits::RangeType tensor;
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
        k.evaluate(eg.entity(),localcenter,tensor);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> js(lfsu.size());
            lfsu.localFiniteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradient to real element
            const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            for (int i=0; i<lfsu.size(); i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // compute K * gradient of shape functions
            std::vector<Dune::FieldVector<RF,dim> > Kgradphi(lfsu.size());
            for (int i=0; i<lfsu.size(); i++)
              tensor.mv(gradphi[i],Kgradphi[i]);
            
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.localFiniteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate Helmholtz term
            typename A0::Traits::RangeType y;
            a0.evaluate(eg.entity(),it->position(),y);

            // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (int j=0; j<lfsu.size(); j++)
              for (int i=0; i<lfsu.size(); i++)
                mat(i,j) += ( Kgradphi[j]*gradphi[i] + y*phi[j]*phi[i] )*factor;
          }
      }

 	  // volume integral depending only on test functions
	  template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
		// domain and range field type
		typedef typename LFSV::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSV::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSV::Traits::LocalFiniteElementType::
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
            std::vector<RangeType> phi(lfsv.size());
            lfsv.localFiniteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate right hand side parameter function
            typename F::Traits::RangeType y;
            f.evaluate(eg.entity(),it->position(),y);

            // integrate f
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (int i=0; i<lfsv.size(); i++)
              r[i] -= y*phi[i]*factor;
          }
      }

      // boundary integral independen of ansatz functions
 	  template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
		// domain and range field type
		typedef typename LFSV::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSV::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSV::Traits::LocalFiniteElementType::
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
            std::vector<RangeType> phi(lfsv.size());
            lfsv.localFiniteElement().localBasis().evaluateFunction(local,phi);
            
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
      const K& k;
      const A0& a0;
      const F& f;
      const B& b;
      const J& j;
	};

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
