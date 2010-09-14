// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_POISSON_HH
#define DUNE_PDELAB_POISSON_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/idefault.hh>

#include"../common/geometrywrapper.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"pattern.hh"
#include"flags.hh"

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** a local operator for solving the Poisson equation
     *
     * \f{align*}{
     *           - \Delta u &=& f \mbox{ in } \Omega,          \\
     *                    u &=& g \mbox{ on } \partial\Omega_D \\
     *  -\nabla u \cdot \nu &=& j \mbox{ on } \partial\Omega_N \\
     * \f}
     * with conforming finite elements on all types of grids in any dimension
     * \tparam F grid function type giving f
     * \tparam B grid function type selecting boundary condition
     * \tparam J grid function type giving j
     */
    template<typename F, typename B, typename J, int qorder=1>
	class Poisson : public NumericalJacobianApplyVolume<Poisson<F,B,J,qorder> >,
                    public NumericalJacobianVolume<Poisson<F,B,J,qorder> >,
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

      Poisson (const F& f_, const B& b_, const J& j_)
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
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> jsu(lfsu.size());
            lfsu.localFiniteElement().localBasis().evaluateJacobian(it->position(),jsu);
            std::vector<JacobianType> jsv(lfsv.size());
            lfsv.localFiniteElement().localBasis().evaluateJacobian(it->position(),jsv);

            // transform gradients to real element
            const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphiu(lfsu.size());
            for (size_t i=0; i<lfsu.size(); i++)
              {
                gradphiu[i] = 0.0;
                jac.umv(jsu[i][0],gradphiu[i]);
              }
            std::vector<Dune::FieldVector<RF,dim> > gradphiv(lfsv.size());
            for (size_t i=0; i<lfsv.size(); i++)
              {
                gradphiv[i] = 0.0;
                jac.umv(jsv[i][0],gradphiv[i]);
              }

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (size_t i=0; i<lfsu.size(); i++)
              gradu.axpy(x[lfsu.localIndex(i)],gradphiu[i]);

            // integrate grad u * grad phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_t i=0; i<lfsv.size(); i++)
              r[lfsv.localIndex(i)] += (gradu*gradphiv[i])*factor;
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
            for (size_t i=0; i<lfsv.size(); i++)
              r[lfsv.localIndex(i)] -= y*phi[i]*factor;
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

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
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
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

            // evaluate test shape functions 
            std::vector<RangeType> phi(lfsv.size());
            lfsv.localFiniteElement().localBasis().evaluateFunction(local,phi);
            
            // evaluate flux boundary condition
            typename J::Traits::RangeType y;
            j.evaluate(*(ig.inside()),local,y);
            
            // integrate J
            RF factor = it->weight()*ig.geometry().integrationElement(it->position());
            for (size_t i=0; i<lfsv.size(); i++)
              r[lfsv.localIndex(i)] += y*phi[i]*factor;
          }
      }

    protected:
      const F& f;
      const B& b;
      const J& j;
	};

    //! \brief a local operator for solving the Poisson equation in
    //!        instationary problems
    /**
     * \f{align*}{
     *           - \Delta u &=& f \mbox{ in } \Omega,          \\
     *                    u &=& g \mbox{ on } \partial\Omega_D \\
     *  -\nabla u \cdot \nu &=& j \mbox{ on } \partial\Omega_N \\
     * \f}
     * with conforming finite elements on all types of grids in any dimension
     *
     * \tparam F Grid function type giving f
     * \tparam B Grid function type selecting boundary condition
     * \tparam J Grid function type giving j
     *
     * \note The grid functions need to support the member function setTime().
     */
    template<typename Time, typename F, typename B, typename J, int qorder=1>
    class InstationaryPoisson
      : public Poisson<F,B,J,qorder>,
        public InstationaryLocalOperatorDefaultMethods<Time>
    {
      typedef Poisson<F,B,J,qorder> Base;
      typedef InstationaryLocalOperatorDefaultMethods<Time> IDefault;

    protected:
      // need non-const references for setTime()
      F& f;
      B& b;
      J& j;

    public:
      //! construct InstationaryPoisson
      InstationaryPoisson(F& f_, B& b_, J& j_)
        : Base(f_, b_, j_), f(f_), b(b_), j(j_)
      {}

      //! set the time for subsequent evaluation on the parameter functions
      void setTime (Time t) {
        f.setTime(t);
        b.setTime(t);
        j.setTime(t);
        IDefault::setTime(t);
      }
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
