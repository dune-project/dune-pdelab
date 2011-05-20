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

#include <dune/localfunctions/common/interfaceswitch.hh>

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

      Poisson (const F& f_, const B& bctype_, const J& j_)
        : f(f_), bctype(bctype_), j(j_)
      {}

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
		// domain and range field type
        typedef FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType
          > FESwitch;
        typedef BasisInterfaceSwitch<
          typename FESwitch::Basis
          > BasisSwitch;
        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::RangeField RF;

        // dimensions
        static const int dimLocal = EG::Geometry::mydimension;
        static const int dimGlobal = EG::Geometry::coorddimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dimLocal>& rule =
          Dune::QuadratureRules<DF,dimLocal>::rule(gt,qorder);

        // loop over quadrature points
        for(typename Dune::QuadratureRule<DF,dimLocal>::const_iterator it =
              rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions
            // (we assume Galerkin method lfsu=lfsv)
            std::vector<Dune::FieldMatrix<RF,1,dimGlobal> >
              gradphiu(lfsu.size());
            BasisSwitch::gradient(FESwitch::basis(lfsu.finiteElement()),
                                  eg.geometry(), it->position(), gradphiu);
            std::vector<Dune::FieldMatrix<RF,1,dimGlobal> >
              gradphiv(lfsv.size());
            BasisSwitch::gradient(FESwitch::basis(lfsv.finiteElement()),
                                  eg.geometry(), it->position(), gradphiv);

            // compute gradient of u
            Dune::FieldVector<RF,dimGlobal> gradu(0.0);
            for (size_t i=0; i<lfsu.size(); i++)
              gradu.axpy(x(lfsu,i),gradphiu[i][0]);

            // integrate grad u * grad phi_i
            RF factor = r.weight() * it->weight() * eg.geometry().integrationElement(it->position());
            for (size_t i=0; i<lfsv.size(); i++)
              r.rawAccumulate(lfsv,i,(gradu*gradphiv[i][0])*factor);
          }
	  }

 	  // volume integral depending only on test functions
	  template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
		// domain and range field type
        typedef FiniteElementInterfaceSwitch<
          typename LFSV::Traits::FiniteElementType
          > FESwitch;
        typedef BasisInterfaceSwitch<
          typename FESwitch::Basis
          > BasisSwitch;
        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::DomainLocal DomainLocal;
        typedef typename BasisSwitch::RangeField RF;
        typedef typename BasisSwitch::Range Range;

        // dimensions
        static const int dimLocal = EG::Geometry::mydimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dimLocal>& rule =
          Dune::QuadratureRules<DF,dimLocal>::rule(gt,qorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dimLocal>::const_iterator it =
               rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate shape functions
            std::vector<Range> phi(lfsv.size());
            FESwitch::basis(lfsv.finiteElement()).
              evaluateFunction(it->position(),phi);

            // evaluate right hand side parameter function
            typename F::Traits::RangeType y(0.0);
            f.evaluate(eg.entity(),it->position(),y);

            // integrate f
            RF factor = - r.weight() * it->weight() * eg.geometry().integrationElement(it->position());
            for (size_t i=0; i<lfsv.size(); i++)
              r.rawAccumulate(lfsv,i,y*phi[i]*factor);
          }
      }

      // boundary integral independen of ansatz functions
 	  template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
		// domain and range field type
        typedef FiniteElementInterfaceSwitch<
          typename LFSV::Traits::FiniteElementType
          > FESwitch;
        typedef BasisInterfaceSwitch<
          typename FESwitch::Basis
          > BasisSwitch;
        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::DomainLocal DomainLocal;
        typedef typename BasisSwitch::RangeField RF;
        typedef typename BasisSwitch::Range Range;

        // dimensions
        static const int dimLocal = IG::Geometry::mydimension;

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dimLocal>& rule =
          Dune::QuadratureRules<DF,dimLocal>::rule(gtface,qorder);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dimLocal>::const_iterator it =
               rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate boundary condition type
            // skip rest if we are on Dirichlet boundary
            if( bctype.isDirichlet( ig,it->position() ) )
              continue;

            // position of quadrature point in local coordinates of element
            const DomainLocal& local =
              ig.geometryInInside().global(it->position());

            // evaluate test shape functions
            std::vector<Range> phi(lfsv.size());
            FESwitch::basis(lfsv.finiteElement()).evaluateFunction(local,phi);

            // evaluate flux boundary condition
            typename J::Traits::RangeType y(0.0);
            j.evaluate(*(ig.inside()),local,y);

            // integrate J
            RF factor = r.weight() * it->weight()*ig.geometry().integrationElement(it->position());
            for (size_t i=0; i<lfsv.size(); i++)
              r.rawAccumulate(lfsv,i,y*phi[i]*factor);
          }
      }

    protected:
      const F& f;
      const B& bctype;
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
      B& bctype;
      J& j;

    public:
      //! construct InstationaryPoisson
      InstationaryPoisson(F& f_, B& bctype_, J& j_)
        : Base(f_, bctype_, j_)
        , f(f_)
        , bctype(bctype_)
        , j(j_)
      {}

      //! set the time for subsequent evaluation on the parameter functions
      void setTime (Time t) {
        f.setTime(t);
        bctype.setTime(t);
        j.setTime(t);
        IDefault::setTime(t);
      }
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
