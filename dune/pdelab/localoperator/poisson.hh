// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_POISSON_HH
#define DUNE_PDELAB_POISSON_HH

#include<vector>

#if HAVE_TBB
#include "tbb/tbb_stddef.h"
#endif

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

#include<dune/geometry/type.hh>
#include<dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/localoperator/laplace.hh>

#include"defaultimp.hh"
#include"idefault.hh"
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
    template<typename F, typename B, typename J>
    class Poisson : public NumericalJacobianApplyVolume<Poisson<F,B,J> >,
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

      /** \brief Constructor
       *
       * \param quadOrder Order of the quadrature rule used for integrating over the element
       *
       * \note For multithreading evaluating f_, bctype_, and j_ must be
       *       thread safe.
       */
      Poisson (const F& f_, const B& bctype_, const J& j_, unsigned int quadOrder)
        : f(f_), bctype(bctype_), j(j_),
        laplace_(quadOrder),
        quadOrder_(quadOrder)
      {}

#if HAVE_TBB
      Poisson(Poisson &other, tbb::split) :
        f(other.f), bctype(other.bctype), j(other.j),
        laplace_(other.laplace_, tbb::split()),
        quadOrder_(other.quadOrder_)
      { }

      void join(Poisson &other)
      {
        laplace_.join(other.laplace_);
      }
#endif // HAVE_TBB

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
        laplace_.alpha_volume(eg, lfsu, x, lfsv, r);
	  }

      /** \brief Compute the Laplace stiffness matrix for the element given in 'eg'
       *
       * \tparam M Type of the element stiffness matrix
       *
       * \param [in]  eg The grid element we are assembling on
       * \param [in]  lfsu Local ansatz function space basis
       * \param [in]  lfsv Local test function space basis
       * \param [in]  x Current configuration; gets ignored for linear problems like this one
       * \param [out] matrix Element stiffness matrix
       */
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & matrix) const
      {
        laplace_.jacobian_volume(eg, lfsu, x, lfsv, matrix);
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
        typedef typename BasisSwitch::RangeField RF;
        typedef typename BasisSwitch::Range Range;

        // dimensions
        static const int dimLocal = EG::Geometry::mydimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dimLocal>& rule =
          Dune::QuadratureRules<DF,dimLocal>::rule(gt,quadOrder_);

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
          Dune::QuadratureRules<DF,dimLocal>::rule(gtface,quadOrder_);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dimLocal>::const_iterator it =
               rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate boundary condition type
            // skip rest if we are on Dirichlet boundary
            if( !bctype.isNeumann( ig,it->position() ) )
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

      // Laplace assembler to handle the matrix assembly
      Laplace laplace_;

      // Quadrature rule order
      unsigned int quadOrder_;
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
    template<typename Time, typename F, typename B, typename J>
    class InstationaryPoisson
      : public Poisson<F,B,J>,
        public InstationaryLocalOperatorDefaultMethods<Time>
    {
      typedef Poisson<F,B,J> Base;
      typedef InstationaryLocalOperatorDefaultMethods<Time> IDefault;

    protected:
      // need non-const references for setTime()
      F& f;
      B& bctype;
      J& j;

    public:
      //! construct InstationaryPoisson
      InstationaryPoisson(F& f_, B& bctype_, J& j_, unsigned int quadOrder)
        : Base(f_, bctype_, j_, quadOrder)
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
