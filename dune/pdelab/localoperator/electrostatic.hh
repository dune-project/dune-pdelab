// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_ELECTROSTATIC_HH
#define DUNE_PDELAB_ELECTROSTATIC_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/geometrytype.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/typetraits.hh>
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

    /** a local operator for solving the electrostatic "wave" equation
     *
     * \f{align*}{
     *    \nabla\times\mu^{-1}\nabla\mathbf E &=0\qquad\text{in }\Omega \\
     *    \mathbf{\hat n}\times(\mathbf{\hat n}\times\mathbf E) &= \mathbf g\qquad\text{on }\partial\Omega
     * \f}
     * or in weak formulation
     * \f[
     *    \int_\Omega\mu^{-1}(\nabla\times\mathbf v)\cdot(\nabla\times\mathbf E)dV
     *    +\oint_{\partial\Omega}\mu^{-1}\mathbf v\cdot(\mathbf{\hat n}\times\nabla\times\mathbf E)dS
     *    = 0 \qquad \forall\mathbf v
     * \f]
     * By expanding with the basis we get
     * \f[
     *    \sum_iE_i\int_\Omega\mu^{-1}(\nabla\times\phi_j)\cdot(\nabla\times\phi_i)dV
     *    +\sum_iE_i\oint_{\partial\Omega}\mu^{-1}\phi_j\cdot(\mathbf{\hat n}\times\nabla\times\phi_i)dS
     *    = 0 \qquad \forall j
     * \f]
     * with edge elements on all grids and any dimension
     * \tparam Mu grid function type giving \f$\mu\f$
     */
    template<typename Mu, int qorder=1>
	class Electrostatic
      : public NumericalJacobianApplyVolume<Electrostatic<Mu, qorder> >
      , public NumericalJacobianVolume<Electrostatic<Mu, qorder> >
      , public NumericalJacobianApplyBoundary<Electrostatic<Mu, qorder> >
      , public NumericalJacobianBoundary<Electrostatic<Mu, qorder> >
      , public FullVolumePattern
	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = false };

	  // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doAlphaSkeleton = false };
      enum { doAlphaBoundary = true };
      enum { doLambdaVolume = false };
      enum { doLambdaSkeleton = false };
      enum { doLambdaBoundary = false };

      Electrostatic (const Mu& mu_)
        : mu(mu_)
      {}

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFS, typename X, typename R>
	  void alpha_volume (const EG& eg, const LFS& lfsu, const X& x, const LFS& lfsv, R& r) const
	  {
        // We support only Galerkin method lfsu==lfsv
		dune_static_assert(LFS::Traits::LocalFiniteElementType::
                           Traits::LocalBasisType::Traits::dimRange == 3,
                           "Works only in 3D");

		// domain and range field type
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
        dune_static_assert((Dune::is_same<typename EG::Geometry::ctype, DF>::value),
                           "Grids ctype and Finite Elements DomainFieldType must match");
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainType D;
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> J(lfsu.size());
            lfsu.localFiniteElement().localBasis().evaluateJacobian(it->position(),J);

            std::vector<D> rotphi(lfsu.size(),D(0));
            for(unsigned i = 0; i < lfsu.size(); ++i)
              for(unsigned j = 0; j < 3; ++j)
                rotphi[i] += J[i][(j+2)%3][(j+1)%3] - J[i][(j+1)%3][(j+2)%3];

            Dune::FieldVector<RF,1> muval;
            mu.evaluate(eg.entity(), it->position(), muval);

            RangeType rotE(0);
            for(unsigned i = 0; i < lfsu.size(); ++i)
              rotE.axpy(x[i], rotphi[i]);
            
            // integrate grad u * grad phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position()) / muval;
            for (size_t j=0; j<lfsu.size(); j++)
              r[j] += (rotphi[j]*rotE)*factor;
          }
	  }

	  // boundary integral depending on test and ansatz functions
	  template<typename IG, typename LFS, typename X, typename R>
	  void alpha_boundary (const IG& ig, const LFS& lfsu_s, const X& x_s,
                           const LFS& lfsv_s, R& r_s) const
	  {
        // We support only Galerkin method lfsu=lfsv
		dune_static_assert(LFS::Traits::LocalFiniteElementType::
                           Traits::LocalBasisType::Traits::dimRange == 3,
                           "Works only in 3D");

		// domain and range field type
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainType D;
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;

        // dimensions
        const int dim = IG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = ig.geometry().type();
        const Dune::QuadratureRule<DF,dim-1>& rule =
          Dune::QuadratureRules<DF,dim-1>::rule(gt,qorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin();
             it!=rule.end(); ++it)
        {
          D pos = ig.geometryInInside().global(it->position());

          // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
          std::vector<JacobianType> J(lfsu_s.size());
          lfsu_s.localFiniteElement().localBasis().evaluateJacobian(pos,J);

          std::vector<D> rotphi(lfsu_s.size(),D(0));
          for(unsigned i = 0; i < lfsu_s.size(); ++i)
            for(unsigned j = 0; j < 3; ++j)
              rotphi[i] += J[i][(j+2)%3][(j+1)%3] - J[i][(j+1)%3][(j+2)%3];

          RangeType rotE(0);
          for(unsigned i = 0; i < lfsu_s.size(); ++i)
            rotE.axpy(x_s[i], rotphi[i]);
            
          RangeType n = ig.unitOuterNormal(it->position());
          RangeType nxrotE(0);
          for(unsigned i = 0; i < 3; ++i)
            nxrotE += n[(i+1)%3]*rotE[(i+2)%3] - n[(i+2)%3]*rotE[(i+1)%3];

          Dune::FieldVector<RF,1> muval;
          mu.evaluate(*ig.inside(), pos, muval);

          std::vector<D> phi(lfsu_s.size());
          lfsu_s.localFiniteElement().localBasis().evaluateFunction(pos,phi);

          // integrate v * (n x rot E)
          RF factor = it->weight() * ig.geometry().integrationElement(it->position()) / muval;
          for (size_t j=0; j<lfsu_s.size(); j++)
            r_s[j] += (phi[j]*nxrotE)*factor;
        }
      }
    private:
      const Mu &mu;
	};

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ELECTROSTATIC_HH
