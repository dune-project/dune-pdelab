// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_TEST_HELMHOLTZOPERATOR_HH
#define DUNE_PDELAB_TEST_HELMHOLTZOPERATOR_HH

#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/pattern.hh>

/** A local operator for solving the complex-valued Helmholtz equation
 *  \f{align*}{
 *          - \Delta u - \omega^2 u &=& f \mbox{in} \Omega \\
 *                                u &=& g \mbox{on} \Gamma_D \\
 *    \nabla u \cdot n - i \omega u &=& 0 \mbox{on} \Gamma_R
 *  \f}
 *
 * with conforming finite elements on all types of grids in any dimension
 *
 * \tparam PARAM Parameter class for this local operator.
 *
 * \author Philipp Stekl, Marian Piatkowski
 *
 */
template<class PARAM>
class HelmholtzLocalOperator :
  public Dune::PDELab::NumericalJacobianApplyVolume<HelmholtzLocalOperator<PARAM> >,
  public Dune::PDELab::NumericalJacobianVolume<HelmholtzLocalOperator<PARAM> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<HelmholtzLocalOperator<PARAM> >,
  public Dune::PDELab::NumericalJacobianBoundary<HelmholtzLocalOperator<PARAM> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doAlphaBoundary = true };

  HelmholtzLocalOperator (const PARAM& param_, unsigned int intorder_=2)
     : param(param_ ), intorder(intorder_)
  {
  }


  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // extract some types
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianType;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename LFSU::Traits::SizeType size_type;

    // dimensions
    const int dim = EG::Geometry::mydimension;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const auto& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // loop over quadrature points
    for (const auto& qp : rule) {
      // evaluate basis functions on reference element
      std::vector<RangeType> phi(lfsu.size());
      lfsu.finiteElement().localBasis().evaluateFunction(qp.position(),phi);

      // compute u at integration point
      RF u=0.0;
      for (size_type i=0; i<lfsu.size(); i++)
        u = u + x(lfsu,i)*phi[i];

      // evaluate gradient of basis functions on reference element
      std::vector<JacobianType> js(lfsu.size());
      lfsu.finiteElement().localBasis().evaluateJacobian(qp.position(),js);

      // transform gradients from reference element to real element
      const typename EG::Geometry::JacobianInverseTransposed
        jac = eg.geometry().jacobianInverseTransposed(qp.position());
      std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
      for (size_type i=0; i<lfsu.size(); i++)
        jac.mv(js[i][0],gradphi[i]);

      // compute gradient of u
      Dune::FieldVector<RF,dim> gradu(0.0);
      for (size_type i=0; i<lfsu.size(); i++)
        gradu.axpy(x(lfsu,i),gradphi[i]);

      // evaluate parameters
      RF f(param.f(eg.entity(),qp.position()));

      // integrate grad u * grad phi_i - omega*omega*u*phi_i - f phi_i
      RF factor = qp.weight()*eg.geometry().integrationElement(qp.position());
      for (size_type i=0; i<lfsu.size(); i++)
        r.accumulate(lfsu,i,( gradu*gradphi[i] - param.omega*param.omega*u*phi[i] - f*phi[i] )*factor);
    }
  }

  // - n \nabla u  = j
  // n \nabla u − i \omega u = 0
  // boundary integral
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
                       const LFSV& lfsv_s, R& r_s) const
  {
    // some types
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType Range;
    typedef typename LFSU::Traits::SizeType size_type;

    // dimensions
    const int dim = IG::dimension;

    Dune::GeometryType gtface = ig.geometryInInside().type();

    // select quadrature rule for face
    const auto& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

    // loop over quadrature points and integrate normal flux
    for (const auto& qp : rule) {
      // skip rest if we are on Dirichlet boundary
      if ( param.isDirichlet(ig, qp.position()))
        continue;

      // position of quadrature point in local coordinates of element
      Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(qp.position());

      // evaluate basis functions at integration point
      std::vector<Range> phi(lfsu_s.size());
      lfsu_s.finiteElement().localBasis().evaluateFunction(local,phi);

      // evaluate u (e.g. flux may depend on u)
      RF u=0.0;
      for (size_type i=0; i<lfsu_s.size(); ++i)
        u = u + x_s(lfsu_s,i)*phi[i];

      //============================================
      // NOTE
      // We treat Robin boundary conditions as Neumann boundary conditions
      // by making the flux also depend on u.
      //============================================
      if ( param.isNeumann(ig.intersection(), qp.position()))
        {
          // evaluate flux boundary condition
          RF j = param.j(ig.intersection(), qp.position(), u);

          // integrate j
          RF factor = qp.weight()*ig.geometry().integrationElement(qp.position());
          for (size_type i=0; i<lfsu_s.size(); ++i)
            r_s.accumulate(lfsu_s,i,j*phi[i]*factor);
        }
    }
  }

private:

  const PARAM& param;
  unsigned int intorder;
};
#endif // DUNE_PDELAB_TEST_HELMHOLTZOPERATOR_HH
