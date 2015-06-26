#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/pattern.hh>

/** a local operator for solving the equation
 *
 *   - \Delta u + a*u = f   in \Omega
 *                  u = g   on \Gamma_D\subseteq\partial\Omega
 *  -\nabla u \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 * with conforming finite elements on all types of grids in any dimension
 *
 * \tparam BCType parameter class indicating the type of boundary condition
 */
template<class BCType>
class Example02LocalOperator :
  public Dune::PDELab::NumericalJacobianApplyVolume<Example02LocalOperator<BCType> >,
  public Dune::PDELab::NumericalJacobianVolume<Example02LocalOperator<BCType> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<Example02LocalOperator<BCType> >,
  public Dune::PDELab::NumericalJacobianBoundary<Example02LocalOperator<BCType> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doAlphaBoundary = true };                                // assemble boundary

  Example02LocalOperator(const BCType& bctype_, // boundary cond.type
                         unsigned int intorder_=2) :
    bctype( bctype_ ), intorder( intorder_ )
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // assume Galerkin: lfsu == lfsv
    // This yields more efficient code since the local functionspace only
    // needs to be evaluated once, but would be incorrect for a finite volume
    // method

    // dimensions
    const int dim = EG::Geometry::dimension;
    const int dimw = EG::Geometry::dimensionworld;

    // extract some types
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType Jacobian;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType Range;
    typedef Dune::FieldVector<RF,dimw> Gradient;
    typedef typename LFSU::Traits::SizeType size_type;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF,dim>&
      rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator
           it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate basis functions on reference element
        std::vector<Range> phi(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

        // compute u at integration point
        RF u=0.0;
        for (size_type i=0; i<lfsu.size(); ++i)
          u += x(lfsu,i)*phi[i];

        // evaluate gradient of basis functions on reference element
        std::vector<Jacobian> js(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);

        // transform gradients from reference element to real element
        const typename EG::Geometry::JacobianInverseTransposed
          jac = eg.geometry().jacobianInverseTransposed(it->position());
        std::vector<Gradient> gradphi(lfsu.size());
        for (size_type i=0; i<lfsu.size(); i++)
          jac.mv(js[i][0],gradphi[i]);

        // compute gradient of u
        Gradient gradu(0.0);
        for (size_type i=0; i<lfsu.size(); ++i)
          gradu.axpy(x(lfsu,i),gradphi[i]);

        // evaluate parameters;
        // Dune::FieldVector<RF,dim>
        //   globalpos = eg.geometry().global(it->position());
        RF f = 0;
        RF a = 0;

        // integrate grad u * grad phi_i + a*u*phi_i - f phi_i
        RF factor = it->weight()*eg.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsu.size(); ++i)
          r.accumulate(lfsu, i, (gradu*gradphi[i] + a*u*phi[i] - f*phi[i]) * factor);
      }
  }

  // boundary integral
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
                       const LFSV& lfsv_s, R& r_s) const
  {
    // assume Galerkin: lfsu_s == lfsv_s
    // This yields more efficient code since the local functionspace only
    // needs to be evaluated once, but would be incorrect for a finite volume
    // method

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

    // select quadrature rule for face
    Dune::GeometryType gtface = ig.geometryInInside().type();
    const Dune::QuadratureRule<DF,dim-1>&
      rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

    // loop over quadrature points and integrate normal flux
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin();
         it!=rule.end(); ++it)
      {
        // skip rest if we are on Dirichlet boundary
        if ( bctype.isDirichlet( ig, it->position() ) )
          continue;

        // position of quadrature point in local coordinates of element
        Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

        // evaluate basis functions at integration point
        std::vector<Range> phi(lfsu_s.size());
        lfsu_s.finiteElement().localBasis().evaluateFunction(local,phi);

        // evaluate u (e.g. flux may depend on u)
        RF u=0.0;
        for (size_type i=0; i<lfsu_s.size(); ++i)
          u += x_s(lfsu_s,i)*phi[i];

        // evaluate flux boundary condition
        Dune::FieldVector<RF,dim>
          globalpos = ig.geometry().global(it->position());
        RF j;
        if (globalpos[1]<0.5)
          j = 1.0;
        else
          j = -1.0; // some outflow

        // integrate j
        RF factor = it->weight()*ig.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsu_s.size(); ++i)
          r_s.accumulate(lfsu_s,i,j*phi[i]*factor);
      }
  }

private:
  const BCType& bctype;
  unsigned int intorder;
};
