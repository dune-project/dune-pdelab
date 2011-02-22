// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_NPDE_CONVECTIONDIFFUSIONFEM_HH
#define DUNE_NPDE_CONVECTIONDIFFUSIONFEM_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/grid/common/quadraturerules.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>

#include"convectiondiffusionparameter.hh"

/** a local operator for solving convection-diffusion equation with standard FEM
 *  
 * \f{align*}{
 *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \\
 *                                              u &=& g \mbox{ on } \partial\Omega_D \\
 *                (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
 *                        -(A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_O
 * \f}
 * Note:
 *  - This formulation is valid for velocity fields which are non-divergence free.
 *  - Outflow boundary conditions should only be set on the outflow boundary
 *
 * \tparam T model of ConvectionDiffusionParameterInterface
 */
template<typename T>
class ConvectionDiffusionFEM : 
  public Dune::PDELab::NumericalJacobianApplyVolume<ConvectionDiffusionFEM<T> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionFEM<T> >,
  public Dune::PDELab::NumericalJacobianVolume<ConvectionDiffusionFEM<T> >,
  public Dune::PDELab::NumericalJacobianBoundary<ConvectionDiffusionFEM<T> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doAlphaBoundary = true };

  ConvectionDiffusionFEM (T& param_, int intorderadd_=0) : param(param_), intorderadd(intorderadd_)
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

    typedef typename LFSU::Traits::SizeType size_type;
        
    // dimensions
    const int dim = EG::Geometry::dimension;
    const int dimw = EG::Geometry::dimensionworld;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const int intorder = intorderadd+2*lfsu.finiteElement().localBasis().order();
    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // evaluate diffusion tensor at cell center, assume it is constant over elements
    typename T::Traits::PermTensorType tensor;
    Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
    tensor = param.A(eg.entity(),localcenter);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate basis functions
        std::vector<RangeType> phi(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

        // evaluate u
        RF u=0.0;
        for (size_type i=0; i<lfsu.size(); i++)
          u += x[lfsu.localIndex(i)]*phi[i];

        // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
        std::vector<JacobianType> js(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);

        // transform gradients of shape functions to real element
        const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
        for (size_type i=0; i<lfsu.size(); i++)
          jac.mv(js[i][0],gradphi[i]);

        // compute gradient of u
        Dune::FieldVector<RF,dim> gradu(0.0);
        for (size_type i=0; i<lfsu.size(); i++)
          gradu.axpy(x[lfsu.localIndex(i)],gradphi[i]);

        // compute A * gradient of u
        Dune::FieldVector<RF,dim> Agradu(0.0);
        tensor.umv(gradu,Agradu);

        // evaluate velocity field, sink term and source te
        typename T::Traits::RangeType b = param.b(eg.entity(),it->position());
        typename T::Traits::RangeFieldType c = param.c(eg.entity(),it->position());
        typename T::Traits::RangeFieldType f = param.f(eg.entity(),it->position());

        // integrate (A grad u)*grad phi_i - u b*grad phi_i + c*u*phi_i
        RF factor = it->weight() * eg.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsu.size(); i++)
          r[lfsu.localIndex(i)] += ( Agradu*gradphi[i] - u*(b*gradphi[i]) + (c*u-f)*phi[i] )*factor;
      }
  }

  // jacobian of volume term
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
                        Dune::PDELab::LocalMatrix<R>& mat) const
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
    typedef typename LFSU::Traits::SizeType size_type;

    // dimensions
    const int dim = EG::Geometry::dimension;
    const int dimw = EG::Geometry::dimensionworld;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const int intorder = intorderadd+2*lfsu.finiteElement().localBasis().order();
    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // evaluate diffusion tensor at cell center, assume it is constant over elements
    typename T::Traits::PermTensorType tensor;
    Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
    tensor = param.A(eg.entity(),localcenter);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
        std::vector<JacobianType> js(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);

        // transform gradient to real element
        const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
        std::vector<Dune::FieldVector<RF,dim> > Agradphi(lfsu.size());
        for (size_type i=0; i<lfsu.size(); i++)
          {
            jac.mv(js[i][0],gradphi[i]);
            tensor.mv(gradphi[i],Agradphi[i]);
          }

        // evaluate basis functions
        std::vector<RangeType> phi(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

        // evaluate velocity field, sink term and source te
        typename T::Traits::RangeType b = param.b(eg.entity(),it->position());
        typename T::Traits::RangeFieldType c = param.c(eg.entity(),it->position());

        // integrate (A grad phi_j)*grad phi_i - phi_j b*grad phi_i + c*phi_j*phi_i
        RF factor = it->weight() * eg.geometry().integrationElement(it->position());
        for (size_type j=0; j<lfsu.size(); j++)
          for (size_type i=0; i<lfsu.size(); i++)
            mat(lfsu.localIndex(i),lfsu.localIndex(j)) += ( Agradphi[j]*gradphi[i]-phi[j]*(b*gradphi[i])+c*phi[j]*phi[i] )*factor;
      }
  }

  // boundary integral
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, 
                       const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                       R& r_s) const
  {
    // domain and range field type
    typedef typename LFSV::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSV::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSV::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;

    typedef typename LFSV::Traits::SizeType size_type;
        
    // dimensions
    const int dim = IG::dimension;
        
    // evaluate boundary condition type
    Dune::GeometryType gtface = ig.geometryInInside().type();
    Dune::FieldVector<DF,dim-1> facecenterlocal = Dune::GenericReferenceElements<DF,dim-1>::general(gtface).position(0,0);
    ConvectionDiffusionBoundaryConditions::Type bctype;
    bctype = param.bctype(ig.intersection(),facecenterlocal);
 
    // skip rest if we are on Dirichlet boundary
    if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet) return;

    // select quadrature rule
    const int intorder = intorderadd+2*lfsu_s.finiteElement().localBasis().order();
    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

    // loop over quadrature points and integrate normal flux
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // position of quadrature point in local coordinates of element 
        Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

        // evaluate shape functions (assume Galerkin method) 
        std::vector<RangeType> phi(lfsu_s.size());
        lfsu_s.finiteElement().localBasis().evaluateFunction(local,phi);

        if (bctype==ConvectionDiffusionBoundaryConditions::Neumann)
          {
            // evaluate flux boundary condition
            typename T::Traits::RangeFieldType j = param.j(ig.intersection(),it->position());
            
            // integrate j
            RF factor = it->weight()*ig.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s[lfsu_s.localIndex(i)] += j*phi[i]*factor;
          }

        if (bctype==ConvectionDiffusionBoundaryConditions::Outflow)
          {
            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              u += x_s[lfsu_s.localIndex(i)]*phi[i];

            // evaluate velocity field and outer unit normal
            typename T::Traits::RangeType b = param.b(*(ig.inside()),local);
            const Dune::FieldVector<DF,dim> n = ig.unitOuterNormal(it->position());

            // evaluate outflow boundary condition
            typename T::Traits::RangeFieldType o = param.o(ig.intersection(),it->position());
            
            // integrate o
            RF factor = it->weight()*ig.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s[lfsu_s.localIndex(i)] += ( (b*n)*u + o)*phi[i]*factor;
          }
      }
  }

  // jacobian contribution from boundary
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_boundary (const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          Dune::PDELab::LocalMatrix<R>& mat_s) const
  {
    // domain and range field type
    typedef typename LFSV::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSV::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSV::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;

    typedef typename LFSV::Traits::SizeType size_type;
        
    // dimensions
    const int dim = IG::dimension;
        
    // evaluate boundary condition type
    Dune::GeometryType gtface = ig.geometryInInside().type();
    Dune::FieldVector<DF,dim-1> facecenterlocal = Dune::GenericReferenceElements<DF,dim-1>::general(gtface).position(0,0);
    ConvectionDiffusionBoundaryConditions::Type bctype;
    bctype = param.bctype(ig.intersection(),facecenterlocal);
 
    // skip rest if we are on Dirichlet boundary
    if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet) return;
    if (bctype==ConvectionDiffusionBoundaryConditions::Neumann) return;

    // select quadrature rule
    const int intorder = intorderadd+2*lfsu_s.finiteElement().localBasis().order();
    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

    // loop over quadrature points and integrate normal flux
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // position of quadrature point in local coordinates of element 
        Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

        // evaluate shape functions (assume Galerkin method) 
        std::vector<RangeType> phi(lfsu_s.size());
        lfsu_s.finiteElement().localBasis().evaluateFunction(local,phi);

        // evaluate velocity field and outer unit normal
        typename T::Traits::RangeType b = param.b(*(ig.inside()),local);
        const Dune::FieldVector<DF,dim> n = ig.unitOuterNormal(it->position());
        
        // integrate 
        RF factor = it->weight()*ig.geometry().integrationElement(it->position());
        for (size_type j=0; j<lfsu_s.size(); j++)
          for (size_type i=0; i<lfsu_s.size(); i++)
            mat_s(lfsu_s.localIndex(i),lfsu_s.localIndex(j)) += (b*n)*phi[j]*phi[i]*factor;
      }
  }


  //! set time in parameter class
  void setTime (double t)
  {
    param.setTime(t);
  }

private:
  T& param;
  int intorderadd;
};

#endif
