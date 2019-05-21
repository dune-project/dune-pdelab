#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/gridfunctionspace/dunefunctionsgridfunctionspace.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>

/**
 * \page recipe-linear-system-assembly Assembling a linear system from a PDE
 *
 * This recipe shows how to assemble a linear system of equations discretizing
 * a PDE.
 *
 * In particular, we start from a LocalOperator representing the PDE's bilinear form. Rather
 * than defining one by ourselves, we use one already provided in PDELab. Further, we assume
 * that a GridFunctionSpace and Constraints have already been set up.
 * \snippet recipe-linear-system-assembly.cc LocalOperator
 *
 * This per-element LocalOperator can now be plugged into a GridOperator, which is capable
 * of assembling a matrix for the entire domain.
 * \snippet recipe-linear-system-assembly.cc Grid operator
 *
 * Now, from the GridOperator, we can assemble our residual vector depending
 * on our initial guess
 * \snippet recipe-linear-system-assembly.cc Residual assembly
 *
 * and print it to console.
 * \snippet recipe-linear-system-assembly.cc Residual output
 *
 * Likewise, we can assemble our matrix which again depends on the initial guess
 * \snippet recipe-linear-system-assembly.cc Matrix assembly
 * and print it as well.
 * \snippet recipe-linear-system-assembly.cc Matrix output

 * Full example code: @ref recipe-linear-system-assembly.cc
 * \example recipe-linear-system-assembly.cc
 * See explanation at @ref recipe-linear-system-assembly
 */


/** Parameter class for the stationary convection-diffusion equation of the following form:
 *
 * \f{align*}{
 *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \ \
 *                                              u &=& g \mbox{ on } \partial\Omega_D (Dirichlet)\ \
 *                (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N (Flux)\ \
 *                        -(A(x)\nabla u) \cdot n &=& o \mbox{ on } \partial\Omega_O (Outflow)
 * \f}
 * Note:
 *  - This formulation is valid for velocity fields which are non-divergence free.
 *  - Outflow boundary conditions should only be set on the outflow boundary
 *
 * The template parameters are:
 *  - GV a model of a GridView
 *  - RF numeric type to represent results
 */
template<typename GV, typename RF>
class GenericEllipticProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1 : 0;
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    return (2*GV::dimension-4*(xglobal*xglobal))*g(e,x);

    // return 0.0;
  }

  //! boundary condition type function
  /* return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet for Dirichlet boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann for flux boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow for outflow boundary conditions
   */
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType x = e.geometry().global(xlocal);
    return exp(-(x*x));
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }
};


int main(int argc, char **argv)
{
  try {
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc,argv);

    // define parameters
    typedef double NumberType;

    // need a grid in order to test grid functions
    constexpr unsigned int dim = 1;
    constexpr unsigned int degree = 1;
    constexpr std::size_t nonzeros = std::pow(2*degree+1,dim);

    Dune::FieldVector<NumberType,dim> L(1.0);
    std::array<int,dim> N(Dune::filledArray<dim,int>(5));

    typedef Dune::YaspGrid<dim> Grid;
    Grid grid(L,N);


    // make grid
    typedef Dune::YaspGrid<dim> GM;


    // make problem parameters
    typedef GenericEllipticProblem<typename GM::LeafGridView,NumberType> Problem;
    Problem problem;
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
    BCType bctype(grid.leafGridView(),problem);


    typedef typename GM::ctype DF;
    typedef Dune::PDELab::QkLocalFiniteElementMap<GM::LeafGridView,DF,NumberType,1> FEM;
    FEM fem(grid.leafGridView());

    // make function space, with overlapping constraints as usual for overlapping Schwarz methods
    typedef Dune::PDELab::GridFunctionSpace<GM::LeafGridView,FEM,
    Dune::PDELab::OverlappingConformingDirichletConstraints,
    Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,1>
    > GFS;
    GFS gfs(grid.leafGridView(),fem);

    typedef typename GFS::template ConstraintsContainer<NumberType>::Type CC;
    CC cc;
    Dune::PDELab::constraints(bctype,gfs,cc);

    // make a degree of freedom vector and initialize it with a function
    typedef Dune::PDELab::Backend::Vector<GFS,NumberType> X;
    X x(gfs,0.0);
    typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
    G g(grid.leafGridView(),problem);
    Dune::PDELab::interpolate(g,gfs,x);


    // assembler for finite elemenent problem
    //! [LocalOperator]
    typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
    LOP lop(problem);
    //! [LocalOperator]

    typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
    //! [Grid operator]
    typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC> GO;
    auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros));
    //! [Grid operator]


    // Assemble residual
    //! [Residual assembly]
    X d(gfs,0.0);
    go.residual(x,d);
    //! [Residual assembly]

    using Dune::PDELab::Backend::native;
    //! [Residual output]
    Dune::printvector(std::cout, native(d), "Residual vector", "");
    //! [Residual output]

    // Assemble jacobian
    //! [Matrix assembly]
    typedef GO::Jacobian M;
    M AF(go);
    go.jacobian(x,AF);
    //! [Matrix assembly]

    //! [Matrix output]
    Dune::printmatrix(std::cout, native(AF), "Stiffness matrix", "");
    //! [Matrix output]


    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
