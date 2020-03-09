#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>

/**
 * \page recipe-linear-system-solution-pdelab Solving a linear system within PDELab
 *
 * Here we show how to solve a linear system of equations originating from a PDE using PDELab.
 *
 * First, we set up a GridOperator as in @ref recipe-linear-system-assembly
 * \snippet recipe-linear-system-solution-pdelab.cc Grid operator
 *
 * Next, we set up our degree of freedom vector
 * \snippet recipe-linear-system-solution-pdelab.cc Make degree of freedom vector
 *
 * and ensure it matches the Dirichlet boundary conditions at constrained degrees of freedom.
 * In addition to specifying Dirichlet constrained degrees of freedom, it also serves as
 * initial guess at unconstrained ones.
 * \snippet recipe-linear-system-solution-pdelab.cc Set it to match boundary conditions
 *
 * Now we choose the preconditioner and solver we want to use
 * \snippet recipe-linear-system-solution-pdelab.cc Linear solver backend
 *
 * and plug it into a StationaryLinearProblemSolver. This takes care of assembling
 * as well as solving the system.
 * \snippet recipe-linear-system-solution-pdelab.cc Assemble and solve linear problem
 *
 * Finally, let's print the result to console via
 * \snippet recipe-linear-system-solution-pdelab.cc Solution output
 *
 * There is a number of alternative solvers and preconditioners available we could use instead, for example this one:
 * \snippet recipe-linear-system-solution-pdelab.cc Linear solver backend 2
 *
 * Full example code: @ref recipe-linear-system-solution-pdelab.cc
 * \example recipe-linear-system-solution-pdelab.cc
 * See explanation at @ref recipe-linear-system-solution-pdelab
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
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1 : 0;
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    return 0.0;
  }

  //! boundary condition type function
  /* return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet for Dirichlet boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann for flux boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow for outflow boundary conditions
   */
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType x = e.geometry().global(xlocal);
    if (x[0] < 0.5)
      return 1.0;

    return 0.0;
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
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
    constexpr std::size_t nonzeros = Dune::power(2*degree+1,dim);

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

    // make function space with constraints
    typedef Dune::PDELab::GridFunctionSpace<GM::LeafGridView,FEM,
    Dune::PDELab::ConformingDirichletConstraints,
    Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,1>> GFS;
    GFS gfs(grid.leafGridView(),fem);

    typedef typename GFS::template ConstraintsContainer<NumberType>::Type CC;
    CC cc;
    Dune::PDELab::constraints(bctype,gfs,cc);

    // assembler for finite elemenent problem
    typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
    LOP lop(problem);

    typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
    // [Grid operator]
    typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC> GO;
    auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros));
    //! [Grid operator]

    // [Make degree of freedom vector]
    typedef Dune::PDELab::Backend::Vector<GFS,NumberType> X;
    X x(gfs,0.0);
    //! [Make degree of freedom vector]
    // [Set it to match boundary conditions]
    typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
    G g(grid.leafGridView(),problem);
    Dune::PDELab::interpolate(g,gfs,x);
    //! [Set it to match boundary conditions]

    // [Linear solver backend]
    typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
    LS ls(100,3);
    //! [Linear solver backend]

    {
    // [Linear solver backend 2]
    typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
    //! [Linear solver backend 2]
    LS ls(100,3);
    }

    // [Assemble and solve linear problem]
    typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,X> SLP;
    SLP slp(go,ls,x,1e-10);
    slp.apply(); // here all the work is done!
    //! [Assemble and solve linear problem]

    using Dune::PDELab::Backend::native;
    // [Solution output]
    Dune::printvector(std::cout, native(x), "Solution", "");
    //! [Solution output]

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
