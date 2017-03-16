// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/common/functionutilities.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>

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
class GenericAdvectionProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  GenericAdvectionProblem ()
  {
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 0 : 0;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      v[i] = 0.0;
    v[0] = 1.0;
    v[1] = 1.0/3.0;
    b0 = 0.25;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
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
    if (is.outerNormal(xlocal)*v<0)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType x = e.geometry().global(xlocal);
    if (x[1] < v[1]*x[0]+b0)
      return 0.0;
    else
      return 1.0;
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

private:
  typename Traits::PermTensorType I;
  typename Traits::RangeType v;
  typename Traits::RangeFieldType b0;
};

// Solve problem
template <typename Grid, typename FS, typename Problem, typename GM, Dune::SolverCategory::Category solvertype, int degree>
void solveProblem(const Grid& grid, FS& fs, typename FS::DOF x, Problem& problem, std::string basename)
{

  // initialize DOF vector it with boundary condition
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(grid->leafGridView(),problem);
  typedef typename Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(grid->leafGridView(),problem);
  Dune::PDELab::interpolate(g,fs.getGFS(),x);

  // assemble constraints
  fs.assembleConstraints(bctype);
  fs.setNonConstrainedDOFS(x,0.0);

  // assembler for finite elemenent problem
  typedef typename Dune::PDELab::ConvectionDiffusionDG<Problem,typename FS::FEM> LOP;
  LOP lop(problem,Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,2.0);
  typedef typename Dune::PDELab::GalerkinGlobalAssembler<FS,LOP,solvertype> ASSEMBLER;
  ASSEMBLER assembler(fs,lop,20);

  // make linear solver and solve problem
  typedef typename Dune::PDELab::ISTLSolverBackend_IterativeDefault<FS,ASSEMBLER,solvertype> SBE;
  SBE sbe(fs,assembler,5000,1);
  typedef typename Dune::PDELab::StationaryLinearProblemSolver<typename ASSEMBLER::GO,typename SBE::LS,typename FS::DOF> SLP;
  SLP slp(*assembler,*sbe,x,1e-6);
  slp.apply();

  // // print statistics about nonzero values per row
  // typename ASSEMBLER::MAT m(assembler.getGO());
  // std::cout << m.patternStatistics() << std::endl;

  // output grid to VTK file
  Dune::SubsamplingVTKWriter<typename GM::LeafGridView> vtkwriter(grid->leafGridView(),degree-1);
  typename FS::DGF xdgf(fs.getGFS(),x);
  vtkwriter.addVertexData(std::make_shared<typename FS::VTKF>(xdgf,"x_h"));
  vtkwriter.write(basename,Dune::VTK::appendedraw);
}


int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // command line args
  int cells=10; if (argc>=2) sscanf(argv[1],"%d",&cells);

  // define parameters
  const unsigned int dim = 2;
  const unsigned int degree = 1;
  const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
  const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::overlapping;
  typedef double NumberType;

  // make grid
  typedef Dune::YaspGrid<dim> GM;
  typedef Dune::PDELab::StructuredGrid<GM> Grid;
  Grid grid(elemtype,cells);
  grid->loadBalance();

  // make problem parameters
  typedef GenericAdvectionProblem<GM::LeafGridView,NumberType> Problem;
  Problem problem;

  // make a finite element space and DOF vector
  typedef Dune::PDELab::DGLegendreSpace<GM,NumberType,degree,elemtype,solvertype> FS;
  FS fs(grid->leafGridView());
  typedef typename FS::DOF X;
  X x(fs.getGFS(),0.0);

  // solve problem
  solveProblem<Grid,FS,Problem,GM,solvertype,degree>(grid,fs,x,problem,"new");

  // make a finite element space and DOF vector
  typedef Dune::PDELab::DGQkSpace<GM,NumberType,degree,elemtype,solvertype> FS2;
  FS2 fs2(grid->leafGridView());
  typedef typename FS2::DOF X2;
  X2 x2(fs2.getGFS(),0.0);

  // solve problem
  solveProblem<Grid,FS2,Problem,GM,solvertype,degree>(grid,fs2,x2,problem,"old");

  // calculate l2 error squared between the two functions
  FS::DGF xdgf(fs.getGFS(),x);
  FS2::DGF xdgf2(fs2.getGFS(),x2);
  typedef Dune::PDELab::DifferenceSquaredAdapter<FS::DGF,FS2::DGF> DifferenceSquared;
  DifferenceSquared differencesquared(xdgf,xdgf2);
  typename DifferenceSquared::Traits::RangeType l2errorsquared(0.0);
  Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,10);

  std::cout << std::endl << "l2 error squared: " << l2errorsquared << std::endl;

  if (l2errorsquared>1e-14){
    return 1;
  }

  // done
  return 0;
}
