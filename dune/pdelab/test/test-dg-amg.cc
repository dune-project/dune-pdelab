#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/istl/matrixmatrix.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include <dune/pdelab/localoperator/l2.hh>

#include <dune/pdelab/backend/istl.hh>

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
    return 0.0;
  }

  //! boundary condition type function
  /* return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet for Dirichlet boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann for flux boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow for outflow boundary conditions
   */
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
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
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // command line args
  int cells=16; if (argc>=2) sscanf(argv[1],"%d",&cells);
  int refinements=0; if (argc>=3) sscanf(argv[2],"%d",&refinements);

  Dune::Timer watch;

  // define parameters
  const unsigned int dim = 2;
  const unsigned int degree = 1;
  const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
  const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::conforming;
  const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::overlapping;
  typedef double NumberType;

  // make grid
  typedef Dune::YaspGrid<dim> GM;
  typedef Dune::PDELab::StructuredGrid<GM> Grid;
  Grid grid(elemtype,cells,1);
  grid->loadBalance();
  grid->refineOptions(false);
  grid->globalRefine(refinements);

  // make problem parameters
  typedef GenericEllipticProblem<GM::LeafGridView,NumberType> Problem;
  Problem problem;
  // make dummy constraints parameters, DG has no essential boundary conditions
  typedef Dune::PDELab::DirichletConstraintsParameters BCType;
  BCType bctype;

  // make DG finite element space
  typedef Dune::PDELab::DGQkSpace<GM,NumberType,degree,elemtype,solvertype> FS;
  FS fs(grid->leafGridView());
  fs.assembleConstraints(bctype);
  //std::cout << "number of constraints is " << fs.getCC().size() << std::endl;

  // assembler for finite elemenent problem
  typedef Dune::PDELab::ConvectionDiffusionDG<Problem,typename FS::FEM> LOP;
  LOP lop(problem,Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,2.0);
  typedef Dune::PDELab::GalerkinGlobalAssemblerNewBackend<FS,LOP,solvertype> ASSEMBLER;
  ASSEMBLER assembler(fs,lop,ASSEMBLER::MBE(5)); // 5 entries per row with cartesian mesh in 2D and blocked DG space

  // allocate solution vector
  typedef FS::DOF V;
  V x(fs.getGFS(),0.0);
  std::cout << "number of elements is " << Dune::PDELab::Backend::native(x).N() << std::endl;

  // use purely Neumann-zero boundary conditions
  // projected matrix from DG space has no essential boundary conditions
  typedef Dune::PDELab::NoDirichletConstraintsParameters CGBCType;
  CGBCType cgbctype;
  // CG space
  typedef Dune::PDELab::CGSpace<GM,NumberType,1,CGBCType,elemtype,meshtype,solvertype> CGFS;
  CGFS cgfs(*grid,cgbctype);
  cgfs.assembleConstraints(cgbctype);

  // need to set up my own grid operator that does not use P0Constraints
  typedef typename FS::GFS GFS;
  typedef typename FS::CC DGCC2;
  DGCC2 dgcc2; // empty: no constraints!
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,DGCC2,DGCC2> DGGO2;
  DGGO2 dggo2(fs.getGFS(),dgcc2,fs.getGFS(),dgcc2,lop,MBE(5));

  /////////////////// SEQUENTIAL
  // make linear solver and solve problem
  {
      typedef Dune::PDELab::ISTLBackend_SEQ_AMG_4_DG<ASSEMBLER::GO,CGFS::GFS,
                                                     Dune::PDELab::CG2DGProlongation,Dune::SeqSSOR,Dune::CGSolver> LS;
      LS ls(assembler.getGO(),cgfs.getGFS(),1000,3);
      // set parameters for AMG in CG-subspace
      Dune::Amg::Parameters params = ls.parameters();
      params.setCoarsenTarget(2000);
      params.setMaxLevel(20);
      params.setProlongationDampingFactor(1.8);
      params.setNoPreSmoothSteps(2);
      params.setNoPostSmoothSteps(2);
      params.setGamma(1);
      params.setAdditive(false);
      ls.setParameters(params);
      typedef Dune::PDELab::StationaryLinearProblemSolver<DGGO2,LS,V> SLP;
      SLP slp(dggo2,ls,x,1e-8);
      slp.setHangingNodeModifications(false);
      slp.apply();
  }

  /////////////////// OVERLAPPING
  // reset initial iterate
  x = 0.0;
  // make linear solver and solve problem
#if HAVE_MPI
  {
      typedef Dune::PDELab::ISTLBackend_OVLP_AMG_4_DG<ASSEMBLER::GO,FS::CC,CGFS::GFS,CGFS::CC,
                                                      Dune::PDELab::CG2DGProlongation,Dune::SeqSSOR,Dune::CGSolver> LS;
      LS ls(assembler.getGO(),fs.getCC(),cgfs.getGFS(),cgfs.getCC(),1000,3);
      // set parameters for AMG in CG-subspace
      Dune::Amg::Parameters params = ls.parameters();
      params.setCoarsenTarget(2000);
      params.setMaxLevel(20);
      params.setProlongationDampingFactor(1.8);
      params.setNoPreSmoothSteps(2);
      params.setNoPostSmoothSteps(2);
      params.setGamma(1);
      params.setAdditive(false);
      ls.setParameters(params);
      typedef Dune::PDELab::StationaryLinearProblemSolver<DGGO2,LS,V> SLP;
      SLP slp(dggo2,ls,x,1e-8);
      slp.setHangingNodeModifications(false);
      slp.apply();
  }
#endif // HAVE_MPI

  // done
  return 0;
}
