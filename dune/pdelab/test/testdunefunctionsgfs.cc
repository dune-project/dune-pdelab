#include "config.h"

#include <iostream>
#include <vector>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>

#include <dune/pdelab/adaptivity/adaptivity.hh>
#include <dune/pdelab/gridfunctionspace/dunefunctionsgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/localoperator/linearelasticity.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/variablefactories.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/newton/newton.hh>


using namespace Dune;

//===============================================================
//===============================================================
// Solve the Poisson equation
//           - \Delta u = f in \Omega,
//                    u = g on \partial\Omega_D
//  -\nabla u \cdot \nu = j on \partial\Omega_N
//===============================================================
//===============================================================

template <class GridView, class RangeType>
class PoissonProblem
: public PDELab::ConvectionDiffusionModelProblem<GridView,RangeType>
{
public:

  template<typename Element, typename Coord>
  auto f(const Element& element, const Coord& x) const
  {
    return 1.0;
  }

  template<typename Element, typename Coord>
  auto bctype(const Element& element, const Coord& x) const
  {
    return PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet extension
  template<typename E, typename X>
  RangeType g (const E& e, const X& x) const
  {
    auto global = e.geometry().global(x);
    RangeType s=0.0;
    for (std::size_t i=0; i<global.size(); i++)
      s+=global[i]*global[i];
    return s;
  }

};

template <int order>
void solvePoissonProblem()
{
  //////////////////////////////////
  //   Generate the grid
  //////////////////////////////////

  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {10, 10};
  GridType grid(l,elements);

  ////////////////////////////////////////////////////////
  //  Assemble the algebraic system
  ////////////////////////////////////////////////////////

  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  MatrixType stiffnessMatrix;

  // Construct Lagrangian finite element space basis
  using GridView = typename GridType::LeafGridView;
  auto gridView = grid.leafGridView();
  using Basis = Functions::LagrangeBasis<GridView,order>;
  auto basis = std::make_shared<Basis>(gridView);

  // What precisely does this do?
  typedef PDELab::ConformingDirichletConstraints Constraints;
  auto con = std::make_shared<Constraints>();

  using VectorBackend = PDELab::ISTL::VectorBackend<>;

  typedef PDELab::Experimental::GridFunctionSpace<Basis,VectorBackend,Constraints> GridFunctionSpace;
  GridFunctionSpace gfs(basis,con);

  // Test the 'size' method
  if (order==1 and gfs.size() != gridView.size(dim))
      DUNE_THROW(Exception, "gfs.size() does not return the correct number!");
  assert(gfs.size() == gfs.globalSize());
  assert(gfs.size() == gfs.blockCount());

  // Container for the Dirichlet boundary conditions
  typedef typename GridFunctionSpace::template ConstraintsContainer<double>::Type C;
  C constraintsContainer;

  PoissonProblem<GridView,double> problem;
  PDELab::ConvectionDiffusionBoundaryConditionAdapter<decltype(problem)> bctype(problem);
  PDELab::constraints(bctype,gfs,constraintsContainer);


  // make grid operator
  typedef PDELab::ConvectionDiffusionFEM<decltype(problem),typename GridFunctionSpace::Traits::FiniteElementMap> LOP;
  LOP lop(problem);

  Dune::PDELab::Backend::Vector<GridFunctionSpace,double> v(gfs,0);

  typedef PDELab::GridOperator<GridFunctionSpace,
                               GridFunctionSpace,
                               LOP,
                               PDELab::ISTL::BCRSMatrixBackend<>,
                               double,double,double,C,C> GO;
  GO go(gfs,constraintsContainer,gfs,constraintsContainer,lop, {9});

  // Pseudo "current" coefficient vector, not actually used
  typedef typename GO::Traits::Domain VectorContainer;
  VectorContainer x0(gfs);
  x0 = 0.0;              // set all entries to zero

  // Make a grid function out of it, just to check that this compiles
  PDELab::DiscreteGridFunction<GridFunctionSpace,VectorContainer> xDiscreteGridFunction(gfs,x0);

  // represent operator as a matrix
  typedef typename GO::Jacobian MatrixContainer;
  MatrixContainer m(go,stiffnessMatrix);    // Use the stiffnessMatrix object for the actual storage
  m = 0.0;             // set all matrix entries to zero

  go.jacobian(x0,m);      // assemble the stiffness matrix

  // evaluate residual w.r.t pseudo "current" iterate
  using VectorType = PDELab::Backend::Native<PDELab::Backend::Vector<GridFunctionSpace,double> >;
  VectorType rhs;
  VectorContainer r(gfs,rhs);  //Use the rhs object for the actual storage
  r = 0.0;

  go.residual(x0,r);

  ///////////////////////////
  //   Compute solution
  ///////////////////////////

  VectorType x(rhs.size());
  x = 0;

  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

  // A preconditioner
  SeqILU<MatrixType,VectorType,VectorType> ilu0(stiffnessMatrix,1.0);

  // A preconditioned conjugate-gradient solver
  CGSolver<VectorType> cg(op,ilu0,1E-4,
                          50,   // maximum number of iterations
                          2);

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  cg.apply(x, rhs, statistics);

  // Test whether we can refine the grid and carry a function along.
  // Of course we cannot: we haven't marked anything, and we are using
  // YaspGrid anyway.  But let's make sure we can at least call the method.
  VectorContainer xContainer(gfs,x);
  PDELab::adapt_grid(grid, gfs, xContainer, 2 );

  // Output result to VTK file
  auto pressureFunction = Functions::makeDiscreteGlobalBasisFunction<double>(*basis,x);

  SubsamplingVTKWriter<GridView> vtkWriter(gridView, refinementLevels(2));
  vtkWriter.addVertexData(pressureFunction, VTK::FieldInfo("pressure", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("testdunefunctionsgfs-poisson");

}

void solveParallelPoissonProblem()
{
  // read ini file
  const int dim = 2;
  const int refinement = 3;
  const int degree = 1;

  // YaspGrid section
  typedef YaspGrid<dim> GridType;
  typedef GridType::ctype DF;
  FieldVector<DF,dim> L = {1.0, 1.0};
  std::array<int,dim> N = {16, 16};
  std::bitset<dim> B(false);
  int overlap=1;
  GridType grid(L,N,B,overlap,MPIHelper::getCollectiveCommunication());
  grid.refineOptions(false); // keep overlap in cells
  grid.globalRefine(refinement);
  using GV = GridType::LeafGridView;
  GV gv=grid.leafGridView();

  // make user functions
  PoissonProblem<GV,double> problem;
  auto glambda = [&](const auto& e, const auto& x){return problem.g(e,x);};
  auto g = PDELab::makeGridFunctionFromCallable(gv,glambda);
  auto blambda = [&](const auto& i, const auto& x){return problem.bctype(i,x);};
  auto b = PDELab::makeBoundaryConditionFromCallable(gv,blambda);

  // Make grid function space
  typedef PDELab::OverlappingConformingDirichletConstraints CON; // NEW IN PARALLEL
  typedef PDELab::ISTL::VectorBackend<> VBE;

  using Basis = Functions::LagrangeBasis<GV,degree>;
  auto basis = std::make_shared<Basis>(gv);

  typedef PDELab::Experimental::GridFunctionSpace<Basis,
                                                  VBE,
                                                  CON> GridFunctionSpace;
  GridFunctionSpace gridFunctionSpace(basis);

  using FEM = GridFunctionSpace::Traits::FEM;

  gridFunctionSpace.name("Vh");

  // Assemble constraints
  typedef typename GridFunctionSpace::template
    ConstraintsContainer<double>::Type CC;
  CC cc;
  PDELab::constraints(b,gridFunctionSpace,cc); // assemble constraints

  // A coefficient vector
  using Z = PDELab::Backend::Vector<GridFunctionSpace,double>;
  Z z(gridFunctionSpace); // initial value

  // Make a grid function out of it
  typedef PDELab::DiscreteGridFunction<GridFunctionSpace,Z> ZDGF;
  ZDGF zdgf(gridFunctionSpace,z);

  // Fill the coefficient vector
  PDELab::interpolate(g,gridFunctionSpace,z);

  // Make a local operator
  typedef PDELab::ConvectionDiffusionFEM<decltype(problem),typename GridFunctionSpace::Traits::FiniteElementMap> LOP;
  LOP lop(problem);

  // Make a global operator
  typedef PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe((int)pow(1+2*degree,dim));
  typedef PDELab::GridOperator<
    GridFunctionSpace,  /* ansatz space */
    GridFunctionSpace,  /* test space */
    LOP,      /* local operator */
    MBE,      /* matrix backend */
    double,double,double, /* domain, range, jacobian field type*/
    CC,CC     /* constraints for ansatz and test space */
    > GO;
  GO go(gridFunctionSpace,cc,gridFunctionSpace,cc,lop,mbe);

  // Select a linear solver backend NEW IN PARALLEL
  typedef PDELab::ISTLBackend_CG_AMG_SSOR<GO> LS;
  int verbose = 0;
  LS ls(gridFunctionSpace,100,verbose);

  // solve nonlinear problem
  PDELab::Newton<GO,LS,Z> newton(go,z,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.setReduction(1e-10);
  newton.setMinLinearReduction(1e-4);
  newton.setMaxIterations(25);
  newton.setLineSearchMaxIterations(10);
  newton.apply();

  // Write VTK output file
  SubsamplingVTKWriter<GV> vtkwriter(gv,refinementIntervals(1));
  typedef PDELab::VTKGridFunctionAdapter<ZDGF> VTKF;
  vtkwriter.addVertexData(std::shared_ptr<VTKF>(new
                                         VTKF(zdgf,"fesol")));
  vtkwriter.write("pdelab-p-laplace-parallel-amd-result",
                  VTK::appendedraw);
}

int main(int argc, char** argv) try
{
  // Set up MPI if available
  MPIHelper::instance(argc, argv);

  // Test simple scalar spaces
  solvePoissonProblem<1>();
  solvePoissonProblem<2>();

  // Test with a parallel setup
  solveParallelPoissonProblem();

  return 0;
}
catch (Exception &e)
{
  std::cerr << "Dune reported error: " << e.what() << std::endl;
  return 1;
}
