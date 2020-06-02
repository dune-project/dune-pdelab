// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune/pdelab.hh"

// Include local operator for nonlinear problem. This header is not part of
// pdelab.hh since it only exists for testing nonlinear problems
#include "dune/pdelab/test/localoperator/nonlinearconvectiondiffusiondg.hh"


//===============================================================
// These compile time definitions are set through CMakeLists.txt!
//===============================================================
// #define FASTDG
// #define SOR
// #define PARTIAL_MATRIX_FREE // This only has effect if SOR is not defined
// #define FULLY_MATRIX_FREE // This only has effect if SOR is not defined


// Poisson problem with nonlinearity where the Dirichlet boundary condition g
// is also the solution.
template <typename GridView, typename RangeType>
class NonlinearPoissonProblem : public Dune::PDELab::ConvectionDiffusionModelProblem<GridView, RangeType>
{
public:
  // Source term
  template<typename Element, typename Coord>
  auto f(const Element& element, const Coord& x) const
  {
    auto global = element.geometry().global(x);
    auto g = global.two_norm2();
    auto f = -4.0 + mu * g * g;
    return f;
  }

  // Nonlinearity
  RangeType q(const RangeType u) const
  {
    return mu*u*u;
  }

  // Derivative of Nonlinearity
  RangeType qprime(const RangeType u) const
  {
    return 2*mu*u;
  }

  // Boundary condition type
  template<typename Element, typename Coord>
  auto bctype(const Element& element, const Coord& x) const
  {
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  // Dirichlet extension
  template<typename Element, typename Coord>
  RangeType g (const Element& element, const Coord& x) const
  {
    auto global = element.geometry().global(x);
    auto g = global.two_norm2();
    return g;
  }

  RangeType mu = 2.0;
};


int main(int argc, char** argv)
{
  try{
     // Maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // Read initree
    Dune::ParameterTree initree;
    Dune::ParameterTreeParser::readINITree("matrix_free_nonlinear.ini", initree);

    // Create grid
    const int dim = 2;
    const int cells = initree.get<int>("grid.cells", 16);
    const int refine = initree.get<int>("grid.refine", 0);
    Dune::FieldVector<double,dim> lowerleft(0.0);
    Dune::FieldVector<double,dim> upperright(1.0);
    auto cellArray = Dune::filledArray<dim,unsigned int>(cells);
    using Grid = Dune::YaspGrid<dim>;
    auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerleft, upperright, cellArray);
    grid -> globalRefine(refine);
    using GridView = Grid::LeafGridView;
    GridView gridView = grid -> leafGridView();

    // Finite element map
    using DomainField = GridView::Grid::ctype;
    using RangeType = double;
    const int degree = 1;
    using FiniteElementMap = Dune::PDELab::QkDGLocalFiniteElementMap<DomainField, RangeType, degree, dim>;
    FiniteElementMap finiteElementMap;

    // Grid function space
    using Constraints = Dune::PDELab::NoConstraints;
    using VectorBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed, Dune::QkStuff::QkSize<degree, dim>::value>;
    using GridFunctionSpace = Dune::PDELab::GridFunctionSpace<GridView, FiniteElementMap, Constraints, VectorBackend>;
    GridFunctionSpace gridFunctionSpace(gridView, finiteElementMap);
    gridFunctionSpace.name("numerical_solution");

    // Local operator
    using Problem = NonlinearPoissonProblem<GridView, RangeType>;
    Problem problem;
    using LocalOperator = Dune::PDELab::NonlinearConvectionDiffusionDG<Problem, FiniteElementMap>;
    LocalOperator localOperator(problem,
                                Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,
                                Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,
                                3.0);

   // Create constraints map
    using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<RangeType>::Type;
    ConstraintsContainer constraintsContainer;
    constraintsContainer.clear();
    Dune::PDELab::constraints(gridFunctionSpace, constraintsContainer);

    // Grid operator
    gridFunctionSpace.update();
    using std::pow;
    const int dofestimate = pow(2, dim) * gridFunctionSpace.maxLocalSize();
    using MatrixBackend = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MatrixBackend matrixBackend(dofestimate);
#ifdef FASTDG
    std::cout << "Info: Using fast DG grid operator." << std::endl;
    using GridOperator = Dune::PDELab::FastDGGridOperator<GridFunctionSpace,
                                                          GridFunctionSpace,
                                                          LocalOperator,
                                                          MatrixBackend,
                                                          DomainField,
                                                          RangeType,
                                                          RangeType,
                                                          ConstraintsContainer,
                                                          ConstraintsContainer>;
#else
    using GridOperator = Dune::PDELab::GridOperator<GridFunctionSpace,
                                                    GridFunctionSpace,
                                                    LocalOperator,
                                                    MatrixBackend,
                                                    DomainField,
                                                    RangeType,
                                                    RangeType,
                                                    ConstraintsContainer,
                                                    ConstraintsContainer>;
#endif
    GridOperator gridOperator(gridFunctionSpace,
                              constraintsContainer,
                              gridFunctionSpace,
                              constraintsContainer,
                              localOperator,
                              matrixBackend);
    std::cout << "gfs with " << gridFunctionSpace.size() << " dofs generated  "<< std::endl;
    std::cout << "cc with " << constraintsContainer.size() << " dofs generated  "<< std::endl;

    // Solution vector
    using CoefficientVector = Dune::PDELab::Backend::Vector<GridFunctionSpace, DomainField>;
    CoefficientVector coefficientVector(gridFunctionSpace);
    coefficientVector = 0.0;
    using DirichletExtension = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>;
    DirichletExtension dirichletExtension(gridView, problem);

    // Solver settings
    const int maxiter = initree.get<int>("linear_solver.maxiter", 5000);
    const int verbosity = initree.get<int>("linear_solver.verbosity", 0);

    // Setup matrix-based solver
    std::cout << "Setup matrix-based solver" << std::endl;
#if defined(SOR)
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_BCGS_SOR;
    LinearSolver linearSolver(maxiter, verbosity);
    using Solver = Dune::PDELab::NewtonMethod<GridOperator, LinearSolver>;
#else
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_BCGS_Jac;
    LinearSolver linearSolver(maxiter, verbosity);
    using Solver = Dune::PDELab::NewtonMethod<GridOperator, LinearSolver>;
#endif

    //=========================
    // Setup matrix free solver
    //=========================
    std::cout << "Setup matrix-free solver" << std::endl;
#if defined(SOR)
    std::cout << "Info: Using matrix-free SOR." << std::endl;

    // Setup local operator for evaluating the block and point diagonal
    using BlockDiagonalLOP = Dune::PDELab::BlockDiagonalLocalOperatorWrapper<LocalOperator>;
    BlockDiagonalLOP blockDiagonalLOP(localOperator);
    using PointDiagonalLOP = Dune::PDELab::PointDiagonalLocalOperatorWrapper<LocalOperator>;
    PointDiagonalLOP pointDiagonalLOP(localOperator);
    using IBJPLOP = Dune::PDELab::IterativeBlockJacobiPreconditionerLocalOperator<BlockDiagonalLOP,
                                                                                  PointDiagonalLOP,
                                                                                  GridFunctionSpace,
                                                                                  DomainField,
                                                                                  Dune::BiCGSTABSolver>;
    Dune::PDELab::SolverStatistics<int> solverStat(gridView.comm());
    Dune::PDELab::BlockSolverOptions blockSolverOptions;
    blockSolverOptions._resreduction = initree.get<double>("block_inverse_solver.reduction", 1e-5);
    blockSolverOptions._maxiter = initree.get<int>("block_inverse_solver.maxiter", 100);
    blockSolverOptions._verbose = initree.get<int>("block_inverse_solver.verbose", 0);
    IBJPLOP ibjplop(blockDiagonalLOP,
                    pointDiagonalLOP,
                    gridFunctionSpace,
                    solverStat,
                    blockSolverOptions,
                    2);

    // Setup local operator for evaluating the block off-diagonals
    using BlockOffDiagonalLOP = Dune::PDELab::BlockOffDiagonalLocalOperatorWrapper<LocalOperator>;
    BlockOffDiagonalLOP blockOffDiagonalLOP(localOperator);

    // Setup local operator for SOR preconditioner
    using BlockSORPreconditionerLOP = Dune::PDELab::BlockSORPreconditionerLocalOperator<IBJPLOP,
                                                                                        BlockOffDiagonalLOP,
                                                                                        GridFunctionSpace>;
    BlockSORPreconditionerLOP blockSORPreconditionerLOP(ibjplop, blockOffDiagonalLOP, gridFunctionSpace, 1.0);

    // Setup grid operator for preconditioner application
#ifdef FASTDG
    using BlockSORPreconditionerGO = Dune::PDELab::FastDGGridOperator<GridFunctionSpace,
                                                                      GridFunctionSpace,
                                                                      BlockSORPreconditionerLOP,
                                                                      MatrixBackend,
                                                                      DomainField,
                                                                      RangeType,
                                                                      RangeType,
                                                                      ConstraintsContainer,
                                                                      ConstraintsContainer>;
#else
    using BlockSORPreconditionerGO = Dune::PDELab::GridOperator<GridFunctionSpace,
                                                                GridFunctionSpace,
                                                                BlockSORPreconditionerLOP,
                                                                MatrixBackend,
                                                                DomainField,
                                                                RangeType,
                                                                RangeType,
                                                                ConstraintsContainer,
                                                                ConstraintsContainer>;
#endif
    BlockSORPreconditionerGO blockSORPreconditionerGO(gridFunctionSpace,
                                                      constraintsContainer,
                                                      gridFunctionSpace,
                                                      constraintsContainer,
                                                      blockSORPreconditionerLOP,
                                                      matrixBackend);

    using LinearSolverMatrixFree = Dune::PDELab::ISTLBackend_SEQ_MatrixFree_Base<GridOperator,
                                                                                 BlockSORPreconditionerGO,
                                                                                 Dune::BiCGSTABSolver>;
    LinearSolverMatrixFree linearSolverMatrixFree(gridOperator, blockSORPreconditionerGO, maxiter, verbosity);
    using SolverMatrixFree = Dune::PDELab::MatrixFreeNewton<GridOperator, LinearSolverMatrixFree, CoefficientVector>;
#elif defined(PARTIAL_MATRIX_FREE) // #if defined(SOR)
    std::cout << "Info: Using partially matrix-free solver." << std::endl;

    using BlockDiagonalLOP = Dune::PDELab::BlockDiagonalLocalOperatorWrapper<LocalOperator>;
    BlockDiagonalLOP blockDiagonalLop(localOperator);
    using ABJPLOP = Dune::PDELab::AssembledBlockJacobiPreconditionerLocalOperator<BlockDiagonalLOP,
                                                                                  GridFunctionSpace,
                                                                                  DomainField>;
    ABJPLOP abjplop(blockDiagonalLop, gridFunctionSpace);
#ifdef FASTDG
    using ABJPLOP_GO = Dune::PDELab::FastDGGridOperator<GridFunctionSpace,
                                                        GridFunctionSpace,
                                                        ABJPLOP,
                                                        MatrixBackend,
                                                        DomainField,
                                                        RangeType,
                                                        RangeType,
                                                        ConstraintsContainer,
                                                        ConstraintsContainer>;
#else
    using ABJPLOP_GO = Dune::PDELab::GridOperator<GridFunctionSpace,
                                                  GridFunctionSpace,
                                                  ABJPLOP,
                                                  MatrixBackend,
                                                  DomainField,
                                                  RangeType,
                                                  RangeType,
                                                  ConstraintsContainer,
                                                  ConstraintsContainer>;
#endif
    ABJPLOP_GO abjplop_go(gridFunctionSpace,
                          constraintsContainer,
                          gridFunctionSpace,
                          constraintsContainer,
                          abjplop,
                          matrixBackend);
    using LinearSolverMatrixFree = Dune::PDELab::ISTLBackend_SEQ_MatrixFree_Base<GridOperator, ABJPLOP_GO, Dune::BiCGSTABSolver>;
    LinearSolverMatrixFree linearSolverMatrixFree(gridOperator, abjplop_go, maxiter, verbosity);
    using SolverMatrixFree = Dune::PDELab::MatrixFreeNewton<GridOperator, LinearSolverMatrixFree, CoefficientVector>;
#elif defined(FULLY_MATRIX_FREE)  // #if defined(SOR)
    std::cout << "Info: Using fully matrix-free solver." << std::endl;

    using BlockDiagonalLOP = Dune::PDELab::BlockDiagonalLocalOperatorWrapper<LocalOperator>;
    BlockDiagonalLOP blockDiagonalLop(localOperator);
    using PointDiagonalLOP = Dune::PDELab::PointDiagonalLocalOperatorWrapper<LocalOperator>;
    PointDiagonalLOP pointDiagonalLop(localOperator);
    using IBJPLOP = Dune::PDELab::IterativeBlockJacobiPreconditionerLocalOperator<BlockDiagonalLOP,
                                                                                  PointDiagonalLOP,
                                                                                  GridFunctionSpace,
                                                                                  DomainField,
                                                                                  Dune::BiCGSTABSolver>;
    Dune::PDELab::SolverStatistics<int> solverStat(gridView.comm());
    Dune::PDELab::BlockSolverOptions blockSolverOptions;
    blockSolverOptions._resreduction = initree.get<double>("block_inverse_solver.reduction", 1e-5);
    blockSolverOptions._maxiter = initree.get<int>("block_inverse_solver.maxiter", 100);
    blockSolverOptions._verbose = initree.get<int>("block_inverse_solver.verbose", 0);
    IBJPLOP ibjplop(blockDiagonalLop,
                    pointDiagonalLop,
                    gridFunctionSpace,
                    solverStat,
                    blockSolverOptions,
                    2);
#ifdef FASTDG
    using IBJPLOP_GO = Dune::PDELab::FastDGGridOperator<GridFunctionSpace,
                                                        GridFunctionSpace,
                                                        IBJPLOP,
                                                        MatrixBackend,
                                                        DomainField,
                                                        RangeType,
                                                        RangeType,
                                                        ConstraintsContainer,
                                                        ConstraintsContainer>;
#else
    using IBJPLOP_GO = Dune::PDELab::GridOperator<GridFunctionSpace,
                                                  GridFunctionSpace,
                                                  IBJPLOP,
                                                  MatrixBackend,
                                                  DomainField,
                                                  RangeType,
                                                  RangeType,
                                                  ConstraintsContainer,
                                                  ConstraintsContainer>;
#endif
    IBJPLOP_GO ibjplop_go(gridFunctionSpace,
                          constraintsContainer,
                          gridFunctionSpace,
                          constraintsContainer,
                          ibjplop,
                          matrixBackend);

    using LinearSolverMatrixFree = Dune::PDELab::ISTLBackend_SEQ_MatrixFree_Base<GridOperator, IBJPLOP_GO, Dune::BiCGSTABSolver>;
    LinearSolverMatrixFree linearSolverMatrixFree(gridOperator, ibjplop_go, maxiter, verbosity);
    using SolverMatrixFree = Dune::PDELab::MatrixFreeNewton<GridOperator, LinearSolverMatrixFree, CoefficientVector>;
#else
    static_assert(false);
    DUNE_THROW(Dune::Exception, "This should not happen");
#endif // #if defined(SOR)


    //==============
    // Solve the PDE
    //==============
    const double reduction = initree.get<double>("solver.reduction", 1e-12);
    Solver solver(gridOperator, linearSolver);
    std::cout << "Solve with matrix-based solver" << std::endl;
    solver.apply(coefficientVector);
    auto result = solver.result();
    coefficientVector = 0.0;
    SolverMatrixFree solverMatrixFree(gridOperator, linearSolverMatrixFree);
    std::cout << "Solve with matrix-free solver" << std::endl;
    solverMatrixFree.apply(coefficientVector);
    auto resultMatrixFree = solverMatrixFree.result();

    // // Visualization
    // //
    // // This is not necessary for testing but might be useful for debugging in
    // // case the test fails.
    // using VTKWriter = Dune::SubsamplingVTKWriter<GridView>;
    // Dune::RefinementIntervals subint(1);
    // VTKWriter vtkwriter(gridView, subint);
    // std::string vtkfile("matrix_free_nonlinear");
    // Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gridFunctionSpace, coefficientVector,
    //                                      Dune::PDELab::vtk::defaultNameScheme());
    // vtkwriter.write(vtkfile, Dune::VTK::ascii);

    // Calculate error
    //
    // Note: The problem is set up in such a way that the Dirichlet boundary
    // condition is the exact solution of the problem.
    using DiscreteGridFunction = Dune::PDELab::DiscreteGridFunction<GridFunctionSpace, CoefficientVector>;
    DiscreteGridFunction discreteGridFunction(gridFunctionSpace, coefficientVector);
    using DifferenceSquaredAdapter = Dune::PDELab::DifferenceSquaredAdapter<DirichletExtension, DiscreteGridFunction>;
    DifferenceSquaredAdapter differenceSquaredAdapder(dirichletExtension, discreteGridFunction);
    DifferenceSquaredAdapter::Traits::RangeType error(0.0);
    Dune::PDELab::integrateGridFunction(differenceSquaredAdapder, error, 10);
    std::cout << "l2errorsquared: " << error << std::endl;

    // Let the test fail if the error is too large
    bool testfail(false);
    using std::abs;
    using std::isnan;
    if (isnan(error) or abs(error)>1e-7)
      testfail = true;
    if (result.linear_solver_iterations != resultMatrixFree.linear_solver_iterations)
      testfail = true;
    return testfail;
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
