// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune/pdelab.hh"

// This should also work if we use the FastDGGridOperator
// #define FASTDG

//======================================
// These are set through CMakeLists.txt!
//======================================
// #define MATRIX_BASED
// #define PARTIAL_MATRIX_FREE
// #define FULLY_MATRIX_FREE
// #define MATRIX_FREE_SOR
// #define MATRIX_BASED_SOR


// Poisson problem where the Dirichlet boundary condition g is also the
// solution.
template <typename GridView, typename RangeType>
class PoissonProblem : public Dune::PDELab::ConvectionDiffusionModelProblem<GridView, RangeType>
{
public:
  // Source term
  template<typename Element, typename Coord>
  auto f(const Element& element, const Coord& x) const
  {
    auto global = element.geometry().global(x);
    auto c = (0.5-global[0])*(0.5-global[0]) + (0.5-global[1])*(0.5-global[1]);
    using std::exp;
    auto g = exp(-1.0*c);
    auto f = 4*(1.0-c)*g;
    return f;
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
    auto c = (0.5-global[0])*(0.5-global[0]) + (0.5-global[1])*(0.5-global[1]);
    using std::exp;
    auto g = exp(-1.0*c);
    return g;
  }
};


int main(int argc, char** argv)
{
  try{
     // Maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // Read initree
    Dune::ParameterTree initree;
    Dune::ParameterTreeParser::readINITree("matrix_free_linear.ini", initree);

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
    using Problem = PoissonProblem<GridView, RangeType>;
    Problem problem;
    using LocalOperator = Dune::PDELab::ConvectionDiffusionDG<Problem, FiniteElementMap>;
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
    //
    // Note: In the DG case we do not need to set the boundary values to the
    // Dirichlet boundary condition values. We create the DirichletExtension
    // anyway since we use it as exact solution in the error calculation. We
    // keep the code here to make it as similar to the non-DG case as possible.
    using CoefficientVector = Dune::PDELab::Backend::Vector<GridFunctionSpace, DomainField>;
    CoefficientVector coefficientVector(gridFunctionSpace);
    coefficientVector = 0.0;
    using DirichletExtension = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>;
    DirichletExtension dirichletExtension(gridView, problem);

    // Setup solver
    //
    // Here we choose between different linear solver backends. Since this is a
    // rather long block we marked the start and end ;)
    //
    //========================
    // {{{
    //========================
    const int maxiter = initree.get<int>("linear_solver.maxiter", 5000);
    const int verbosity = initree.get<int>("linear_solver.verbosity", 0);
#if defined(MATRIX_BASED)
    std::cout << "Info: Using matrix-based solver." << std::endl;
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_BCGS_Jac;
    LinearSolver linearSolver(maxiter, verbosity);
    using Solver = Dune::PDELab::StationaryLinearProblemSolver<GridOperator, LinearSolver, CoefficientVector>;
#elif defined(PARTIAL_MATRIX_FREE)
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
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_MatrixFree_Base<GridOperator, ABJPLOP_GO, Dune::BiCGSTABSolver>;
    LinearSolver linearSolver(gridOperator, abjplop_go, maxiter, verbosity);
    using Solver = Dune::PDELab::MatrixFreeStationaryLinearProblemSolver<GridOperator, LinearSolver, CoefficientVector>;
#elif defined(FULLY_MATRIX_FREE)
    std::cout << "Info: Using fully matrix-free solver." << std::endl;

    Dune::PDELab::SolverStatistics<int> solverStat(gridView.comm());
    Dune::PDELab::BlockSolverOptions blockSolverOptions;
    blockSolverOptions._resreduction = initree.get<double>("block_inverse_solver.reduction", 1e-5);
    blockSolverOptions._maxiter = initree.get<int>("block_inverse_solver.maxiter", 100);
    blockSolverOptions._verbose = initree.get<int>("block_inverse_solver.verbose", 0);

    // Setup local operator for evaluating the block and point diagonal
    using BlockDiagonalLOP = Dune::PDELab::BlockDiagonalLocalOperatorWrapper<LocalOperator>;
    BlockDiagonalLOP blockDiagonalLop(localOperator);
    using PointDiagonalLOP = Dune::PDELab::PointDiagonalLocalOperatorWrapper<LocalOperator>;
    PointDiagonalLOP pointDiagonalLop(localOperator);

    // Setup preconditioner local operator and grid operator
    using IBJPLOP = Dune::PDELab::IterativeBlockJacobiPreconditionerLocalOperator<BlockDiagonalLOP,
                                                                                  PointDiagonalLOP,
                                                                                  GridFunctionSpace,
                                                                                  DomainField,
                                                                                  Dune::BiCGSTABSolver>;
    IBJPLOP ibjplop(blockDiagonalLop, pointDiagonalLop, gridFunctionSpace, solverStat, blockSolverOptions, 2);
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

    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_MatrixFree_Base<GridOperator, IBJPLOP_GO, Dune::BiCGSTABSolver>;
    LinearSolver linearSolver(gridOperator, ibjplop_go, maxiter, verbosity);
    using Solver = Dune::PDELab::MatrixFreeStationaryLinearProblemSolver<GridOperator, LinearSolver, CoefficientVector>;
#elif defined(MATRIX_FREE_SOR)
    std::cout << "Info: Using matrix-free SOR." << std::endl;

    // Setup local operator for evaluating the block and point diagonal
    using BlockDiagonalLOP = Dune::PDELab::BlockDiagonalLocalOperatorWrapper<LocalOperator>;
    BlockDiagonalLOP blockDiagonalLop(localOperator);
    using PointDiagonalLOP = Dune::PDELab::PointDiagonalLocalOperatorWrapper<LocalOperator>;
    PointDiagonalLOP pointDiagonalLop(localOperator);

    // Setup local operator for jacobian preconditioner
    Dune::PDELab::SolverStatistics<int> solverStat(gridView.comm());
    Dune::PDELab::BlockSolverOptions blockSolverOptions;
    blockSolverOptions._resreduction = initree.get<double>("block_inverse_solver.reduction", 1e-5);
    blockSolverOptions._maxiter = initree.get<int>("block_inverse_solver.maxiter", 100);
    blockSolverOptions._verbose = initree.get<int>("block_inverse_solver.verbose", 0);
    using IBJPLOP = Dune::PDELab::IterativeBlockJacobiPreconditionerLocalOperator<BlockDiagonalLOP,
                                                                                  PointDiagonalLOP,
                                                                                  GridFunctionSpace,
                                                                                  DomainField,
                                                                                  Dune::BiCGSTABSolver>;
    IBJPLOP ibjplop(blockDiagonalLop, pointDiagonalLop, gridFunctionSpace, solverStat, blockSolverOptions, 2);

    // Setup local operator for evaluating the block off-diagonals
    using BlockOffDiagonalLOP = Dune::PDELab::BlockOffDiagonalLocalOperatorWrapper<LocalOperator>;
    BlockOffDiagonalLOP blockOffDiagonalLOP(localOperator);

    // Setup local operator for SOR preconditioner
    using BlockSORPreconditionerLOP = Dune::PDELab::BlockSORPreconditionerLocalOperator<IBJPLOP,
                                                                                        BlockOffDiagonalLOP,
                                                                                        GridFunctionSpace>;
    BlockSORPreconditionerLOP blockSORPreconditionerLOP(ibjplop, blockOffDiagonalLOP, gridFunctionSpace, 1.0);

    // Setup grid operator for preconditioner application
    using BlockSORPreconditionerGO = Dune::PDELab::GridOperator<GridFunctionSpace,
                                                                GridFunctionSpace,
                                                                BlockSORPreconditionerLOP,
                                                                MatrixBackend,
                                                                DomainField,
                                                                RangeType,
                                                                RangeType,
                                                                ConstraintsContainer,
                                                                ConstraintsContainer>;
    BlockSORPreconditionerGO blockSORPreconditionerGO(gridFunctionSpace,
                                                      constraintsContainer,
                                                      gridFunctionSpace,
                                                      constraintsContainer,
                                                      blockSORPreconditionerLOP,
                                                      matrixBackend);

    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_MatrixFree_Base<GridOperator,
                                                                       BlockSORPreconditionerGO,
                                                                       Dune::BiCGSTABSolver>;
    LinearSolver linearSolver(gridOperator, blockSORPreconditionerGO, maxiter, verbosity);
    using Solver = Dune::PDELab::MatrixFreeStationaryLinearProblemSolver<GridOperator, LinearSolver, CoefficientVector>;
#elif defined(MATRIX_BASED_SOR)
    std::cout << "Info: Using matrix-based SOR solver." << std::endl;
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_BCGS_SOR;
    LinearSolver linearSolver(maxiter, verbosity);
    using Solver = Dune::PDELab::StationaryLinearProblemSolver<GridOperator, LinearSolver, CoefficientVector>;
#else
    static_assert(false);
    DUNE_THROW(Dune::Exception, "This should not happen");
#endif
    //========================
    // }}}
    //========================


    // Solve the PDE
    const double reduction = initree.get<double>("solver.reduction", 1e-12);
    Solver solver(gridOperator, linearSolver, coefficientVector, reduction);
    solver.apply();

    // // Visualization
    // //
    // // This is not necessary for testing but might be useful for debugging in
    // // case the test fails.
    // using VTKWriter = Dune::SubsamplingVTKWriter<GridView>;
    // Dune::RefinementIntervals subint(1);
    // VTKWriter vtkwriter(gridView, subint);
    // std::string vtkfile("matrix_free_linear");
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
