#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune/pdelab.hh"

template <class GridView, class RangeType>
class PoissonProblem
  : public Dune::PDELab::ConvectionDiffusionModelProblem<GridView, RangeType>
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

    // Create grid
    const int dim = 2;
    Dune::FieldVector<double,dim> lowerleft(0.0);
    Dune::FieldVector<double,dim> upperright(1.0);
    auto cells = Dune::filledArray<dim,unsigned int>(4);
    using Grid = Dune::YaspGrid<dim>;
    auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerleft, upperright, cells);
    grid -> globalRefine(3);
    using GridView = Grid::LeafGridView;
    GridView gridView = grid -> leafGridView();

    // Finite element map
    using DomainField = GridView::Grid::ctype;
    using RangeType = double;
    const int degree = 2;
    using FiniteElementMap = Dune::PDELab::QkLocalFiniteElementMap<GridView, DomainField, RangeType, degree>;
    FiniteElementMap finiteElementMap(gridView);

    // Grid function space
    using Constraints = Dune::PDELab::ConformingDirichletConstraints;
    using VectorBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
    using GridFunctionSpace = Dune::PDELab::GridFunctionSpace<GridView, FiniteElementMap, Constraints, VectorBackend>;
    GridFunctionSpace gridFunctionSpace(gridView, finiteElementMap);
    gridFunctionSpace.name("numerical_solution");

    // Local operator
    using Problem = PoissonProblem<GridView, RangeType>;
    Problem problem;
    using LocalOperator = Dune::PDELab::ConvectionDiffusionFEM<Problem,FiniteElementMap>;
    LocalOperator localOperator(problem);

    // Create constraints map
    using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<RangeType>::Type;
    ConstraintsContainer constraintsContainer;
    constraintsContainer.clear();
    Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> bctype(gridView, problem);
    Dune::PDELab::constraints(bctype, gridFunctionSpace, constraintsContainer);

    // Grid operator
    using MatrixBackend = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    const int dofestimate = 4 * gridFunctionSpace.maxLocalSize();
    MatrixBackend matrixBackend(dofestimate);
    using GridOperator = Dune::PDELab::GridOperator<GridFunctionSpace,
                                                    GridFunctionSpace,
                                                    LocalOperator,
                                                    MatrixBackend,
                                                    DomainField,
                                                    RangeType,
                                                    RangeType,
                                                    ConstraintsContainer,
                                                    ConstraintsContainer>;
    GridOperator gridOperator(gridFunctionSpace,
                              constraintsContainer,
                              gridFunctionSpace,
                              constraintsContainer,
                              localOperator,
                              matrixBackend);

    // Solution vector
    using CoefficientVector = Dune::PDELab::Backend::Vector<GridFunctionSpace, DomainField>;
    CoefficientVector coefficientVector(gridFunctionSpace);
    using DirichletExtension = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>;
    DirichletExtension dirichletExtension(gridView, problem);
    Dune::PDELab::interpolate(dirichletExtension, gridFunctionSpace, coefficientVector);
    Dune::PDELab::set_nonconstrained_dofs(constraintsContainer, 0.0, coefficientVector);

    // Copy for matrix free solution
    CoefficientVector coefficientVectorMatrixFree(coefficientVector);

    // Solver matrix based
    // using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
    // LinearSolver linearSolver(false);
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_BCGS_Richardson;
    LinearSolver linearSolver;
    using Solver = Dune::PDELab::StationaryLinearProblemSolver<GridOperator, LinearSolver, CoefficientVector>;
    const double reduction = 1e-10;
    Solver solver(gridOperator, linearSolver, coefficientVector, reduction);
    solver.apply();

    // Solve matrix free
    using LinearSolverMatrixFree = Dune::PDELab::ISTLBackend_SEQ_MatrixFree_BCGS_Richardson<GridOperator>;
    LinearSolverMatrixFree linearSolverMatrixFree(gridOperator);
    using TrialGridFunctionSpace = typename GridOperator::Traits::TrialGridFunctionSpace;
    using W = Dune::PDELab::Backend::Vector<TrialGridFunctionSpace,typename CoefficientVector::ElementType>;
    W residual(gridOperator.testGridFunctionSpace(), 0.0);
    gridOperator.residual(coefficientVectorMatrixFree, residual);
    CoefficientVector update(gridOperator.trialGridFunctionSpace(), 0.0);
    linearSolverMatrixFree.apply(update, residual, reduction);
    coefficientVectorMatrixFree -= update;

    // Visualization
    using VTKWriter = Dune::SubsamplingVTKWriter<GridView>;
    Dune::RefinementIntervals subint(2);
    VTKWriter vtkwriter(gridView, subint);
    std::string vtkfile("testmatrixfree");
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gridFunctionSpace, coefficientVector,
                                         Dune::PDELab::vtk::defaultNameScheme());
    gridFunctionSpace.name("numerical_solution_matrix_free");
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gridFunctionSpace, coefficientVectorMatrixFree,
                                         Dune::PDELab::vtk::defaultNameScheme());
    vtkwriter.write(vtkfile, Dune::VTK::ascii);

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
    std::cout << "l2errorsquared matrix based: " << error << std::endl;

    // Error for matrix free version
    DiscreteGridFunction discreteGridFunctionMatrixFree(gridFunctionSpace, coefficientVectorMatrixFree);
    DifferenceSquaredAdapter differenceSquaredAdapderMatrixFree(dirichletExtension, discreteGridFunctionMatrixFree);
    DifferenceSquaredAdapter::Traits::RangeType errorMatrixFree(0.0);
    Dune::PDELab::integrateGridFunction(differenceSquaredAdapderMatrixFree, errorMatrixFree, 10);
    std::cout << "l2errorsquared matrix free: " << errorMatrixFree << std::endl;

    // Let the test fail if the error is too large
    bool testfail(false);
    using std::abs;
    using std::isnan;
    if (isnan(error) or abs(error)>1e-7)
      testfail = true;
    if (isnan(error) or abs(errorMatrixFree)>1e-7)
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
