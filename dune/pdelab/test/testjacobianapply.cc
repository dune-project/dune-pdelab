// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune/pdelab.hh"

#include "nonlinearpoissonfem.hh"


template<typename Number>
class NonlinearPoissonProblem
{
  Number eta;
public:
  typedef Number value_type;

  //! Constructor without arg sets nonlinear term to zero
  NonlinearPoissonProblem () : eta(0.0) {}

  //! Constructor takes eta parameter
  NonlinearPoissonProblem (const Number& eta_) : eta(eta_) {}

  //! nonlinearity
  Number q (Number u) const
  {
    return eta*u*u;
  }

  //! derivative of nonlinearity
  Number qprime (Number u) const
  {
    return 2*eta*u;
  }

  //! right hand side
  template<typename E, typename X>
  Number f (const E& e, const X& x) const
  {
    auto global = e.geometry().global(x);
    return -2.0*x.size() + eta*global.two_norm2()*global.two_norm2();
  }

  //! boundary condition type function (true = Dirichlet)
  template<typename I, typename X>
  bool b (const I& i, const X& x) const
  {
    return true;
  }

  //! Dirichlet extension
  template<typename E, typename X>
  Number g (const E& e, const X& x) const
  {
    auto global = e.geometry().global(x);
    return global.two_norm2();
  }

  //! Neumann boundary condition
  template<typename I, typename X>
  Number j (const I& i, const X& x) const
  {
    return 0.0;
  }
};


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


bool run_linear_test()
{
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

  // Palpo
  using Range = typename GridOperator::Range;
  Range residual(gridOperator.testGridFunctionSpace(), 0.0);
  using Domain = typename GridOperator::Domain;
  Domain linearizationPoint(gridOperator.trialGridFunctionSpace(), 0.0);
  // gridOperator.nonlinear_jacobian_apply(coefficientVector, linearizationPoint, residual);
  // gridOperator.jacobian_apply(coefficientVector, linearizationPoint, residual);
  gridOperator.jacobian_apply(coefficientVector, residual);

  // Solver
  using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
  LinearSolver linearSolver(false);
  using Solver = Dune::PDELab::StationaryLinearProblemSolver<GridOperator, LinearSolver, CoefficientVector>;
  const double reduction = 1e-7;
  Solver solver(gridOperator, linearSolver, coefficientVector, reduction);
  solver.apply();

  // Visualization
  using VTKWriter = Dune::SubsamplingVTKWriter<GridView>;
  Dune::RefinementIntervals subint(2);
  VTKWriter vtkwriter(gridView, subint);
  std::string vtkfile("testjacobianapply_linear");
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gridFunctionSpace, coefficientVector,
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
  std::cout << "l2errorsquared: " << error << std::endl;

  // Let the test fail if the error is too large
  bool testfail(false);
  using std::abs;
  using std::isnan;
  if (isnan(error) or abs(error)>1e-7)
    testfail = true;
  return testfail;
}


bool run_nonlinear_test()
{
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

  // Solution vector
  using CoefficientVector = Dune::PDELab::Backend::Vector<GridFunctionSpace, DomainField>;
  CoefficientVector coefficientVector(gridFunctionSpace);

  // Discrete grid function of solution vetor
  using DiscreteGridFunction = Dune::PDELab::DiscreteGridFunction<GridFunctionSpace, CoefficientVector>;
  DiscreteGridFunction discreteGridFunction(gridFunctionSpace, coefficientVector);

  // Local operator (problem is nonlinear and problem class depends on the solution)
  using Problem = NonlinearPoissonProblem<RangeType>;
  Problem problem(2.0);
  using LocalOperator = NonlinearPoissonFEM<Problem, FiniteElementMap>;
  LocalOperator localOperator(problem);

  // Create constraints map
  using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<RangeType>::Type;
  ConstraintsContainer constraintsContainer;
  auto blambda = [&](const auto& i, const auto& x){return problem.b(i,x);};
  auto bctype = Dune::PDELab::makeBoundaryConditionFromCallable(gridView, blambda);
  Dune::PDELab::constraints(bctype, gridFunctionSpace, constraintsContainer);

  // Print number of DOFs
  std::cout << "gfs with " << gridFunctionSpace.size() << " dofs generated  "<< std::endl;
  std::cout << "cc with " << constraintsContainer.size() << " dofs generated  "<< std::endl;

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

  // Get grid function from Dirichlet boundary condition
  auto glambda = [&](const auto& e, const auto& x){return problem.g(e,x);};
  auto boundaryCondition = Dune::PDELab::makeGridFunctionFromCallable(gridView, glambda);
  Dune::PDELab::interpolate(boundaryCondition, gridFunctionSpace, coefficientVector);
  Dune::PDELab::set_nonconstrained_dofs(constraintsContainer, 0.0, coefficientVector);

  // Palpo
  using Range = typename GridOperator::Range;
  Range residual(gridOperator.testGridFunctionSpace(), 0.0);
  using Domain = typename GridOperator::Domain;
  Domain linearizationPoint(gridOperator.trialGridFunctionSpace(), 0.0);
  // gridOperator.nonlinear_jacobian_apply(coefficientVector, linearizationPoint, residual);
  gridOperator.jacobian_apply(coefficientVector, linearizationPoint, residual);
  // gridOperator.jacobian_apply(coefficientVector, residual);

  // Create Netwon solver
  using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
  LinearSolver linearSolver(false);
  const double reduction = 1e-7;
  using Solver = Dune::PDELab::Newton<GridOperator, LinearSolver, CoefficientVector>;
  Solver solver(gridOperator, coefficientVector, linearSolver);
  solver.apply(coefficientVector);

  // Visualization
  using VTKWriter = Dune::SubsamplingVTKWriter<GridView>;
  Dune::RefinementIntervals subint(2);
  VTKWriter vtkwriter(gridView, subint);
  std::string vtkfile("testjacobianapply_nonlinear");
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gridFunctionSpace, coefficientVector,
                                       Dune::PDELab::vtk::defaultNameScheme());
  vtkwriter.write(vtkfile, Dune::VTK::ascii);

  // Calculate error
  //
  // Note: The problem is set up in such a way that the Dirichlet boundary
  // condition is the exact solution of the problem.
  using DifferenceSquaredAdapter = Dune::PDELab::DifferenceSquaredAdapter<decltype(boundaryCondition),
                                                                          DiscreteGridFunction>;
  DifferenceSquaredAdapter differenceSquaredAdapder(boundaryCondition, discreteGridFunction);
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

int main(int argc, char** argv)
{
  try{
    // Maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    run_linear_test();
    run_nonlinear_test();

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
