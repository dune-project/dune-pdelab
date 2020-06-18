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

    // Create Netwon solver
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
    LinearSolver linearSolver(false);
    const double reduction = 1e-7;
    using Solver = Dune::PDELab::NewtonMethod<GridOperator, LinearSolver>;
    Solver solver(gridOperator, linearSolver);

    // Set some parameters without loading parameter tree from ini file (just to show that it works)
    Dune::ParameterTree ptree;
    ptree["Verbosity"] = "4";
    ptree["Terminate.MaxIterations"] = "39";
    ptree["LineSearch.DampingFactor"] = "0.3";
    ptree["UseMaxNorm"] = "1";
    solver.setParameters(ptree);

    // Solve PDE
    solver.apply(coefficientVector);

    // Visualization
    using VTKWriter = Dune::SubsamplingVTKWriter<GridView>;
    Dune::RefinementIntervals subint(2);
    VTKWriter vtkwriter(gridView, subint);
    std::string vtkfile("testnewton");
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
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
  }
}
