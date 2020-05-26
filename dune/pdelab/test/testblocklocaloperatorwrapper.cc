// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune/pdelab.hh"


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

    // Create grid
    const int dim = 2;
    const int cells = 16;
    const int refine = 0;
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
    std::cout << "gfs with " << gridFunctionSpace.size() << " dofs generated  "<< std::endl;
    std::cout << "cc with " << constraintsContainer.size() << " dofs generated  "<< std::endl;

    // Coefficient vector
    using CoefficientVector = Dune::PDELab::Backend::Vector<GridFunctionSpace, DomainField>;
    CoefficientVector coefficientVector(gridFunctionSpace);
    coefficientVector = 0.0;
    using DirichletExtension = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>;
    DirichletExtension dirichletExtension(gridView, problem);

    // Assemble Jacobian
    using Jacobian = typename GridOperator::Traits::Jacobian;
    Jacobian jacobian(gridOperator);
    gridOperator.jacobian(coefficientVector, jacobian);
    using Dune::PDELab::Backend::native;
    // Dune::printmatrix(std::cout, native(jacobian), "global stiffness matrix","row", 9, 1);


    // Assemble block diagonal of Jacobian
    using BlockDiagonalLocalOperator = Dune::PDELab::BlockDiagonalLocalOperatorWrapper<LocalOperator>;
    BlockDiagonalLocalOperator blockDiagonalLocalOperator(localOperator);
    using BlockDiagonalGridOperator = Dune::PDELab::GridOperator<GridFunctionSpace,
                                                                 GridFunctionSpace,
                                                                 BlockDiagonalLocalOperator,
                                                                 MatrixBackend,
                                                                 DomainField,
                                                                 RangeType,
                                                                 RangeType,
                                                                 ConstraintsContainer,
                                                                 ConstraintsContainer>;
    BlockDiagonalGridOperator blockDiagonalGridOperator(gridFunctionSpace,
                                                        constraintsContainer,
                                                        gridFunctionSpace,
                                                        constraintsContainer,
                                                        blockDiagonalLocalOperator,
                                                        matrixBackend);
    typename BlockDiagonalGridOperator::Traits::Jacobian blockDiagonalJacobian(blockDiagonalGridOperator);
    blockDiagonalGridOperator.jacobian(coefficientVector, blockDiagonalJacobian);
    // Dune::printmatrix(std::cout, native(blockDiagonalJacobian),"block diagonal","row", 9, 1);


    // Assemble block off diagonal of Jacobian
    using BlockOffDiagonalLocalOperator = Dune::PDELab::BlockOffDiagonalLocalOperatorWrapper<LocalOperator>;
    BlockOffDiagonalLocalOperator blockOffDiagonalLocalOperator(localOperator);
    using BlockOffDiagonalGridOperator = Dune::PDELab::GridOperator<GridFunctionSpace,
                                                                    GridFunctionSpace,
                                                                    BlockOffDiagonalLocalOperator,
                                                                    MatrixBackend,
                                                                    DomainField,
                                                                    RangeType,
                                                                    RangeType,
                                                                    ConstraintsContainer,
                                                                    ConstraintsContainer>;
    BlockOffDiagonalGridOperator blockOffDiagonalGridOperator(gridFunctionSpace,
                                                              constraintsContainer,
                                                              gridFunctionSpace,
                                                              constraintsContainer,
                                                              blockOffDiagonalLocalOperator,
                                                              matrixBackend);
    typename BlockOffDiagonalGridOperator::Traits::Jacobian blockOffDiagonalJacobian(blockOffDiagonalGridOperator);
    blockOffDiagonalGridOperator.jacobian(coefficientVector, blockOffDiagonalJacobian);
    // Dune::printmatrix(std::cout, native(blockOffDiagonalJacobian),"block diagonal","row", 9, 1);


    bool testfail = false;

    // Test if the Block Diagonal Matrix Looks as expected:
    //
    // 1. The diagonal blocks are the same as the diagonal blocks of the jacobian
    // 2. Off-diagonal blocks do not exist
    for (std::size_t i=0; i<native(blockDiagonalJacobian).N(); ++i){
      for (std::size_t j=0; j<native(blockDiagonalJacobian).M(); ++j){
        if (i == j){
          auto block = native(blockDiagonalJacobian)[i][j] - native(jacobian)[i][j];
          if (block.infinity_norm() > 1e-14)
            testfail = true;
          // std::cout << block.infinity_norm() << std::endl;
        }
        else{
          if (native(blockDiagonalJacobian).exists(i, j))
            testfail = true;
          // std::cout << native(blockDiagonalJacobian).exists(i, j) << std::endl;
        }
      }
    }

    // Test if the Block Diagonal Matrix Looks as expected:
    //
    // 1. Diagonal blocks do not exists
    // 2. If an off diagonal block does not exist in the jacobian it does not exist
    // 3. If a diagonal block exists in the jacobian it is identical
    for (std::size_t i=0; i<native(blockOffDiagonalJacobian).N(); ++i){
      for (std::size_t j=0; j<native(blockOffDiagonalJacobian).M(); ++j){
        if (i == j){
          if (native(blockOffDiagonalJacobian).exists(i, j))
            testfail = true;
          // std::cout << native(blockOffDiagonalJacobian).exists(i, j) << std::endl;
        }
        else if (not native(jacobian).exists(i, j)){
          if (native(blockOffDiagonalJacobian).exists(i, j))
            testfail = true;
          // std::cout << native(blockOffDiagonalJacobian).exists(i, j) << std::endl;
        }
        else{
          auto block = native(blockOffDiagonalJacobian)[i][j] - native(jacobian)[i][j];
          if (block.infinity_norm() > 1e-14)
            testfail = true;
          // std::cout << block.infinity_norm() << std::endl;
        }
      }
    }

    if (testfail == true)
      std::cout << "This test failed" << std::endl;

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
