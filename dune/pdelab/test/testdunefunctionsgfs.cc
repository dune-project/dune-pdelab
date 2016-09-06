#include "config.h"

#include <iostream>
#include <vector>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>

#include <dune/pdelab/gridfunctionspace/dunefunctionsgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/dunefunctionslfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/localoperator/linearelasticity.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

#include<dune/pdelab/finiteelementmap/qkfem.hh>

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

  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  MatrixType stiffnessMatrix;
  VectorType rhs;

  // Construct Lagrangian finite element space basis
  using GridView_ = typename GridType::LeafGridView;
  auto gridView_ = grid.leafGridView();
  using GridView = Dune::PDELab::AllEntitySet<GridView_>;
  auto gridView = GridView(gridView_,GridView::allCodims());
  using Basis = Functions::PQkNodalBasis<GridView,order>;
  auto basis = std::make_shared<Basis>(gridView);

  // What precisely does this do?
  typedef PDELab::ConformingDirichletConstraints Constraints;
  auto con = std::make_shared<Constraints>();

  typedef PDELab::Experimental::GridFunctionSpace<Basis,VectorType,Constraints> GridFunctionSpace;
  GridFunctionSpace gfs(basis,con);

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
                               PDELab::istl::BCRSMatrixBackend<>,
                               double,double,double,C,C> GO;
  GO go(gfs,constraintsContainer,gfs,constraintsContainer,lop, {9});

  // Pseudo "current" coefficient vector, not actually used
  typedef typename GO::Traits::Domain VectorContainer;
  VectorContainer x0(gfs);
  x0 = 0.0;              // set all entries to zero

  // represent operator as a matrix
  typedef typename GO::Jacobian MatrixContainer;
  MatrixContainer m(go,stiffnessMatrix);    // Use the stiffnessMatrix object for the actual storage
  m = 0.0;             // set all matrix entries to zero

  go.jacobian(x0,m);      // assemble the stiffness matrix

  // evaluate residual w.r.t pseudo "current" iterate
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
  SeqILU0<MatrixType,VectorType,VectorType> ilu0(stiffnessMatrix,1.0);

  // A preconditioned conjugate-gradient solver
  CGSolver<VectorType> cg(op,ilu0,1E-4,
                          50,   // maximum number of iterations
                          2);

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  cg.apply(x, rhs, statistics);

  // Output result to VTK file
  auto pressureFunction = Functions::makeDiscreteGlobalBasisFunction<double>(*basis,x);

  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
  vtkWriter.addVertexData(pressureFunction, VTK::FieldInfo("pressure", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("testdunefunctionsgfs-poisson");

}

template <class GridView>
class StVenantKirchhoffParameters
  : public Dune::PDELab::LinearElasticityParameterInterface<
  Dune::PDELab::LinearElasticityParameterTraits<GridView, double>,
  StVenantKirchhoffParameters<GridView> >
{
public:
  typedef Dune::PDELab::LinearElasticityParameterTraits<GridView, double> Traits;

  StVenantKirchhoffParameters(typename Traits::RangeFieldType l,
                              typename Traits::RangeFieldType m) :
    lambda_(l), mu_(m)
  {}

  void f (const typename Traits::ElementType& e,
          const typename Traits::DomainType& x,
          typename Traits::RangeType & y) const
  {
    std::fill(y.begin(), y.end(), 1e5);
  }

  template<typename I>
  bool isDirichlet(const I & ig,
                   const typename Traits::IntersectionDomainType & coord) const
  {
    return true;
  }

  typename Traits::RangeFieldType
  lambda (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return lambda_;
  }

  typename Traits::RangeFieldType
  mu (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return mu_;
  }

private:
  typename Traits::RangeFieldType lambda_;
  typename Traits::RangeFieldType mu_;
};


template <int order>
void solveElasticityProblem()
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

  // Construct Lagrangian finite element space basis
  using GridView = typename GridType::LeafGridView;
  auto gridView = grid.leafGridView();

  typedef BlockVector<FieldVector<double,dim> > VectorType;

  using namespace Functions::BasisBuilder;

  auto basis = makeBasis(
    gridView,
    power<dim>(
      lagrange<order>(),
      flatInterleaved())
    );

  using Basis = decltype(basis);

  typedef PDELab::ConformingDirichletConstraints Constraints;
  Constraints con;

  typedef PDELab::Experimental::GridFunctionSpace<Basis,VectorType,Constraints> GFS;
  GFS gfs(basis,con);

  gfs.name("displacement");

  // Container for the Dirichlet boundary conditions
  typedef typename GFS::template ConstraintsContainer<double>::Type C;
  C constraintsContainer;

  // create the model describing our problem
  auto lambda = 1e6;
  auto mu     = 1e6;
  StVenantKirchhoffParameters<GridView> model(lambda, mu);

  PDELab::constraints(model,gfs,constraintsContainer);

  // make grid operator
  PDELab::LinearElasticity<StVenantKirchhoffParameters<GridView> > lop(model);

  // set up linear operator acting on the FEM space
  typedef PDELab::GridOperator<GFS,
                               GFS,
                               PDELab::LinearElasticity<StVenantKirchhoffParameters<GridView> >,
                               PDELab::istl::BCRSMatrixBackend<>,
                               double,double,double,C,C> GridOperator;

  GridOperator gridOperator(gfs, constraintsContainer, gfs, constraintsContainer, lop, {9});

  typedef typename GridOperator::Traits::Domain V;
  typedef typename GridOperator::Jacobian M;
  using MatrixType = typename M::Container;   //  BCRSMatrix<FieldMatrix<double, dim, dim> >
  //using VectorType = typename V::Container;   //  BlockVector<FieldVector<double,dim> >

  // Dummy coefficient vector
  V x0(gfs);
  x0 = 0.0;

  // represent operator as a matrix
  MatrixType stiffnessMatrix;
  M m(gridOperator,stiffnessMatrix);
  m = 0.0;

  // Compute stiffness matrix
  gridOperator.jacobian(x0, m);

  // evaluate residual w.r.t pseudo "current" iterate
  VectorType rhs;
  V r(gfs,rhs);  //Use the rhs object for the actual storage
  r = 0.0;

  gridOperator.residual(x0,r);

  ///////////////////////////
  //   Compute solution
  ///////////////////////////

  VectorType x(rhs.size());
  x = 0;

  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

  // A preconditioner
  SeqILU0<MatrixType,VectorType,VectorType> ilu0(stiffnessMatrix,1.0);

  // A preconditioned conjugate-gradient solver
  CGSolver<VectorType> cg(op,ilu0,1E-4,
                          50,   // maximum number of iterations
                          2);

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  cg.apply(x, rhs, statistics);

  // Output result to VTK file
  auto pressureFunction = Functions::makeDiscreteGlobalBasisFunction<FieldVector<double,dim> >(basis,x);

  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
  vtkWriter.addVertexData(pressureFunction, VTK::FieldInfo("displacement", VTK::FieldInfo::Type::vector, dim));
  vtkWriter.write("testdunefunctionsgfs-elasticity");
}


int main(int argc, char** argv) try
{
  // Set up MPI if available
  MPIHelper::instance(argc, argv);

  // Test simple scalar spaces
  solvePoissonProblem<1>();
  solvePoissonProblem<2>();

  // Test a vector-valued space
  solveElasticityProblem<1>();
  solveElasticityProblem<2>();

  return 0;
}
catch (Exception &e)
{
  std::cerr << "Dune reported error: " << e.what() << std::endl;
  return 1;
}
