#include "config.h"

#include <iostream>
#include <vector>

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

  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  MatrixType stiffnessMatrix;

  // Construct Lagrangian finite element space basis
  using GridView = typename GridType::LeafGridView;
  auto gridView = grid.leafGridView();
  using Basis = Functions::PQkNodalBasis<GridView,order>;
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
  SeqILU0<MatrixType,VectorType,VectorType> ilu0(stiffnessMatrix,1.0);

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


int main(int argc, char** argv) try
{
  // Set up MPI if available
  MPIHelper::instance(argc, argv);

  // Test simple scalar spaces
  solvePoissonProblem<1>();
  solvePoissonProblem<2>();

  return 0;
}
catch (Exception &e)
{
  std::cerr << "Dune reported error: " << e.what() << std::endl;
  return 1;
}
