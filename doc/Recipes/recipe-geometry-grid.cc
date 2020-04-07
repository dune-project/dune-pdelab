#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/geometrygrid/grid.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/pdelab.hh>

/**
 * \page recipe-geometry-grid Transforming a cartesian mesh
 *
 * Using GeometryGrid, one can transform an existing mesh into a different shape.
 * This is particularly useful for modelling complex geometries based on efficient
 * structured grids like YaspGrid.
 * 
 * A grid transformation is a user-defined function which gives the transformation
 * of a point x in the base cartesian grid to a final geometry:
 *
 * \f$ y = F(x) : \Omega -> \hat \Omega. \f$
 *
 * First, let's set up a rectangular grid:
 * \snippet recipe-geometry-grid.cc Setting up grid
 *
 * Then, lets define a function that maps the grid to a new geometry:
 * \snippet recipe-geometry-grid.cc Define function
 *
 * Finally, lets map our initial grid using the GridTransformation we defined:
 * \snippet recipe-geometry-grid.cc Mapping grid
 *
 * Full example code: @ref recipe-geometry-grid.cc
 * \example recipe-geometry-grid.cc
 * See explanation at @ref recipe-geometry-grid
 */

template<typename GV, typename RF>
class GenericEllipticProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1 : 0;
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    return e.geometry().global(xlocal)[0] < 0.7 ? 5.0 : 1.0;
  }

  //! boundary condition type function
  /* return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet for Dirichlet boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann for flux boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow for outflow boundary conditions
   */
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    return 0.0;
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    return 0.0;
  }
};

// [Define function]
template <int dim>
class GridTransformation
: public Dune :: AnalyticalCoordFunction< double, dim, dim, GridTransformation <dim> >{
    typedef GridTransformation This;
    typedef Dune :: AnalyticalCoordFunction< double, dim, dim, This > Base;

    public:
    typedef typename Base :: DomainVector DomainVector;
    typedef typename Base :: RangeVector RangeVector;

    GridTransformation(){}

    void evaluate(const DomainVector &x, RangeVector &y) const{
        y = x;
        if(x[0] < 0.8)
            y[1] = (1.0 + 5.0/4.0 * (sin(M_PI/18.0) - 1.0) * x[0]) * (x[1] - 1.0);
        else
            y[1] = sin((x[0] - 0.6)/3.6 * M_PI) * (x[1] - 1.0);
        if(x[0] > 3.8)
            y[0] += 0.5*(x[0] - 3.8) * (1.0 - pow(x[1] - 1.0, 2.0));
    }
};
//! [Define function]

int main(int argc, char** argv)
{
  try{
    // Initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // [Setting up grid]
    const unsigned int dim = 2;
    Dune::FieldVector<double,dim> L = {4.0,2.0};
    std::array<int,dim> N ={64,32};

    typedef Dune::YaspGrid<dim> SquareGrid;
    SquareGrid sgrid(L,N);
    //! [Setting up grid]

    // [Mapping grid]
    typedef GridTransformation<dim> GridTransformation;
    GridTransformation gTrafo;

    typedef typename Dune::GeometryGrid<SquareGrid,GridTransformation> Grid;
    Grid grid(sgrid,gTrafo);
    //! [Mapping grid]

    // define parameters
    typedef double NumberType;

    // need a grid in order to test grid functions
    constexpr unsigned int degree = 1;
    constexpr std::size_t nonzeros = Dune::power(2*degree+1,dim);

    // make problem parameters
    typedef GenericEllipticProblem<typename Grid::LeafGridView,NumberType> Problem;
    Problem problem;
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
    BCType bctype(grid.leafGridView(),problem);

    // Make FEM space
    typedef typename Grid::ctype DF;
    typedef Dune::PDELab::QkLocalFiniteElementMap<Grid::LeafGridView,DF,NumberType,1> FEM;
    FEM fem(grid.leafGridView());

    // make function space with constraints
    typedef Dune::PDELab::GridFunctionSpace<Grid::LeafGridView,FEM,
    Dune::PDELab::ConformingDirichletConstraints,
    Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,1>> GFS;
    GFS gfs(grid.leafGridView(),fem); gfs.name("solution");

    // Make constraints
    typedef typename GFS::template ConstraintsContainer<NumberType>::Type CC;
    CC cc;
    Dune::PDELab::constraints(bctype,gfs,cc);

    // Set up local operator
    typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
    LOP lop(problem);

    // Set up grid operator
    typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
    typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC> GO;
    auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros));

    // Define solution
    typedef Dune::PDELab::Backend::Vector<GFS,NumberType> X;
    X x(gfs,0.0);
    typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
    G g(grid.leafGridView(),problem);
    Dune::PDELab::interpolate(g,gfs,x);

    // Solve Poisson equation using AMG
    typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
    LS ls(100,3);
    typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,X> SLP;
    SLP slp(go,ls,x,1e-10);
    slp.apply(); // here all the work is done!

    // Plot the mesh
    Dune::SubsamplingVTKWriter<decltype(grid.leafGridView())> vtkwriter(grid.leafGridView(),Dune::refinementLevels(0));
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x);
    vtkwriter.write("mesh",Dune::VTK::appendedraw);

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
