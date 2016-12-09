// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Solve problem in parallel on non-overlapping grids using conforming linear finite elements
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<string>

#if GRID_UG
#include<dune/grid/uggrid.hh>
#else
#include<dune/grid/yaspgrid.hh>
#endif

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include<dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include "testnonoverlappingsinglephaseflow-problem.hh"

//===============================================================
// set up diffusion problem and solve it
//===============================================================

template< typename PROBLEM, typename GV, typename FEM>
void driver(PROBLEM& problem, GV gv, const FEM& fem,
            std::string filename)
{
  using ES = Dune::PDELab::NonOverlappingEntitySet<GV>;

  // constants and types and global variables
  typedef typename FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;
  Dune::Timer watch;

  // make function space
  //typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints<GV> CON;
  using CON = Dune::PDELab::ConformingDirichletConstraints;
  typedef Dune::PDELab::GridFunctionSpace
    <ES,FEM,CON,Dune::PDELab::ISTL::VectorBackend<> > GFS;
  CON con; // (gv);
  GFS gfs(gv,fem,con);
  gfs.name("solution");
  // con.compute_ghosts(gfs);

  // make constraints map and initialize it from a function and ghost
  typedef typename GFS::template ConstraintsContainer<R>::Type CC;
  CC cc;
  cc.clear();
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<PROBLEM> bctype(gv,problem);

  // make grid operator
  typedef Dune::PDELab::ConvectionDiffusionFEM<PROBLEM,FEM> LOP;
  LOP lop(problem);
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::GridOperator
      <GFS,GFS,LOP,
       MBE,
       R,R,R,CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // make coefficent Vector and initialize it from a function
  typedef typename GO::Traits::Domain V;
  V x(gfs,0.0);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PROBLEM> G;
  G g(gv,problem);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::constraints(bctype,gfs,cc,false);
  Dune::PDELab::set_nonconstrained_dofs(cc,0.0,x);

  typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<GO> LS;
  LS ls (go,5000,3,2);
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> SLP;
  SLP slp(go,ls,x,1e-12);
  slp.apply();

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x);
  vtkwriter.write(filename,Dune::VTK::ascii);
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
    {
      if(helper.rank()==0)
        std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
    }

    // make grid
    const int dim = 2;
#if GRID_UG
    typedef Dune::UGGrid<dim> GridType;
#else
    typedef Dune::YaspGrid<dim> GridType;
#endif

    // Build grid with uneven row/col number to provoke a reentrant corner in parallel case with UG
    Dune::FieldVector<typename GridType::ctype,dim> lowerLeft(0);
    Dune::FieldVector<typename GridType::ctype,dim> upperRight(1);
    std::array<unsigned int,dim> elements;
    std::fill(elements.begin(), elements.end(), 17);

    auto grid = Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
    grid->loadBalance();

    typedef GridType::LeafGridView GV;
    GV gv = grid->leafGridView();

    typedef Parameter<GV,double> PROBLEM;
    PROBLEM problem;

    typedef GridType::ctype DF;
    const int degree = 1;
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,degree> FEM;
    FEM fem(gv);

#if GRID_UG
    driver(problem,gv,fem,"nonoverlappingsinglephaseflow-ug");
#else
    driver(problem,gv,fem,"nonoverlappingsinglephaseflow-yasp");
#endif

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
