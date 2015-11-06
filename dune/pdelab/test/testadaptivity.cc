// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file

    \brief Solve elliptic problem with adaptive conforming finite element method
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/typetraits.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/uggrid.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include<dune/pdelab/adaptivity/adaptivity.hh>

#include<dune/pdelab/test/testadaptivity-bctype.hh>
#include<dune/pdelab/test/testadaptivity-bcextension.hh>
#include<dune/pdelab/test/testadaptivity-operator.hh>
#include<dune/pdelab/test/testadaptivity-error-indicator.hh>
#include<dune/pdelab/test/testadaptivity-adaptivity.hh>

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    int startLevel = 4;
    int maxLevel = 6;

    // make grid
    const int dim = 2;
    typedef Dune::UGGrid<dim> Grid;
    Dune::FieldVector<Grid::ctype, dim> ll(0.0);
    Dune::FieldVector<Grid::ctype, dim> ur(1.0);
    std::array<unsigned int, dim> elements;
    std::fill(elements.begin(), elements.end(), 1);

    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, elements);
    grid->globalRefine(startLevel);

    typedef Grid::LeafGridView GV;
    GV gv = grid->leafGridView();
    adaptivity(*grid,gv,startLevel,maxLevel);
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
