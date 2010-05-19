// -*- tab-width: 4; indent-tabs-mode: nil -*-

#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#include <string>
#include <sstream>

#include <dune/common/fvector.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#ifdef HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/gridfactory.hh>
#endif
#ifdef HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#ifdef HAVE_UG
#include <dune/grid/uggrid.hh>
#include <dune/grid/uggrid/uggridfactory.hh>
#endif

#include "../common/vtkexport.hh"
#include "../gridfunctionspace/gridfunctionspace.hh"
#include "../gridfunctionspace/gridfunctionspaceutilities.hh"
#include "../finiteelementmap/rt02dfem.hh"

template<typename GV>
void rt02DGridFunctionSpace (const GV& gv, const std::string &suffix = "")
{
  typedef typename GV::Grid::ctype D; // domain type
  typedef double R;                   // range type

  std::ostringstream filename;
  filename << "rt02dgridfunctionspace";
  if(suffix != "") filename << "-" << suffix;

  Dune::PDELab::RT02DLocalFiniteElementMap<GV,D,R> fem(gv);   // maps entity to finite element

  typedef Dune::PDELab::GridFunctionSpace<
    GV,
    Dune::PDELab::RT02DLocalFiniteElementMap<GV,D,R>
    > GFS;
  GFS gfs(gv,fem);                    // make grid function space

  typedef typename GFS::template VectorContainer<R>::Type X;
  X x(gfs,0.0);                       // make coefficient vector
  x[2] = 1.0;                         // set a component

  typedef Dune::PDELab::DiscreteGridFunctionPiola<GFS,X> DGF;
  DGF dgf(gfs,x);                     // make a grid function

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);  // plot result
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"rt02d"));
  vtkwriter.write(filename.str(),Dune::VTKOptions::ascii);
}


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // default exitcode 77 (=skipped); returned in case none of the supported
    // Grids were found
    int result = 77;

#ifdef HAVE_ALBERTA
    std::cout << "Alberta" << std::endl;
    {
      typedef Dune::AlbertaGrid<2, 2> Grid;
      Dune::GridFactory<Grid> gf;
      Dune::FieldVector<Grid::ctype, 2> pos;

      pos[0] = 0; pos[1] = 0; gf.insertVertex(pos);
      pos[0] = 1; pos[1] = 0; gf.insertVertex(pos);
      pos[0] = 0; pos[1] = 1; gf.insertVertex(pos);
      //pos[0] = 1; pos[1] = 1; gf.insertVertex(pos);

      Dune::GeometryType type;
      type.makeTriangle();
      std::vector<unsigned int> vid(3);

      vid[0] = 0; vid[1] = 1; vid[2] = 2; gf.insertElement(type, vid);
      //vid[0] = 1; vid[1] = 3; vid[2] = 2; gf.insertElement(type, vid);

      Grid *grid = gf.createGrid();
      //grid->globalRefine(1);

      rt02DGridFunctionSpace(grid->leafView(), "alberta");

      Dune::GridFactory<Grid>::destroyGrid(grid);
    }
    result = 0;
#endif // HAVE_ALBERTA


#ifdef HAVE_ALUGRID
    std::cout << "ALU" << std::endl;
    {
      typedef Dune::ALUSimplexGrid<2, 2> Grid;

      Grid grid("grids/2dtriangle.alu");
      //grid->globalRefine(1);

      rt02DGridFunctionSpace(grid.leafView(), "alu");
    }
    result = 0;
#endif // HAVE_ALUGRID

#ifdef HAVE_UG
    std::cout << "UG" << std::endl;
    {
      typedef Dune::UGGrid<2> Grid;
      Dune::GridFactory<Grid> gf;
      Dune::FieldVector<Grid::ctype, 2> pos;

      pos[0] = 0; pos[1] = 0; gf.insertVertex(pos);
      pos[0] = 1; pos[1] = 0; gf.insertVertex(pos);
      pos[0] = 0; pos[1] = 1; gf.insertVertex(pos);
      //pos[0] = 1; pos[1] = 1; gf.insertVertex(pos);

      Dune::GeometryType type;
      type.makeTriangle();
      std::vector<unsigned int> vid(3);

      vid[0] = 0; vid[1] = 1; vid[2] = 2; gf.insertElement(type, vid);
      //vid[0] = 1; vid[1] = 3; vid[2] = 2; gf.insertElement(type, vid);

      Grid *grid = gf.createGrid();
      //grid->globalRefine(1);

      rt02DGridFunctionSpace(grid->leafView(), "ug");

      delete grid;
    }
    result = 0;
#endif // HAVE_ALBERTA

    return result;
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
