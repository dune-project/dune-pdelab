// -*- tab-width: 4; indent-tabs-mode: nil -*-

#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#include <dune/common/fvector.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/albertagrid/agrid.hh>
#include <dune/grid/albertagrid/gridfactory.hh>

#include "../common/vtkexport.hh"
#include "../gridfunctionspace/gridfunctionspace.hh"
#include "../gridfunctionspace/gridfunctionspaceutilities.hh"
#include "../finiteelementmap/edges02dfem.hh"

template<typename GV>
void edgeS02DGridFunctionSpace (const GV& gv)
{
  typedef typename GV::Grid::ctype D; // domain type
  typedef double R;                   // range type

  Dune::PDELab::EdgeS02DLocalFiniteElementMap<GV,D,R> fem(gv);   // maps entity to finite element

  typedef Dune::PDELab::GridFunctionSpace<
    GV,
    Dune::PDELab::EdgeS02DLocalFiniteElementMap<GV,D,R>
    > GFS;
  GFS gfs(gv,fem);                    // make grid function space

  typedef typename GFS::template VectorContainer<R>::Type X;
  X x(gfs,0.0);                       // make coefficient vector
  x[2] = 1.0;                         // set a component

  typedef Dune::PDELab::DiscreteGridFunctionPiola<GFS,X> DGF;
  DGF dgf(gfs,x);                     // make a grid function

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);  // plot result
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"edges02d"));
  vtkwriter.write("edges02dgridfunctionspace",Dune::VTKOptions::ascii);
}


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    // make grid
    typedef Dune::AlbertaGrid<2, 2> Grid;
    Dune::GridFactory<Grid> gf;
    {
      Dune::FieldVector<Grid::ctype, 2> pos;

      pos[0] = 0; pos[1] = 0; gf.insertVertex(pos);
      pos[0] = 1; pos[1] = 0; gf.insertVertex(pos);
      pos[0] = 0; pos[1] = 1; gf.insertVertex(pos);
      //pos[0] = 1; pos[1] = 1; gf.insertVertex(pos);

      Dune::GeometryType type;
      type.makeTriangle();
      std::vector<unsigned int> vid(3);

      vid[0] = 2; vid[1] = 1; vid[2] = 0; gf.insertElement(type, vid);
      //vid[0] = 1; vid[1] = 3; vid[2] = 2; gf.insertElement(type, vid);
    }

    Grid *grid = gf.createGrid("AlbertaGrid", true);
    //grid->globalRefine(1);

    edgeS02DGridFunctionSpace(grid->leafView());

    Dune::GridFactory<Grid>::destroyGrid(grid);

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
