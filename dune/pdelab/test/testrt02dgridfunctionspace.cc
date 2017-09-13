// -*- tab-width: 4; indent-tabs-mode: nil -*-

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/gridfactory.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#include <dune/grid/uggrid/uggridfactory.hh>
#endif

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/finiteelementmap/raviartthomasfem.hh>

template<typename GV>
void rt02DGridFunctionSpace (const GV& gv, const std::string &suffix = "")
{
  typedef typename GV::Grid::ctype D; // domain type
  typedef double R;                   // range type

  std::ostringstream filename;
  filename << "rt02dgridfunctionspace";
  if(suffix != "") filename << "-" << suffix;

  Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV,D,R,0,Dune::GeometryType::simplex> fem(gv);   // maps entity to finite element

  typedef Dune::PDELab::GridFunctionSpace<
    GV,
    Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV,D,R,0,Dune::GeometryType::simplex>
    > GFS;
  GFS gfs(gv,fem);                    // make grid function space

  using X = Dune::PDELab::Backend::Vector<GFS, R>;
  X x(gfs,0.0);                       // make coefficient vector
  Dune::PDELab::Backend::native(x)[2] = 1.0;                         // set a component

  typedef Dune::PDELab::DiscreteGridFunctionPiola<GFS,X> DGF;
  DGF dgf(gfs,x);                     // make a grid function

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);  // plot result
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgf,"rt02d"));
  vtkwriter.write(filename.str(),Dune::VTK::ascii);
}


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // default exitcode 77 (=skipped); returned in case none of the supported
    // Grids were found
    int result = 77;

#if HAVE_ALBERTA
    std::cout << "Alberta" << std::endl;
    {
      typedef Dune::AlbertaGrid<2, 2> Grid;
      Dune::GridFactory<Grid> gf;
      Dune::FieldVector<Grid::ctype, 2> pos;

      pos[0] = 0; pos[1] = 0; gf.insertVertex(pos);
      pos[0] = 1; pos[1] = 0; gf.insertVertex(pos);
      pos[0] = 0; pos[1] = 1; gf.insertVertex(pos);
      //pos[0] = 1; pos[1] = 1; gf.insertVertex(pos);

      auto type = Dune::GeometryTypes::triangle;
      std::vector<unsigned int> vid(3);

      vid[0] = 0; vid[1] = 1; vid[2] = 2; gf.insertElement(type, vid);
      //vid[0] = 1; vid[1] = 3; vid[2] = 2; gf.insertElement(type, vid);

      Grid *grid = gf.createGrid();
      //grid->globalRefine(1);

      rt02DGridFunctionSpace(grid->leafGridView(), "alberta");

      Dune::GridFactory<Grid>::destroyGrid(grid);
    }
    result = 0;
#endif // HAVE_ALBERTA


#if HAVE_DUNE_ALUGRID
    std::cout << "ALU" << std::endl;
    {
      using ALUType = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>;
      auto alugrid = Dune::StructuredGridFactory<ALUType>::createSimplexGrid(Dune::FieldVector<ALUType::ctype, 2>(0.0), Dune::FieldVector<ALUType::ctype, 2>(1.0), Dune::make_array(1u, 1u));
      alugrid->globalRefine(4);

      rt02DGridFunctionSpace(alugrid->leafGridView(), "alu");
    }
    result = 0;
#endif // HAVE_DUNE_ALUGRID

#if HAVE_UG
    std::cout << "UG" << std::endl;
    {
      typedef Dune::UGGrid<2> Grid;
      Dune::GridFactory<Grid> gf;
      Dune::FieldVector<Grid::ctype, 2> pos;

      pos[0] = 0; pos[1] = 0; gf.insertVertex(pos);
      pos[0] = 1; pos[1] = 0; gf.insertVertex(pos);
      pos[0] = 0; pos[1] = 1; gf.insertVertex(pos);
      //pos[0] = 1; pos[1] = 1; gf.insertVertex(pos);

      auto type = Dune::GeometryTypes::triangle;
      std::vector<unsigned int> vid(3);

      vid[0] = 0; vid[1] = 1; vid[2] = 2; gf.insertElement(type, vid);
      //vid[0] = 1; vid[1] = 3; vid[2] = 2; gf.insertElement(type, vid);

      Grid *grid = gf.createGrid();
      //grid->globalRefine(1);

      rt02DGridFunctionSpace(grid->leafGridView(), "ug");

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
