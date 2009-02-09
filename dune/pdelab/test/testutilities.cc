// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<vector>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/grid/yaspgrid.hh>
#include"../finiteelementmap/q22dfem.hh"
#include"../finiteelementmap/q12dfem.hh"
#include"../finiteelementmap/edger12dfem.hh"
#include"../gridfunctionspace/gridfunctionspace.hh"
#include"../gridfunctionspace/gridfunctionspaceutilities.hh"
#include"../common/vtkexport.hh"

// generate a Q1 function and output it
template<class GV> 
void testq1 (const GV& gv)
{
  // instantiate finite element maps
  typedef Dune::PDELab::Q12DLocalFiniteElementMap<typename GV::Grid::ctype,double> Q12DFEM;
  Q12DFEM q12dfem;
  
  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> Q1GFS;
  Q1GFS q1gfs(gv,q12dfem);

  // make coefficent Vectors
  typedef typename Q1GFS::template VectorContainer<double>::Type V;
  V x(q1gfs);
  x = 0.0;
  x[3] = 1.0;

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<Q1GFS,V> DGF;
  DGF dgf(q1gfs,x);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"test"));
  vtkwriter.write("q1",Dune::VTKOptions::ascii);
}

// generate an edge element function and output it
template<class GV> 
void testedger (const GV& gv)
{
  // instantiate finite element maps
  typedef Dune::PDELab::EdgeR12DLocalFiniteElementMap<typename GV::Grid::ctype,double> EdgeR12DFEM;
  EdgeR12DFEM edger12dfem;
  
  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,EdgeR12DFEM> EGFS;
  EGFS egfs(gv,edger12dfem);

  // make coefficent Vectors
  typedef typename EGFS::template VectorContainer<double>::Type V;
  V x(egfs);
  x = 1.0;

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<EGFS,V> DGF;
  DGF dgf(egfs,x);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"test"));
  vtkwriter.write("edger",Dune::VTKOptions::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

	// need a grid in order to test grid functions
	Dune::FieldVector<double,2> L(1.0);
	Dune::FieldVector<int,2> N(1);
	Dune::FieldVector<bool,2> B(false);
	Dune::YaspGrid<2,2> grid(L,N,B,0);
    grid.globalRefine(3);

	testq1(grid.leafView());
	testedger(grid.leafView());

	// test passed
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
