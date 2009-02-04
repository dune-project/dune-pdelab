// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/grid/yaspgrid.hh>
#include"../finiteelementmap/q12dfem.hh"
#include"../finiteelementmap/q22dfem.hh"
#include"../finiteelementmap/edger12dfem.hh"
#include"../gridfunctionspace/gridfunctionspace.hh"

// test function trees
template<class GV> 
void testleafgridfunction (const GV& gv)
{
  // instantiate finite element maps
  typedef Dune::PDELab::Q12DLocalFiniteElementMap<float,double> Q12DFEM;
  Q12DFEM q12dfem;
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<float,double> Q22DFEM;
  Q22DFEM q22dfem;
  typedef Dune::PDELab::EdgeR12DLocalFiniteElementMap<float,double> EdgeR12DFEM;
  EdgeR12DFEM edger12dfem;
  
  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> GFS1; 
  GFS1 gfs1(gv,q12dfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> GFS2;
  GFS2 gfs2(gv,q22dfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM,Dune::PDELab::StdVectorBackend,
	Dune::PDELab::GridFunctionRestrictedMapper> GFS3;
  GFS3 gfs3(gv,q22dfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,EdgeR12DFEM> GFS4;
  GFS4 gfs4(gv,edger12dfem);

  // make coefficent Vectors
  typedef typename GFS1::template VectorContainer<double>::Type V1;
  V1 x1(gfs1);
  x1 = 0.0;
  typedef typename GFS2::template VectorContainer<double>::Type V2;
  V2 x2(gfs2);
  x2 = 0.0;
  typedef typename GFS4::template VectorContainer<double>::Type V4;
  V4 x4(gfs4);
  x4 = 0.0;
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
    grid.globalRefine(1);

	testleafgridfunction(grid.leafView());

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
