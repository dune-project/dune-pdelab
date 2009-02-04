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
#include"../gridfunctionspace/leafgridfunctionspace.hh"

// test function trees
template<class GV> 
void testleafgridfunction (const GV& gv)
{
  // instantiate finite element maps
  typedef Dune::PDELab::Q12DLocalFiniteElementMap<float,double> Q12DFEM;
  Q12DFEM q12dfem;
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<float,double> Q22DFEM;
  Q22DFEM q22dfem;
  
  // make a grid function space
  Dune::PDELab::LeafGridFunctionSpace<GV,Q12DFEM> gfs1(gv,q12dfem);
  Dune::PDELab::LeafGridFunctionSpace<GV,Q22DFEM> gfs2(gv,q22dfem);
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
