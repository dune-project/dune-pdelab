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
#include"../gridfunctionspace/gridfunctionspace.hh"

// test function trees
template<class GV> 
void test (const GV& gv)
{
  // instantiate finite element maps
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<float,double> Q22DFEM;
  Q22DFEM q22dfem;
  
  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> GFS;
  GFS gfs(gv,q22dfem);

  // make coefficent Vectors
  typedef typename GFS::template VectorContainer<double>::Type V;
  V x(gfs);
  x = 0.0;

  // make local function space object
  typename GFS::LocalFunctionSpace lfs(gfs);
  std::vector<double> xl(lfs.maxSize());

  // loop over elements
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  for (ElementIterator it = gv.template begin<0>();
	   it!=gv.template end<0>(); ++it)
	{
      lfs.bind(*it);
      lfs.debug();
      lfs.vread(x,xl);
	}
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

	test(grid.leafView());

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
