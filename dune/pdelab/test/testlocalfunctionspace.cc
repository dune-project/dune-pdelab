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
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> P2GFS;
  P2GFS p2gfs(gv,q22dfem);

  // power grid function space
  typedef Dune::PDELab::PowerGridFunctionSpace<P2GFS,2> PowerGFS;
  PowerGFS powergfs(p2gfs);

  // make coefficent Vectors
  typedef typename P2GFS::template VectorContainer<double>::Type V;
  V x(p2gfs);
  x = 0.0;

  // make local function space object
  typename P2GFS::LocalFunctionSpace p2lfs(p2gfs);
  std::vector<double> xl(p2lfs.maxSize());

  // loop over elements
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  for (ElementIterator it = gv.template begin<0>();
	   it!=gv.template end<0>(); ++it)
	{
      p2lfs.bind(*it);
      p2lfs.debug();
      p2lfs.vread(x,xl);
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
