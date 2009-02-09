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
#include"../gridfunctionspace/gridfunctionspace.hh"

// test function trees
template<class GV> 
void test (const GV& gv)
{
  // instantiate finite element maps
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<float,double> Q22DFEM;
  Q22DFEM q22dfem;
  typedef Dune::PDELab::Q12DLocalFiniteElementMap<float,double> Q12DFEM;
  Q12DFEM q12dfem;
  
  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> P2GFS;
  P2GFS p2gfs(gv,q22dfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> Q1GFS;
  Q1GFS q1gfs(gv,q12dfem);

  // power grid function space
  typedef Dune::PDELab::PowerGridFunctionSpace<Q1GFS,2,
    Dune::PDELab::GridFunctionSpaceLexicographicMapper> PowerGFS;
  PowerGFS powergfs(q1gfs);

  // make coefficent Vectors
  typedef typename P2GFS::template VectorContainer<double>::Type V;
  V x(p2gfs);
  x = 0.0;
  typedef typename PowerGFS::template VectorContainer<double>::Type VP;
  VP xp(powergfs);
  xp = 0.0;

  // make local function space object
  typename P2GFS::LocalFunctionSpace p2lfs(p2gfs);
  std::vector<double> xl(p2lfs.maxSize());
  typename PowerGFS::LocalFunctionSpace powerlfs(powergfs);
  std::vector<double> xlp(powerlfs.maxSize());

  // loop over elements
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  for (ElementIterator it = gv.template begin<0>();
	   it!=gv.template end<0>(); ++it)
	{
      p2lfs.bind(*it);
      p2lfs.debug();
      p2lfs.vread(x,xl);

      powerlfs.bind(*it);
      powerlfs.debug();
      powerlfs.vread(xp,xlp);
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
