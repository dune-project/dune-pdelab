#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include"../finiteelementmap/p12dfem.hh"


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

	// instantiate finite element maps
	Dune::PDELab::P12DLocalFiniteElementMap<float,double> p12dmap;

	// test passed
	return 0;

  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
} 
