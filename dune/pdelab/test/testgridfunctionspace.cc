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

  // test power
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,2> PGFS2;
  PGFS2 pgfs2(gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,3> PGFS3;
  PGFS3 pgfs3(gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,4> PGFS4;
  PGFS4 pgfs4(gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,5> PGFS5;
  PGFS5 pgfs5(gfs2,gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,6> PGFS6;
  PGFS6 pgfs6(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,7> PGFS7;
  PGFS7 pgfs7(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,8> PGFS8;
  PGFS8 pgfs8(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,9> PGFS9;
  PGFS9 pgfs9(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,10> PGFS10;
  PGFS10 pgfs10(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,17> PGFS17;
  PGFS17 pgfs17(gfs2);
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,17,Dune::PDELab::GridFunctionSpaceBlockwiseMapper> PGFS17B;
  PGFS17B pgfs17b(gfs2);

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

  // test composite
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
    GFS1,PGFS2> CGFS2;
  CGFS2 cgfs2(gfs1,pgfs2);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
    GFS1,PGFS2,CGFS2> CGFS3;
  CGFS3 cgfs3(gfs1,pgfs2,cgfs2);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
    GFS1,PGFS2,CGFS2,CGFS3> CGFS4;
  CGFS4 cgfs4(gfs1,pgfs2,cgfs2,cgfs3);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
    GFS1,PGFS2,CGFS2,CGFS3,CGFS4> CGFS5;
  CGFS5 cgfs5(gfs1,pgfs2,cgfs2,cgfs3,cgfs4);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
    GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5> CGFS6;
  CGFS6 cgfs6(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
    GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5,CGFS6> CGFS7;
  CGFS7 cgfs7(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5,cgfs6);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
    GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5,CGFS6,CGFS7> CGFS8;
  CGFS8 cgfs8(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5,cgfs6,cgfs7);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
    GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5,CGFS6,CGFS7,CGFS8> CGFS9;
  CGFS9 cgfs9(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5,cgfs6,cgfs7,cgfs8);
  
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
