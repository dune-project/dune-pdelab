// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<vector>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/yaspgrid.hh>
#include"../finiteelementmap/p0fem.hh"
#include"../finiteelementmap/edges02dfem.hh"
#include"../gridfunctionspace/gridfunctionspace.hh"
#include"../gridfunctionspace/gridfunctionspaceutilities.hh"
#include"../gridfunctionspace/interpolate.hh"
#include"../common/function.hh"
#include"../common/vtkexport.hh"
#include"gridexamples.hh"


// define some grid functions to interpolate from
template<typename GV, typename RF>
class F
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  F<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F<GV,RF> > BaseT;

  F (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
    typename Traits::DomainType center;
    for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
    center -= x;
	y = exp(-center.two_norm2());
  }
};

template<typename GV, typename RF>
class V 
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,2>,
													  V<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,2> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,V<GV,RF> > BaseT;

  V (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {  
    if (x[0]<0.5)
      {
        y[0] = x[0]/sqrt(2.0);
        y[1] = x[1]/sqrt(2.0);
      }
    else
      {
        y[0] = x[0]/sqrt(2.0);
        y[1] = -x[1]/sqrt(2.0);
      }
  }
};


template<class GV> 
void testedges02d (const GV& gv)
{
  typedef typename GV::Grid::ctype DF;
  const int dim = GV::dimension;

  // instantiate finite element maps
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,double,dim> P0FEM;
  P0FEM p0fem(Dune::GeometryType::simplex);
  typedef Dune::PDELab::EdgeS02DLocalFiniteElementMap<GV,double> EdgeS02DFEM;
  EdgeS02DFEM edges02dfem(gv);
  
  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM> P0GFS; 
  P0GFS p0gfs(gv,p0fem);
  typedef Dune::PDELab::GridFunctionSpace<GV,EdgeS02DFEM> EdgeS02DGFS; 
  EdgeS02DGFS edges02dgfs(gv,edges02dfem);

  // make coefficent Vectors
  typedef typename P0GFS::template VectorContainer<double>::Type P0V;
  P0V p0xg(p0gfs);
  p0xg = 0.0;
  typedef typename EdgeS02DGFS::template VectorContainer<double>::Type EdgeS02DV;
  EdgeS02DV edges02dxg(edges02dgfs);
  edges02dxg = 0.0;

  // construct a grid function
  typedef F<GV,double> FType;
  FType f(gv);
  typedef V<GV,double> VType;
  VType v(gv);

  // do interpolation
  Dune::PDELab::interpolate(f,p0gfs,p0xg);
  Dune::PDELab::interpolateEdge(v,edges02dgfs,edges02dxg);

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<P0GFS,P0V> P0DGF;
  P0DGF p0dgf(p0gfs,p0xg);
  typedef Dune::PDELab::DiscreteGridFunctionEdge<EdgeS02DGFS,EdgeS02DV> EdgeS02DDGF;
  EdgeS02DDGF edges02ddgf(edges02dgfs,edges02dxg);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P0DGF>(p0dgf,"p0"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<EdgeS02DDGF>(edges02ddgf,"edges02d"));
  vtkwriter.write("testedges02d",Dune::VTKOptions::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

#if HAVE_UG
    Dune::SmartPointer<Dune::UGGrid<2> > uggrid(TriangulatedLDomainMaker<Dune::UGGrid<2> >::create());
  	uggrid->globalRefine(4);
    testedges02d(uggrid->leafView());
#endif

#if HAVE_ALBERTA
 	AlbertaLDomain albertagrid;
  	albertagrid.globalRefine(4);
    testedges02d(albertagrid.leafView());
#endif

#if HAVE_ALUGRID
 	ALUUnitSquare alugrid;
  	alugrid.globalRefine(5);
    testedges02d(alugrid.leafView());
#endif

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
