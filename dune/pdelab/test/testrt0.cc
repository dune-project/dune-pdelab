// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<vector>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/yaspgrid.hh>
#include"../finiteelementmap/p0fem.hh"
#include"../finiteelementmap/rt02dfem.hh"
#include"../finiteelementmap/rt0q2dfem.hh"
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
void testrt0 (const GV& gv)
{
  typedef typename GV::Grid::ctype DF;
  const int dim = GV::dimension;

  // instantiate finite element maps
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,double,dim> P0FEM;
  P0FEM p0fem(Dune::GeometryType::simplex);
  //P0FEM p0fem(Dune::GeometryType::cube);
  typedef Dune::PDELab::RT02DLocalFiniteElementMap<GV,DF,double> RT0FEM;
  //typedef Dune::PDELab::RT0Q2DLocalFiniteElementMap<GV,DF,double> RT0FEM;
  RT0FEM rt0fem(gv);
  
  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM> P0GFS; 
  P0GFS p0gfs(gv,p0fem);
  typedef Dune::PDELab::GridFunctionSpace<GV,RT0FEM> RT0GFS; 
  RT0GFS rt0gfs(gv,rt0fem);

  // make coefficent Vectors
  typedef typename P0GFS::template VectorContainer<double>::Type P0V;
  P0V p0xg(p0gfs);
  p0xg = 0.0;
  typedef typename RT0GFS::template VectorContainer<double>::Type RT0V;
  RT0V rt0xg(rt0gfs);
  rt0xg = 0.0;

  // construct a grid function
  typedef F<GV,double> FType;
  FType f(gv);
  typedef V<GV,double> VType;
  VType v(gv);
  typedef Dune::PDELab::PiolaBackwardAdapter<VType> RVType;
  RVType rv(v);

  // do interpolation
  Dune::PDELab::interpolate(f,p0gfs,p0xg);
  Dune::PDELab::interpolate(rv,rt0gfs,rt0xg);

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<P0GFS,P0V> P0DGF;
  P0DGF p0dgf(p0gfs,p0xg);
  typedef Dune::PDELab::DiscreteGridFunctionPiola<RT0GFS,RT0V> RT0DGF;
  RT0DGF rt0dgf(rt0gfs,rt0xg);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P0DGF>(p0dgf,"p0"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<RT0DGF>(rt0dgf,"rt0"));
  vtkwriter.write("testrt0",Dune::VTKOptions::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    if (false)
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(1);
      Dune::FieldVector<bool,2> B(false);
      Dune::YaspGrid<2> grid(L,N,B,0);
      grid.globalRefine(5);
      testrt0(grid.leafView());
      return 0;
    }

#if HAVE_UG
    Dune::shared_ptr<Dune::UGGrid<2> > uggrid(TriangulatedLDomainMaker<Dune::UGGrid<2> >::create());
  	uggrid->globalRefine(4);
    testrt0(uggrid->leafView());
#endif

#if HAVE_ALBERTA
 	AlbertaLDomain albertagrid;
  	albertagrid.globalRefine(4);
    testrt0(albertagrid.leafView());
#endif

#if HAVE_ALUGRID
 	ALUUnitSquare alugrid;
  	alugrid.globalRefine(5);
    testrt0(alugrid.leafView());
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
