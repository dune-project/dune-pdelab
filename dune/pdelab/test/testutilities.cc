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
#include"../finiteelementmap/edger02dfem.hh"
#include"../gridfunctionspace/gridfunctionspace.hh"
#include"../gridfunctionspace/gridfunctionspaceutilities.hh"
#include"../gridfunctionspace/interpolate.hh"
#include"../common/function.hh"
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
  typedef Dune::PDELab::EdgeR02DLocalFiniteElementMap<typename GV::Grid::ctype,double> EdgeR02DFEM;
  EdgeR02DFEM edger02dfem;
  
  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,EdgeR02DFEM> EGFS;
  EGFS egfs(gv,edger02dfem);

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
class G
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  G<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G<GV,RF> > BaseT;

  G (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
    y = sin(3.1415*x[0])*cos(3*3.1415*x[1]);
  }
};

// generate a Q1 function and output it
template<class GV> 
void testinterpolate (const GV& gv)
{
  // instantiate finite element maps
  typedef Dune::PDELab::Q12DLocalFiniteElementMap<typename GV::Grid::ctype,double> Q12DFEM;
  Q12DFEM q12dfem;
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<typename GV::Grid::ctype,double> Q22DFEM;
  Q22DFEM q22dfem;
  
  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> Q1GFS;
  Q1GFS q1gfs(gv,q12dfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> Q2GFS;
  Q2GFS q2gfs(gv,q22dfem);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
    Q1GFS,Q2GFS> CGFS;
  CGFS cgfs(q1gfs,q2gfs);
  typedef Dune::PDELab::PowerGridFunctionSpace<Q2GFS,2> PGFS;
  PGFS pgfs(q2gfs,q2gfs);

  // make coefficent Vectors
  typedef typename Q1GFS::template VectorContainer<double>::Type V;
  V xg(q1gfs);
  xg = 0.0;
  typedef typename CGFS::template VectorContainer<double>::Type CV;
  CV cxg(cgfs);
  cxg = 0.0;
  typedef typename PGFS::template VectorContainer<double>::Type PV;
  PV pxg(pgfs);
  pxg = 0.0;

  // construct a grid function
  typedef F<GV,double> FType;
  FType f(gv);
  typedef G<GV,double> GType;
  GType g(gv);
  typedef Dune::PDELab::CompositeGridFunction<FType,GType> HType;
  HType h(f,g);

  // do interpolation
  Dune::PDELab::interpolate(f,q1gfs,xg);
  Dune::PDELab::interpolate(h,cgfs,cxg); // krass !
  Dune::PDELab::interpolate(h,pgfs,pxg); // krass !

  // subspaces
  typedef Dune::PDELab::GridFunctionSubSpace<CGFS,0> SUBGFS0;
  SUBGFS0 subgfs0(cgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<CGFS,1> SUBGFS1;
  SUBGFS1 subgfs1(cgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<PGFS,0> PSUBGFS0;
  PSUBGFS0 psubgfs0(pgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<PGFS,1> PSUBGFS1;
  PSUBGFS1 psubgfs1(pgfs);

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<SUBGFS0,CV> DGF0;
  DGF0 dgf0(subgfs0,cxg);
  typedef Dune::PDELab::DiscreteGridFunction<SUBGFS1,CV> DGF1;
  DGF1 dgf1(subgfs1,cxg);
  typedef Dune::PDELab::DiscreteGridFunction<PSUBGFS0,PV> PDGF0;
  PDGF0 pdgf0(psubgfs0,pxg);
  typedef Dune::PDELab::DiscreteGridFunction<PSUBGFS1,PV> PDGF1;
  PDGF1 pdgf1(psubgfs1,pxg);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF0>(dgf0,"comp 0"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF1>(dgf1,"comp 1"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF0>(pdgf0,"comp 3"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF1>(pdgf1,"comp 4"));
  vtkwriter.write("interpolated",Dune::VTKOptions::ascii);
}

template<typename GV, typename RF>
class One
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  One<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,One<GV,RF> > BaseT;

  One (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
    y = 1.0;
  }
};

template<typename GV, typename RF>
class Two
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  Two<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Two<GV,RF> > BaseT;

  Two (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
    y = 2.0;
  }
};

template<typename GV, typename RF>
class Three
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  Three<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Three<GV,RF> > BaseT;

  Three (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
    y = 3.0;
  }
};

template<typename GV, typename RF>
class Velocity 
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,2>,
													  Velocity<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,2> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Velocity<GV,RF> > BaseT;

  Velocity (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {  
	y[0] = 1.0;
	y[1] = 2.0;
  }
};

// generate a Q1 function and output it
template<class GV> 
void testtaylorhood (const GV& gv)
{
  // instantiate finite element maps
  typedef Dune::PDELab::Q12DLocalFiniteElementMap<typename GV::Grid::ctype,double> Q12DFEM;
  Q12DFEM q12dfem;
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<typename GV::Grid::ctype,double> Q22DFEM;
  Q22DFEM q22dfem;
  
  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> Q1GFS;
  Q1GFS q1gfs(gv,q12dfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> Q2GFS;
  Q2GFS q2gfs(gv,q22dfem);
  typedef Dune::PDELab::PowerGridFunctionSpace<Q2GFS,GV::dimension> VGFS;
  VGFS vgfs(q2gfs);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,
    VGFS,Q1GFS> THGFS;
  THGFS thgfs(vgfs,q1gfs);

  // make coefficent Vector
  typedef typename THGFS::template VectorContainer<double>::Type V;
  V xg(thgfs);
  xg = 0.0;

  // construct a grid function
  typedef One<GV,double> OneType;
  OneType v0(gv);
  typedef Two<GV,double> TwoType;
  TwoType v1(gv);
  typedef Three<GV,double> ThreeType;
  ThreeType p(gv);
  typedef Dune::PDELab::CompositeGridFunction<OneType,TwoType> VType;
  VType v(v0,v1);
  typedef Dune::PDELab::CompositeGridFunction<VType,ThreeType> THType;
  THType th(v,p);
  typedef Velocity<GV,double> VelocityType;
  VelocityType velocity(gv);
  typedef Dune::PDELab::CompositeGridFunction<VelocityType,ThreeType> AlternativeTHType;
  AlternativeTHType alternativeth(velocity,p);
  
  // do interpolation
  Dune::PDELab::interpolate(th,thgfs,xg);
  Dune::PDELab::interpolate(alternativeth,thgfs,xg);

  // check entries of global vector
  for (typename V::size_type i=0; i<xg.size(); i++)
    std::cout << "[" << i << ":" << xg[i] << "] ";
  std::cout << std::endl;

  // check entries
  for (int i=0; i<25; i++)
    if (xg[i]!=1.0)
      exit(1);
  for (int i=25; i<50; i++)
    if (xg[i]!=2.0)
      exit(1);
  for (int i=50; i<59; i++)
    if (xg[i]!=3.0)
      exit(1);
  std::cout << "all entries correct" << std::endl;

  // subspaces
  typedef Dune::PDELab::GridFunctionSubSpace<THGFS,1> SUBP;
  SUBP subp(thgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<THGFS,0> SUBV;
  SUBV subv(thgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<SUBV,0> SUBV0;
  SUBV0 subv0(subv);
  typedef Dune::PDELab::GridFunctionSubSpace<SUBV,1> SUBV1;
  SUBV1 subv1(subv);


  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<SUBV0,V> DGFV0;
  DGFV0 dgfv0(subv0,xg);
  typedef Dune::PDELab::DiscreteGridFunction<SUBV1,V> DGFV1;
  DGFV1 dgfv1(subv1,xg);
  typedef Dune::PDELab::VectorDiscreteGridFunction<SUBV,V> DGFV;
  DGFV dgfv(subv,xg);
  typedef Dune::PDELab::DiscreteGridFunction<SUBP,V> DGFP;
  DGFP dgfp(subp,xg);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGFV0>(dgfv0,"v0"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGFV1>(dgfv1,"v1"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGFV>(dgfv,"v"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGFP>(dgfp,"p"));
  vtkwriter.write("taylorhood",Dune::VTKOptions::ascii);
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
	Dune::YaspGrid<2> grid(L,N,B,0);
    grid.globalRefine(5);

	testq1(grid.leafView());
	testedger(grid.leafView());
    testinterpolate(grid.leafView());
    testtaylorhood(grid.levelView(1));

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
