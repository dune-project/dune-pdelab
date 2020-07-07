// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <vector>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>

// generate a Q1 function and output it
template<class GV>
void testq1 (const GV& gv)
{
  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,double,double,1> Q12DFEM;
  Q12DFEM q12dfem(gv);

  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> Q1GFS;
  Q1GFS q1gfs(gv,q12dfem);

  // make coefficent Vectors
  using V = Dune::PDELab::Backend::Vector<Q1GFS, double>;
  V x(q1gfs);
  x = 0.0;
  // Don't do this at home: access raw vector
  Dune::PDELab::Backend::native(x)[3] = 1.0;

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<Q1GFS,V> DGF;
  DGF dgf(q1gfs,x);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgf,"test"));
  vtkwriter.write("q1",Dune::VTK::ascii);
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
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,double,double,1> Q12DFEM;
  Q12DFEM q12dfem(gv);
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,double,double,2> Q22DFEM;
  Q22DFEM q22dfem(gv);

  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> Q1GFS;
  Q1GFS q1gfs(gv,q12dfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> Q2GFS;
  Q2GFS q2gfs(gv,q22dfem);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTL::VectorBackend<>,
    Dune::PDELab::LexicographicOrderingTag,Q1GFS,Q2GFS> CGFS;
  CGFS cgfs(q1gfs,q2gfs);
  typedef Dune::PDELab::PowerGridFunctionSpace<Q2GFS,2,Dune::PDELab::ISTL::VectorBackend<> > PGFS;
  PGFS pgfs(q2gfs,q2gfs);

  // make coefficent Vectors
  using V = Dune::PDELab::Backend::Vector<Q1GFS, double>;
  Q1GFS q1gfs2(gv,q12dfem);
  V xg(q1gfs2);
  xg = 0.0;
  using CV = Dune::PDELab::Backend::Vector<CGFS, double>;
  CV cxg(cgfs);
  cxg = 0.0;
  using PV = Dune::PDELab::Backend::Vector<PGFS, double>;
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
  Dune::PDELab::interpolate(f,q1gfs2,xg);
  Dune::PDELab::interpolate(h,cgfs,cxg); // krass !
  Dune::PDELab::interpolate(h,pgfs,pxg); // krass !

  // subspaces
  typedef Dune::PDELab::GridFunctionSubSpace<CGFS,Dune::TypeTree::StaticTreePath<0> > SUBGFS0;
  SUBGFS0 subgfs0(cgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<CGFS,Dune::TypeTree::StaticTreePath<1> > SUBGFS1;
  SUBGFS1 subgfs1(cgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<PGFS,Dune::TypeTree::StaticTreePath<0> > PSUBGFS0;
  PSUBGFS0 psubgfs0(pgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<PGFS,Dune::TypeTree::StaticTreePath<1> > PSUBGFS1;
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
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF0> >(dgf0,"comp 0"));
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF1> >(dgf1,"comp 1"));
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PDGF0> >(pdgf0,"comp 3"));
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PDGF1> >(pdgf1,"comp 4"));
  vtkwriter.write("interpolated",Dune::VTK::ascii);
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

template<typename GV, typename RF>
class VelocityLinear
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,2>,
    VelocityLinear<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,2> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,VelocityLinear<GV,RF> > BaseT;

  VelocityLinear (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                typename Traits::RangeType& y) const
  {
  y[0] = 1.0 * x[0];
  y[1] = 2.0 * x[1];
  }
};

// generate a Q1 function and output it
template<class GV>
void testtaylorhood (const GV& gv)
{
  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,double,double,1> Q12DFEM;
  Q12DFEM q12dfem(gv);
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,typename GV::ctype,double,2> Q22DFEM;
  Q22DFEM q22dfem(gv);

  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> Q1GFS;
  Q1GFS q1gfs(gv,q12dfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> Q2GFS;
  Q2GFS q2gfs(gv,q22dfem);
  typedef Dune::PDELab::PowerGridFunctionSpace<Q2GFS,GV::dimension,Dune::PDELab::ISTL::VectorBackend<> > VGFS;
  VGFS vgfs(q2gfs);
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTL::VectorBackend<>,
    Dune::PDELab::LexicographicOrderingTag,VGFS,Q1GFS> THGFS;
  THGFS thgfs(vgfs,q1gfs);

  // make coefficent Vector
  using V = Dune::PDELab::Backend::Vector<THGFS, double>;
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
  for (typename V::size_type i=0; i<xg.flatsize(); i++)
    std::cout << "[" << i << ":" << Dune::PDELab::Backend::native(xg)[i] << "] ";
  std::cout << std::endl;

  // check entries
  for (int i=0; i<25; i++)
    if (Dune::PDELab::Backend::native(xg)[i]!=1.0)
      exit(1);
  for (int i=25; i<50; i++)
    if (Dune::PDELab::Backend::native(xg)[i]!=2.0)
      exit(1);
  for (int i=50; i<59; i++)
    if (Dune::PDELab::Backend::native(xg)[i]!=3.0)
      exit(1);
  std::cout << "all entries correct" << std::endl;

  // subspaces
  typedef Dune::PDELab::GridFunctionSubSpace<THGFS,Dune::TypeTree::StaticTreePath<1> > SUBP;
  SUBP subp(thgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<THGFS,Dune::TypeTree::StaticTreePath<0> > SUBV;
  SUBV subv(thgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<SUBV,Dune::TypeTree::StaticTreePath<0> > SUBV0;
  SUBV0 subv0(subv);
  typedef Dune::PDELab::GridFunctionSubSpace<SUBV,Dune::TypeTree::StaticTreePath<1> > SUBV1;
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
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGFV0> >(dgfv0,"v0"));
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGFV1> >(dgfv1,"v1"));
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGFV> >(dgfv,"v"));
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGFP> >(dgfp,"p"));
  vtkwriter.write("taylorhood",Dune::VTK::ascii);
}

template<class GV>
void testgridfunctions (const GV& gv)
{
  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,typename GV::ctype,double,2> Q22DFEM;
  Q22DFEM q22dfem(gv);

  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> Q2GFS;
  Q2GFS q2gfs(gv,q22dfem);
  typedef Dune::PDELab::PowerGridFunctionSpace<Q2GFS,GV::dimension,Dune::PDELab::ISTL::VectorBackend<> > VGFS;
  VGFS vgfs(q2gfs);

  // make coefficent Vector
  using V = Dune::PDELab::Backend::Vector<VGFS, double>;
  V xv(vgfs);
  xv = 0.0;

  // construct a grid function
  typedef VelocityLinear<GV,double> VelocityLinearType;
  VelocityLinearType velocity_lin(gv);
  Dune::PDELab::interpolate(velocity_lin,vgfs,xv);

  // subspaces
  typedef Dune::PDELab::GridFunctionSubSpace<VGFS,Dune::TypeTree::StaticTreePath<0> > SUBV0;
  SUBV0 subv0(vgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<VGFS,Dune::TypeTree::StaticTreePath<1> > SUBV1;
  SUBV1 subv1(vgfs);

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<SUBV0,V> DGFV0;
  DGFV0 dgfv0(subv0,xv);
  typedef Dune::PDELab::DiscreteGridFunction<SUBV1,V> DGFV1;
  DGFV1 dgfv1(subv1,xv);
  typedef Dune::PDELab::VectorDiscreteGridFunction<VGFS,V> DGFV;
  DGFV dgfv(vgfs,xv);

  // scalar gradient gridfunction
  typedef Dune::PDELab::DiscreteGridFunctionGradient<SUBV0,V> DGFV0G;
  DGFV0G dgfv0g(subv0,xv);
  typedef Dune::PDELab::DiscreteGridFunctionGradient<SUBV1,V> DGFV1G;
  DGFV1G dgfv1g(subv1,xv);

  // vector gradient gridfunction
  typedef Dune::PDELab::VectorDiscreteGridFunctionGradient<VGFS,V> DGFVG;
  DGFVG dgfvg(vgfs,xv);

  // divergence gridfunction of vector
  typedef Dune::PDELab::VectorDiscreteGridFunctionDiv<VGFS,V> DGFVD;
  DGFVD dgfvd(vgfs,xv);

  // curl gridfunction of vector
  typedef Dune::PDELab::VectorDiscreteGridFunctionCurl<VGFS,V> DGFVC;
  DGFVC dgfvc(vgfs,xv);

  // check entries of velocity vector
  for (typename V::size_type i=0; i<xv.flatsize(); i++)
    std::cout << "[" << i << ":" << Dune::PDELab::Backend::native(xv)[i] << "] ";
  std::cout << std::endl;

  // values at element centers
  typename DGFV1G::Traits::DomainType x(0.5);
  typename DGFV0::Traits::RangeType v0;
  typename DGFV1::Traits::RangeType v1;
  typename DGFV0G::Traits::RangeType v0grad;
  typename DGFV1G::Traits::RangeType v1grad;
  typename DGFVG::Traits::RangeType vgrad;
  typename DGFVD::Traits::RangeType vdiv;
  typename DGFVC::Traits::RangeType vcurl;

  // evaluate gridfunctions
  for(typename GV::template Codim<0>::Iterator eit = gv.template begin<0>();
      eit != gv.template end<0>(); ++eit)
  {
    dgfv0.evaluate(*eit, x, v0);
    dgfv1.evaluate(*eit, x, v1);
    dgfv0g.evaluate(*eit, x, v0grad);
    dgfv1g.evaluate(*eit, x, v1grad);
    dgfvg.evaluate(*eit, x, vgrad);
    dgfvd.evaluate(*eit, x, vdiv);
    dgfvc.evaluate(*eit, x, vcurl);

    // check matching components of gradients
    if (v0grad[0]!=1.0)
      exit(1);
    if (v1grad[1]!=2.0)
      exit(1);

    // check gradients on the diagonal of vgrad
    if (vgrad[0][0]!=1.0)
      exit(1);
    if (vgrad[1][1]!=2.0)
      exit(1);

    // check divergence
    if(vdiv[0]!=3.0)
      exit(1);

    // check curl
    if(vcurl[0]!=0.0)
      exit(1);
  }
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // need a grid in order to test grid functions
    Dune::FieldVector<double,2> L(1.0);
    std::array<int,2> N(Dune::filledArray<2,int>(1));
    Dune::YaspGrid<2> grid(L,N);
    grid.globalRefine(2);

    testq1(grid.leafGridView());
    testinterpolate(grid.leafGridView());
    testtaylorhood(grid.levelGridView(1));
    testgridfunctions(grid.levelGridView(1));

    std::cout << "All testutilities tests passed." << std::endl;

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
