// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<memory>
#include<vector>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include <dune/common/std/make_array.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/yaspgrid.hh>
#include"../finiteelementmap/p0fem.hh"
#include"../finiteelementmap/pkfem.hh"
#include"../gridfunctionspace/gridfunctionspace.hh"
#include"../backend/istl.hh"
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

// generate a Q1 function and output it
template<class GV>
void testpk (const GV& gv)
{
  typedef typename GV::Grid::ctype DF;
  const int dim = GV::dimension;
  const int k = 5;

  // instantiate finite element maps
  auto gt = Dune::GeometryTypes::simplex(dim);
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,double,dim> P0FEM;
  P0FEM p0fem(gt);
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,1> P1FEM;
  P1FEM p1fem(gv);
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,k> PkFEM;
  PkFEM pkfem(gv);

  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
  typedef Dune::PDELab::NoConstraints CON;

  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,CON,VBE> P0GFS;
  P0GFS p0gfs(gv,p0fem);
  typedef Dune::PDELab::GridFunctionSpace<GV,P1FEM,CON,VBE> P1GFS;
  P1GFS p1gfs(gv,p1fem);
  typedef Dune::PDELab::GridFunctionSpace<GV,PkFEM,CON,VBE> PkGFS;
  PkGFS pkgfs(gv,pkfem);

  // make coefficent Vectors
  using P0V = Dune::PDELab::Backend::Vector<P0GFS, double>;
  P0V p0xg(p0gfs);
  p0xg = 0.0;
  using P1V = Dune::PDELab::Backend::Vector<P1GFS, double>;
  P1V p1xg(p1gfs);
  p1xg = 0.0;
  using PkV = Dune::PDELab::Backend::Vector<PkGFS, double>;
  PkV pkxg(pkgfs);
  pkxg = 0.0;

  // construct a grid function
  typedef F<GV,double> FType;
  FType f(gv);

  // do interpolation
  Dune::PDELab::interpolate(f,p0gfs,p0xg);
  Dune::PDELab::interpolate(f,p1gfs,p1xg);
  Dune::PDELab::interpolate(f,pkgfs,pkxg);

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<P0GFS,P0V> P0DGF;
  P0DGF p0dgf(p0gfs,p0xg);
  typedef Dune::PDELab::DiscreteGridFunction<P1GFS,P1V> P1DGF;
  P1DGF p1dgf(p1gfs,p1xg);
  typedef Dune::PDELab::DiscreteGridFunction<PkGFS,PkV> PkDGF;
  PkDGF pkdgf(pkgfs,pkxg);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<P0DGF> >(p0dgf,"p0"));
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<P1DGF> >(p1dgf,"p1"));
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PkDGF> >(pkdgf,"pk"));
  vtkwriter.write("testpk",Dune::VTK::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

#if HAVE_UG
    std::shared_ptr<Dune::UGGrid<2> > uggrid(TriangulatedLDomainMaker<Dune::UGGrid<2> >::create());
  	uggrid->globalRefine(4);
    testpk(uggrid->leafGridView());
#endif

#if HAVE_ALBERTA
 	AlbertaLDomain albertagrid;
  	albertagrid.globalRefine(4);
    testpk(albertagrid.leafGridView());

 	AlbertaReentrantCorner albertagridr;
  	albertagridr->globalRefine(4);
    testpk(albertagridr->leafGridView());
#endif

#if HAVE_DUNE_ALUGRID
    using ALUType = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>;
    auto alugrid = Dune::StructuredGridFactory<ALUType>::createSimplexGrid(Dune::FieldVector<ALUType::ctype, 2>(0.0), Dune::FieldVector<ALUType::ctype, 2>(1.0), Dune::Std::make_array(1u, 1u));
    alugrid->globalRefine(4);
    testpk(alugrid->leafGridView());
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
