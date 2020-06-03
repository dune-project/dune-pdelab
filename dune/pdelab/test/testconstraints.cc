// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <vector>
#include <map>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab.hh>
#include "gridexamples.hh"


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

// define some boundary grid functions to define boundary conditions
template<typename GV>
class B
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,
                                                                                           Dune::FieldVector<int,1> >,
                                                  B<GV> >
{
  GV gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B<GV> > BaseT;

  B (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 1; // all is Dirichlet boundary
  }

  //! get a reference to the GridView
  inline const GV& getGridView () const
  {
    return gv;
  }
};


// generate a P1 function and output it
template<class GV>
void testp1 (const GV& gv)
{
  typedef typename GV::Grid::ctype DF;

  // instantiate finite element maps
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,1> P1FEM;
  P1FEM p1fem(gv);

  // make constrained space
  typedef Dune::PDELab::GridFunctionSpace<GV,P1FEM,Dune::PDELab::ConformingDirichletConstraints> P1GFS;
  P1GFS p1gfs(gv,p1fem);

  // make coefficent Vectors
  using P1V = Dune::PDELab::Backend::Vector<P1GFS,double>;
  P1V p1xg(p1gfs);
  p1xg = 0.0;

  // make constraints map
  typedef typename P1GFS::template ConstraintsContainer<double>::Type P1C;
  P1C p1cg;
  p1cg.clear();

  // interpolate from grid function
  typedef F<GV,double> FType;
  FType f(gv);
  Dune::PDELab::interpolate(f,p1gfs,p1xg);

  // set up constraints from boundary condition function
  typedef B<GV> BType;
  BType b(gv);
  Dune::PDELab::constraints(b,p1gfs,p1cg);

  // set Dirichlet nodes to zero
  Dune::PDELab::set_nonconstrained_dofs(p1cg,0.0,p1xg);

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<P1GFS,P1V> P1DGF;
  P1DGF p1dgf(p1gfs,p1xg);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<P1DGF> >(p1dgf,"p1"));
  vtkwriter.write("testconstraintsp1",Dune::VTK::ascii);
}


// define m component function
template<typename GV, typename RF, int m>
class Fm
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,m>,
                                                  Fm<GV,RF,m> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,m> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Fm<GV,RF,m> > BaseT;

  Fm (const GV& gv) : BaseT(gv) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
							  typename Traits::RangeType& y) const
  {
    for (int i=0; i<m; i++)
      y[i] = i*x.two_norm();
  }
};

// define m component boundary condition function
template<typename GV, int m>
class Bm
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::BoundaryGridFunctionTraits<GV,int,m,
                                                                                           Dune::FieldVector<int,m> >,
                                                  Bm<GV,m> >
{
  GV gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,m,Dune::FieldVector<int,m> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,Bm<GV,m> > BaseT;

  Bm (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typedef Dune::PDELab::IntersectionGeometry<I> IG;

    // map from local coordinates in intersection to global coordinates
    Dune::FieldVector<typename IG::ctype,IG::Geometry::coorddimension>
      xg = ig.geometry().global(x);

    // set boundary condition according to coordinates
    // here we could also use boundaryid etc.
    for (int i=0; i<m; i++)
      if (xg[i%2]>1.0-1E-6)
        y[i] = 1; // Dirichlet boundary
      else
        y[i] = 0;
  }

  //! get a reference to the GridView
  inline const GV& getGridView () const
  {
    return gv;
  }
};


// generate a P1 function and output it
template<class GV>
void testpowerp1 (const GV& gv)
{
  typedef typename GV::Grid::ctype DF;
  const int m=5;

  // instantiate finite element map
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,1> P1FEM;
  P1FEM p1fem(gv);

  // make constrained space
  typedef Dune::PDELab::GridFunctionSpace<GV,P1FEM,Dune::PDELab::ConformingDirichletConstraints> P1GFS;
  P1GFS p1gfs(gv,p1fem);

  // make m components of type P1
  typedef Dune::PDELab::PowerGridFunctionSpace<P1GFS,m,Dune::PDELab::ISTL::VectorBackend<> > P1mGFS;
  P1mGFS p1mgfs(p1gfs);

  // make coefficent Vector
  using P1mV = Dune::PDELab::Backend::Vector<P1mGFS,double>;
  P1mV p1mxg(p1mgfs);
  p1mxg = 0.0;

  // make constraints map
  typedef typename P1mGFS::template ConstraintsContainer<double>::Type P1mC;
  P1mC p1mcg;

  // interpolate from grid function
  typedef Fm<GV,double,m> FmType;
  FmType fm(gv);
  Dune::PDELab::interpolate(fm,p1mgfs,p1mxg);

  // set up constraints from boundary condition function
  typedef Bm<GV,m> BmType;
  BmType bm(gv);
  Dune::PDELab::constraints(bm,p1mgfs,p1mcg);

  // set Dirichlet nodes to zero
  Dune::PDELab::set_constrained_dofs(p1mcg,0.0,p1mxg);

  // subspaces
  typedef Dune::PDELab::GridFunctionSubSpace<P1mGFS,Dune::TypeTree::StaticTreePath<0> > SUB0GFS;
  SUB0GFS sub0gfs(p1mgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<P1mGFS,Dune::TypeTree::StaticTreePath<1> > SUB1GFS;
  SUB1GFS sub1gfs(p1mgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<P1mGFS,Dune::TypeTree::StaticTreePath<2> > SUB2GFS;
  SUB2GFS sub2gfs(p1mgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<P1mGFS,Dune::TypeTree::StaticTreePath<3> > SUB3GFS;
  SUB3GFS sub3gfs(p1mgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<P1mGFS,Dune::TypeTree::StaticTreePath<4> > SUB4GFS;
  SUB4GFS sub4gfs(p1mgfs);

  // make discrete function objects (this is not yet generic enough
  typedef Dune::PDELab::DiscreteGridFunction<SUB0GFS,P1mV> SUB0DGF;
  SUB0DGF sub0dgf(sub0gfs,p1mxg);
  typedef Dune::PDELab::DiscreteGridFunction<SUB1GFS,P1mV> SUB1DGF;
  SUB1DGF sub1dgf(sub1gfs,p1mxg);
  typedef Dune::PDELab::DiscreteGridFunction<SUB2GFS,P1mV> SUB2DGF;
  SUB2DGF sub2dgf(sub2gfs,p1mxg);
  typedef Dune::PDELab::DiscreteGridFunction<SUB3GFS,P1mV> SUB3DGF;
  SUB3DGF sub3dgf(sub3gfs,p1mxg);
  typedef Dune::PDELab::DiscreteGridFunction<SUB4GFS,P1mV> SUB4DGF;
  SUB4DGF sub4dgf(sub4gfs,p1mxg);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addVertexData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<SUB0DGF> >(sub0dgf,"comp 0"));
  vtkwriter.addVertexData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<SUB1DGF> >(sub1dgf,"comp 1"));
  vtkwriter.addVertexData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<SUB2DGF> >(sub2dgf,"comp 2"));
  vtkwriter.addVertexData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<SUB3DGF> >(sub3dgf,"comp 3"));
  vtkwriter.addVertexData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<SUB4DGF> >(sub4dgf,"comp 4"));
  vtkwriter.write("testconstraintspowerp1",Dune::VTK::ascii);
}



int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

#if HAVE_UG
    std::shared_ptr<Dune::UGGrid<2> > uggrid(TriangulatedUnitSquareMaker<Dune::UGGrid<2> >::create());
  	uggrid->globalRefine(4);
    testp1(uggrid->leafGridView());
    testpowerp1(uggrid->leafGridView());
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
