// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/pdelab.hh>

double x_component_A(const Dune::FieldVector<double,2> & x)
{
    return x[0];
}

template<typename Coord>
double x_component_B(const Coord & x)
{
    return x[0];
}

template<typename GV, std::size_t range_dim>
struct interpolation_function
  : public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,double,range_dim>,
  interpolation_function<GV,range_dim> >
{

public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,double,range_dim> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,interpolation_function<GV,range_dim> > BaseT;

  interpolation_function (const GV& gv) : BaseT(gv) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    y = x[0];
  }
};


template<class GV>
static void test_interpolate_old_interface(const GV& gv)
{
  // instantiate finite element maps
  auto gt = Dune::GeometryTypes::quadrilateral;
  typedef Dune::PDELab::P0LocalFiniteElementMap<double,double,GV::dimension> P0FEM;
  P0FEM p0fem(gt);
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,double,double,1> Q12DFEM;
  Q12DFEM q12dfem(gv);
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,double,double,2> Q22DFEM;
  Q22DFEM q22dfem(gv);

  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM> P0GFS;
  P0GFS p0gfs(gv,p0fem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> GFS1;
  GFS1 gfs1(gv,q12dfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> GFS2;
  GFS2 gfs2(gv,q22dfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM,Dune::PDELab::NoConstraints,
                                          Dune::PDELab::ISTL::VectorBackend<> > GFS3;
  GFS3 gfs3(gv,q22dfem);

  // test power
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,3,Dune::PDELab::ISTL::VectorBackend<> > PGFS;
  PGFS pgfs(gfs2,gfs2,gfs2);

  // test composite
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    Dune::PDELab::ISTL::VectorBackend<>,
    Dune::PDELab::LexicographicOrderingTag,
    P0GFS,GFS1,GFS2,GFS3
    > CGFS;

  CGFS cgfs(p0gfs,gfs1,gfs2,gfs3);

  // master space - contains power and composite twice, once for interpolation from
  // scalar function, once for interpolation from vector function
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    Dune::PDELab::ISTL::VectorBackend<>,
    Dune::PDELab::LexicographicOrderingTag,
    GFS1,PGFS,CGFS,PGFS,CGFS
    > GFS;

  GFS gfs(gfs1,pgfs,cgfs,pgfs,cgfs);

  // make coefficent Vector
  using V = Dune::PDELab::Backend::Vector<GFS, double>;
  V x(gfs,0.0);

  // build interpolation function tree
  typedef interpolation_function<GV,1> scalar_interpolation_function;
  typedef interpolation_function<GV,3> interpolation_function_3d;
  typedef interpolation_function<GV,4> interpolation_function_4d;

  scalar_interpolation_function sif(gv);
  interpolation_function_3d if3d(gv);
  interpolation_function_4d if4d(gv);

  typedef Dune::PDELab::CompositeGridFunction<
    scalar_interpolation_function,
    scalar_interpolation_function,
    scalar_interpolation_function,
    interpolation_function_3d,
    interpolation_function_4d
    > F;

  F f(sif,sif,sif,if3d,if4d);

  // interpolate
  Dune::PDELab::interpolate(f,gfs,x);
}

template<class GV>
static void test_interpolate(const GV& gv)
{
  // instantiate finite element maps
  auto gt = Dune::GeometryTypes::quadrilateral;
  using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<double,double,GV::dimension>;
  P0FEM p0fem(gt);

  // make a grid function space
  using P0GFS = Dune::PDELab::GridFunctionSpace<GV,P0FEM>;
  P0GFS p0gfs(gv,p0fem);

  // make coefficent Vector
  using V = Dune::PDELab::Backend::Vector<P0GFS, double>;
  V x(p0gfs,0.0);

  // interpolate from old-style scalar function
  {
      using scalar_interpolation_function = interpolation_function<GV,1>;
      scalar_interpolation_function f(gv);
      auto f2 = f;
      Dune::PDELab::interpolate(f,p0gfs,x);
      Dune::PDELab::interpolate(f2,p0gfs,x);
  }
  // interpolate from global function
  {
      auto f = x_component_A;
      auto lf = Dune::Functions::makeAnalyticGridViewFunction(f, gv);
      Dune::PDELab::interpolate(lf,p0gfs,x);
      Dune::PDELab::interpolate(x_component_A,p0gfs,x);
      Dune::PDELab::interpolate(f,p0gfs,x);
  }
  // interpolate from global template function
  {
      using Domain = Dune::FieldVector<typename GV::ctype, GV::dimension>;
      auto f = x_component_B<Domain>;
      auto lf = Dune::Functions::makeAnalyticGridViewFunction(f, gv);
      Dune::PDELab::interpolate(lf,p0gfs,x);
      Dune::PDELab::interpolate(x_component_B<Domain>,p0gfs,x);
      Dune::PDELab::interpolate(f,p0gfs,x);
  }
  // interpolate from lambda
  {
      using Domain = Dune::FieldVector<typename GV::ctype, GV::dimension>;
      auto f = [](const Domain& x) {return x[0];};
      auto lf = Dune::Functions::makeAnalyticGridViewFunction(f, gv);
      Dune::PDELab::interpolate(lf,p0gfs,x);
      Dune::PDELab::interpolate(f,p0gfs,x);
  }
  // interpolate from generic lambda
  {
      using Domain = Dune::FieldVector<typename GV::ctype, GV::dimension>;
      auto f = [](const Domain& x) {return x[0];};
      auto lf = Dune::Functions::makeAnalyticGridViewFunction(f, gv);
      Dune::PDELab::interpolate(lf,p0gfs,x);
      Dune::PDELab::interpolate(f,p0gfs,x);
  }
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    std::cout << "interpolation tests (2D)" << std::endl;
    // need a grid in order to test grid functions
    Dune::FieldVector<double,2> L(1.0);
    std::array<int,2> N(Dune::filledArray<2,int>(1));
    Dune::YaspGrid<2> grid(L,N);
    grid.globalRefine(1);

    test_interpolate_old_interface(grid.leafGridView());
    test_interpolate(grid.leafGridView());

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
