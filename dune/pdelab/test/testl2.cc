// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include<iostream>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dune/pdelab.hh>

/* This actually tests a number of things:
 *
 * - The L2 operator and its application to systems
 * - subspaces
 * - Interpolating a vector-valued function into trees with depth > 1
 * - Interpolating into AnalyticGridViewFunctions from dune-functions
 * - DiscreteGridViewFunction (also for subspaces)
 *
 */

// Test that L2 operator runs and that it gives identical results for scalar
// and system spaces
template<typename GV>
void test_l2(const GV& gv)
{
  using RF = typename GV::Grid::ctype;
  constexpr int dim = GV::dimension;

  // instantiate finite element maps
  using Q1FEM = Dune::PDELab::QkLocalFiniteElementMap<GV,RF,RF,1>;
  Q1FEM q1fem(gv);
  using Q2FEM = Dune::PDELab::QkLocalFiniteElementMap<GV,RF,RF,2>;
  Q2FEM q2fem(gv);

  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using Ordering = Dune::PDELab::LexicographicOrderingTag;

  using Q1GFS = Dune::PDELab::GridFunctionSpace<GV,Q1FEM>;
  using Q2GFS = Dune::PDELab::GridFunctionSpace<GV,Q2FEM>;

  // build scalar space
  using ScalarGFS = Q1GFS;
  ScalarGFS scalar_gfs(gv,q1fem);
  scalar_gfs.name("scalar");

  // build system as composite of a power and a scalar space
  Q1GFS q1_gfs(gv,q1fem);
  Q2GFS q2_gfs(gv,q2fem);

  using PowerGFS = Dune::PDELab::PowerGridFunctionSpace<Q1GFS,dim,VBE>;
  PowerGFS power_gfs(q1_gfs);

  using CompositeGFS = Dune::PDELab::CompositeGridFunctionSpace<VBE,Ordering,PowerGFS,Q2GFS>;
  using SystemGFS = CompositeGFS;
  SystemGFS system_gfs(power_gfs,q2_gfs);

  // create scalar vector and interpolate some test data into it
  using ScalarV = Dune::PDELab::Backend::Vector<Q1GFS,RF>;
  ScalarV scalar_v(scalar_gfs,0.0);
  auto f = Dune::Functions::makeAnalyticGridViewFunction(
    [](auto x)
    {
      using std::exp;
      return exp(x.two_norm2());
    },
    gv);
  Dune::PDELab::interpolate(f,scalar_gfs,scalar_v);

  // create system vector and interpolate identical test data into each component
  using SystemV = Dune::PDELab::Backend::Vector<SystemGFS,RF>;
  SystemV system_v(system_gfs,0.0);
  auto g = Dune::Functions::makeAnalyticGridViewFunction(
    [](auto x)
    {
      using std::exp;
      auto r = exp(x.two_norm2());
      return Dune::FieldVector<RF,Dune::TypeTree::TreeInfo<SystemGFS>::leafCount>(r);
    },
    gv);
  Dune::PDELab::interpolate(g,system_gfs,system_v);


  using LOP = Dune::PDELab::L2;
  LOP lop;

  // don't care about correct sizing here...
  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  MBE mbe(10);

  // Calculate L2 residual for scalar space
  using ScalarGO = Dune::PDELab::GridOperator<ScalarGFS,ScalarGFS,LOP,MBE,RF,RF,RF>;
  ScalarGO scalar_go(scalar_gfs,scalar_gfs,lop,mbe);
  ScalarV scalar_r(scalar_gfs,0.0);
  scalar_go.residual(scalar_v,scalar_r);

  // Calculate L2 residual for system space
  using SystemGO = Dune::PDELab::GridOperator<SystemGFS,SystemGFS,LOP,MBE,RF,RF,RF>;
  SystemGO system_go(system_gfs,system_gfs,lop,mbe);
  SystemV system_r(system_gfs,0.0);
  system_go.residual(system_v,system_r);

  // Extract first component from system space
  using SubSpace = Dune::PDELab::GridFunctionSubSpace<SystemGFS,Dune::TypeTree::StaticTreePath<0,0>>;
  SubSpace sub_space(system_gfs);

  // Integrate L2 norm of difference between scalar space and first component of system space
  auto scalar_rf = Dune::PDELab::DiscreteGridViewFunction<ScalarGFS,ScalarV>(scalar_gfs,scalar_v);
  auto sub_rf = Dune::PDELab::DiscreteGridViewFunction<SubSpace,SystemV>(sub_space,system_v);

  RF diff_l2 = 0.0;
  auto scalar_lf = localFunction(scalar_rf);
  auto sub_lf = localFunction(sub_rf);
  for (auto cell : elements(gv))
    {
      scalar_lf.bind(cell);
      sub_lf.bind(cell);
      auto geo = cell.geometry();
      for (auto qp : Dune::PDELab::quadratureRule(geo,4))
        {
          auto diff = scalar_lf(qp.position()) - sub_lf(qp.position());
          diff_l2 += qp.weight() * geo.integrationElement(qp.position()) * diff * diff;
        }
    }

  // should be zero
  assert(diff_l2 < 1e-12);
}


int main(int argc, char** argv)
{
  try {

    Dune::MPIHelper::instance(argc, argv);

    Dune::YaspGrid<2> grid({{1.0,1.0}},{{32,32}});
    test_l2(grid.leafGridView());

    return 0;

  }
  catch (std::exception &e){
    std::cerr << "Dune reported error: " << e.what() << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
