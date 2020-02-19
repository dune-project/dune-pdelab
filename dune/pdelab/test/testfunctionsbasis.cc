// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/uggrid.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/pdelab.hh>

#include "gridexamples.hh"

namespace Dune::PDELab {
template<typename FSB, typename B>
struct FunctionsBasisInfo
{
  FunctionsBasisInfo(FSB & b) :
    _basis(b) {}

  struct Traits {
    using Basis = FSB;
    using Backend = B;
  };

  FSB & _basis;
};
}


template<typename GridView, typename Basis>
void backend_test(GridView gridView, Basis basis)
{
  // using Vec = Dune::PDELab::Backend::Vector<BasisInfo,double>;
  // Vec x(basis);
}

template<typename GridView>
void test_lagrange_gfs(GridView gridView)
{
  using QkFEM = Dune::PDELab::QkLocalFiniteElementMap<GridView,float,double,1>;
  QkFEM qkfem(gridView);

  using Backend = Dune::PDELab::ISTL::VectorBackend<>;
  using GFS = Dune::PDELab::GridFunctionSpace<GridView,QkFEM,Backend>;
  GFS gfs(gridView,qkfem);

  static_assert(Dune::models< Dune::PDELab::Concept::GridFunctionSpace, GFS>(), "not a grid function space");

  backend_test(gridView, gfs);
}

template<typename GridView>
void test_combined_gfs(GridView gridView)
{
  using QkFEM = Dune::PDELab::QkLocalFiniteElementMap<GridView,float,double,1>;
  QkFEM qkfem(gridView);

  using GFS = Dune::PDELab::GridFunctionSpace<GridView,QkFEM>;
  GFS gfs(gridView,qkfem);

  using Backend = Dune::PDELab::ISTL::VectorBackend<>;
  using PGFS = Dune::PDELab::PowerGridFunctionSpace<GFS,2,Backend>;
  PGFS pgfs(gfs);

  static_assert(Dune::models< Dune::PDELab::Concept::GridFunctionSpace, PGFS>(), "not a grid function space");

  backend_test(gridView, pgfs);
}

template<typename GridView>
void test_functions_lagrange_basis(GridView gridView)
{
  using namespace Dune::Functions::BasisFactory;

  auto basis = makeBasis(gridView, lagrange<3>());

  static_assert(Dune::models< Dune::Functions::Concept::GlobalBasis<std::decay_t<decltype(gridView)>>, decltype(basis)>(), "not a basis");

  using Backend = Dune::PDELab::ISTL::VectorBackend<>;
  using BasisInfo = Dune::PDELab::FunctionsBasisInfo<decltype(basis),Backend>;

  BasisInfo basisInfo(basis);
  backend_test(gridView, basisInfo);
}

template<typename GridView>
void test_functions_combined_basis(GridView gridView)
{
  using namespace Dune::Functions::BasisFactory;

  auto basis = makeBasis(gridView,
    power<2>(
      lagrange<1>(),
      blockedInterleaved())
    );

  // auto basis = makeBasis(gridView,
  //     power<N>(
  //       power<M>(
  //         composite(
  //           lagrange<3>(),
  //           lagrange<1>(),
  //           flatLexicographic()),
  //         flatLexicographic()),
  //       blockedInterleaved())
  //     );

  static_assert(Dune::models< Dune::Functions::Concept::GlobalBasis<GridView>, decltype(basis)>(), "not a basis");

  using Backend = Dune::PDELab::ISTL::VectorBackend<>;
  using BasisInfo = Dune::PDELab::FunctionsBasisInfo<decltype(basis),Backend>;

  BasisInfo basisInfo(basis);
  backend_test(gridView, basisInfo);
}

template<typename Grid>
void test(const std::unique_ptr<Grid>& grid)
{
  auto gv = grid->leafGridView();

  // first make sure we don't break PDELab
  test_lagrange_gfs(gv);
  test_combined_gfs(gv);
  // now try to use dune-functions
  test_functions_lagrange_basis(gv);
  test_functions_combined_basis(gv);
}

int main(int argc, char** argv)
{
  try{

    Dune::MPIHelper::instance(argc,argv);

#if HAVE_UG
    test(UnitTriangleMaker<Dune::UGGrid<2>>::create());
    // test(TriangulatedUnitSquareMaker<Dune::UGGrid<2>>::create());
#endif // HAVE_UG

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
