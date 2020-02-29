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

//! \brief thin wrapper around a dune-functions basis to add information about the backend
template<typename FSB, typename B>
struct FunctionsBasisInfo
{
  FunctionsBasisInfo(FSB & b) :
    _basis(b) {}

  using Basis = FSB;
  using GridView = typename FSB::GridView;

  // try to get around without Traits
  struct Traits {
    using Backend = B;
  };

  FSB & basis() { return _basis; }
  const FSB & basis() const { return _basis; }

  FSB & _basis;
};
}

using namespace Dune;

template<typename GridView, typename Basis, int expectedBlockSize>
void backend_test(GridView gridView, Basis basis,
  std::integral_constant<int,expectedBlockSize>)
{
  using Vec = PDELab::Backend::Vector<Basis,double>;
  // Vec x(basis);
  // std::cout << className<Vec>() << std::endl;
  std::cout << className<typename Vec::Container>() << std::endl;
}

template<typename GridView>
void test_lagrange_gfs(GridView gridView)
{
  using QkFEM = PDELab::QkLocalFiniteElementMap<GridView,float,double,1>;
  QkFEM qkfem(gridView);

  using Backend = PDELab::ISTL::VectorBackend<>;
  using GFS = PDELab::GridFunctionSpace<GridView,QkFEM,Backend>;
  GFS gfs(gridView,qkfem);

  static_assert(models< PDELab::Concept::GridFunctionSpace, GFS>(), "not a grid function space");
  backend_test(gridView, gfs, std::integral_constant<int,1>());
}

template<typename GridView>
void test_combined_gfs(GridView gridView)
{
  using QkFEM = PDELab::QkLocalFiniteElementMap<GridView,float,double,1>;
  QkFEM qkfem(gridView);

  using GFS = PDELab::GridFunctionSpace<GridView,QkFEM>;
  GFS gfs(gridView,qkfem);

  using Backend = PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;
  using PGFS = PDELab::PowerGridFunctionSpace<GFS,4,Backend,
                                              // Dune::PDELab::InterleavedOrderingTag
                                              Dune::PDELab::EntityBlockedOrderingTag
                                              >;
  PGFS pgfs(gfs);

  static_assert(models< PDELab::Concept::GridFunctionSpace, PGFS>(), "not a grid function space");
  backend_test(gridView, pgfs, std::integral_constant<int,2>());
}

template<typename GridView>
void test_functions_lagrange_basis(GridView gridView)
{
  using namespace Functions::BasisFactory;

  auto basis = makeBasis(gridView, lagrange<3>());

  static_assert(models< Functions::Concept::GlobalBasis<std::decay_t<decltype(gridView)>>, decltype(basis)>(), "not a basis");

  using Backend = PDELab::ISTL::VectorBackend<>;
  using BasisInfo = PDELab::FunctionsBasisInfo<decltype(basis),Backend>;

  BasisInfo basisInfo(basis);
  backend_test(gridView, basisInfo, std::integral_constant<int,1>());
}

template<typename GridView>
void test_functions_combined_basis(GridView gridView)
{
  using namespace Functions::BasisFactory;

#if 0
  auto basis = makeBasis(gridView,
    power<4>(
      lagrange<1>(),
      blockedInterleaved())
    );
#else
  auto basis = makeBasis(gridView,
    power<4>(lagrange<1>(), blockedLexicographic()) );
  // -> Dune::Blocked::tag::blocked<Dune::Blocked::tag::flat, Dune::Blocked::tag::flat, Dune::Blocked::tag::flat, Dune::Blocked::tag::flat>
  // -> Dune::BlockVector<Dune::BlockVector<Dune::FieldVector<double, 1>>>
  //    bzw. Dune::BlockVector<Dune::BlockVector<double>>


  // -> Dune::Blocked::tag::leafBlocked<4>
  // -> Dune::BlockVector<Dune::FieldVector<double, 4>>

  // auto basis = makeBasis(gridView,
  //   composite(
  //     power<3>(
  //       lagrange<1>(),
  //       flatInterleaved()),
  //     lagrange<1>(),
  //     blockedInterleaved())
  //   );
#endif

  static_assert(models< Functions::Concept::GlobalBasis<GridView>, decltype(basis)>(), "not a basis");

  using Backend = PDELab::ISTL::VectorBackend<>;
  using BasisInfo = PDELab::FunctionsBasisInfo<decltype(basis),Backend>;

  BasisInfo basisInfo(basis);
  backend_test(gridView, basisInfo, std::integral_constant<int,2>());
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
