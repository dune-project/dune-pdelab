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

using namespace Dune;
using namespace Dune::PDELab;

template<typename GridView, typename Basis, int expectedBlockSize>
void backend_test(GridView gridView, Basis basis,
  std::integral_constant<int,expectedBlockSize>)
{
  using Vec = Backend::Vector<Basis,double>;
  std::cout << "----------------------------\n";

  // create vector
  Vec x(basis);
  std::cout << className(Backend::native(x)) << std::endl;

  // assign value
  x = 1.0;

#if 0
  // create discrete function
  typedef Dune::PDELab::DiscreteGridFunction<Basis,Vec> DGF;
  DGF dgf(basis,x);

  // create discrete gridview-function
  typedef Dune::PDELab::DiscreteGridViewFunction<Basis,Vec> DGVF;
  DGVF dgvf(gfs,x);

  // Write solution to VTK
  Dune::VTKWriter<GV> vtkwriter(gfs.gridView());
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  auto adapt = std::make_shared<ADAPT>(dgf,"solution");
  vtkwriter.addVertexData(adapt);
  vtkwriter.write("testgeneo_basis_" + basis_type + "_part_unity_" + part_unity_type);
#endif
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
void test_functions_lagrange_basis_eigen(GridView gridView)
{
  using namespace Functions::BasisFactory;

  auto basis = makeBasis(gridView, lagrange<3>());

  static_assert(models< Functions::Concept::GlobalBasis<std::decay_t<decltype(gridView)>>, decltype(basis)>(), "not a basis");

  using Backend = PDELab::Eigen::VectorBackend;
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
  // try with other backends
  test_functions_lagrange_basis_eigen(gv);
}

int main(int argc, char** argv)
{
  try{

    Dune::MPIHelper::instance(argc,argv);

    // make grid
    const int dim = 2;
    // typedef Dune::UGGrid<dim> GridType;
    typedef Dune::YaspGrid<dim> GridType;

    // Build grid with uneven row/col number to provoke a reentrant corner in parallel case with UG
    Dune::FieldVector<typename GridType::ctype,dim> lowerLeft(0);
    Dune::FieldVector<typename GridType::ctype,dim> upperRight(1);
    std::array<unsigned int,dim> elements;
    std::fill(elements.begin(), elements.end(), 17);

    auto grid = Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);

    test(grid);

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
