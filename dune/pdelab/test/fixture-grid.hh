#ifndef DUNE_PDELAB_TEST_FIXTURE_GRID_HH
#define DUNE_PDELAB_TEST_FIXTURE_GRID_HH

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/grid/concepts/gridview.hh>

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dune/common/fvector.hh>

#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <tuple>


template<std::size_t dim>
class StructuredGridFixture : public ::testing::Test {
 protected:
  using Grid = Dune::YaspGrid<dim>;

  void SetUp() override {
    Dune::FieldVector<double,dim> L(1.0);
    std::array<int,dim> N;
    std::fill(begin(N), end(N), 16);
    _grid = std::make_unique<Grid>(L,N);
  }

  std::shared_ptr<Grid> _grid;
};

#if HAVE_DUNE_UGGRID
class UnstructuredGridFixture2D : public ::testing::Test {

 protected:
  using Grid = Dune::UGGrid<2>;

  void SetUp() override {
    Dune::GridFactory<Grid> gf;
    Dune::FieldVector<typename Grid::ctype, 2> pos;
    pos[0] =-1.0;  pos[1] =-1.0; gf.insertVertex(pos);
    pos[0] = 0.0;  pos[1] =-1.0; gf.insertVertex(pos);
    pos[0] =-1.0;  pos[1] = 0.0; gf.insertVertex(pos);
    pos[0] = 0.0;  pos[1] = 0.0; gf.insertVertex(pos);
    pos[0] = 1.0;  pos[1] = 0.0; gf.insertVertex(pos);
    pos[0] =-1.0;  pos[1] = 1.0; gf.insertVertex(pos);
    pos[0] = 0.0;  pos[1] = 1.0; gf.insertVertex(pos);
    pos[0] = 1.0;  pos[1] = 1.0; gf.insertVertex(pos);

    auto type = Dune::GeometryTypes::triangle;
    std::vector<unsigned int> vid(3);
    vid[0] = 0;  vid[1] = 1;  vid[2] = 2; gf.insertElement(type, vid);
    vid[0] = 2;  vid[1] = 1;  vid[2] = 3; gf.insertElement(type, vid);
    vid[0] = 2;  vid[1] = 3;  vid[2] = 5; gf.insertElement(type, vid);
    vid[0] = 5;  vid[1] = 3;  vid[2] = 6; gf.insertElement(type, vid);
    vid[0] = 3;  vid[1] = 4;  vid[2] = 6; gf.insertElement(type, vid);
    vid[0] = 6;  vid[1] = 4;  vid[2] = 7; gf.insertElement(type, vid);

    _grid = std::unique_ptr<Grid>(gf.createGrid());
  }

  std::shared_ptr<Grid> _grid;
};
#endif

using StructuredGridFixture1D = StructuredGridFixture<1>;
using StructuredGridFixture2D = StructuredGridFixture<2>;
using StructuredGridFixture3D = StructuredGridFixture<3>;
using StructuredGridFixture4D = StructuredGridFixture<4>;

// list of grids to test
using GridFixtures = ::testing::Types<
    StructuredGridFixture1D
  , StructuredGridFixture2D
  , StructuredGridFixture3D
  , StructuredGridFixture4D
#if HAVE_DUNE_UGGRID
  , UnstructuredGridFixture2D
#endif
>;

#endif // DUNE_PDELAB_TEST_FIXTURE_GRID_HH
