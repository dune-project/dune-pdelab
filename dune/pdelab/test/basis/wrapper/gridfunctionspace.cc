#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "../../fixture-grid.hh"
#include "../ordering.hh"

#include <dune/pdelab/basis/wrapper/gridfunctionspace.hh>
#include <dune/pdelab/basis/constraints/unconstrained.hh>
#include <dune/pdelab/basis/constraints/composite.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/backend/istl.hh>

TEST_F(StructuredGridFixture2D, TestPDELabQ1Ordering) {
  using GV =  typename Grid::LeafGridView;
  using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV, double, double, 1>;

  auto fem = std::make_shared<FEM>(_grid->leafGridView());
  using CON = Dune::PDELab::NoConstraints;

  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed> > GFS;
  auto q1_gfs = std::make_shared<GFS>(_grid->leafGridView(),fem);
  q1_gfs->update();

  Dune::PDELab::test_ordering( Dune::PDELab::Legacy::Basis{q1_gfs, Dune::PDELab::Unconstrained{}} );
}

TEST_F(StructuredGridFixture2D, TestPDELabQ2Ordering) {
  using GV =  typename Grid::LeafGridView;
  using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV, double, double, 2>;

  auto fem = std::make_shared<FEM>(_grid->leafGridView());
  using CON = Dune::PDELab::NoConstraints;

  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed> > GFS;
  auto q2_gfs = std::make_shared<GFS>(_grid->leafGridView(),fem);
  q2_gfs->update();

  Dune::PDELab::test_ordering( Dune::PDELab::Legacy::Basis{q2_gfs, Dune::PDELab::Unconstrained{}} );
}


TEST_F(StructuredGridFixture2D, TestPDELabQ1x3Ordering) {
  using GV =  typename Grid::LeafGridView;
  using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV, double, double, 1>;

  auto fem = std::make_shared<FEM>(_grid->leafGridView());
  using CON = Dune::PDELab::NoConstraints;

  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none> > GFS;
  GFS q1_gfs(_grid->leafGridView(),fem);

  using PGFS = Dune::PDELab::PowerGridFunctionSpace<GFS,3,Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>>;
  auto q1x3_gfs = std::make_shared<PGFS>(q1_gfs, q1_gfs, q1_gfs);
  q1x3_gfs->update();

  auto con = Dune::PDELab::makeCompositeConstraints(std::array<Dune::PDELab::Unconstrained,3>());
  Dune::PDELab::test_ordering( Dune::PDELab::Legacy::Basis{q1x3_gfs, con} );
}


TEST_F(StructuredGridFixture2D, TestPDELabQ1x3BlockedOrdering) {
  using GV =  typename Grid::LeafGridView;
  using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV, double, double, 1>;

  auto fem = std::make_shared<FEM>(_grid->leafGridView());
  using CON = Dune::PDELab::NoConstraints;

  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed> > GFS;
  GFS q1_gfs(_grid->leafGridView(),fem);

  using PGFS = Dune::PDELab::PowerGridFunctionSpace<GFS,3,Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>>;
  auto q1x3_gfs = std::make_shared<PGFS>(q1_gfs, q1_gfs, q1_gfs);
  q1x3_gfs->update();

  auto con = Dune::PDELab::makeCompositeConstraints(std::array<Dune::PDELab::Unconstrained,3>());
  Dune::PDELab::test_ordering( Dune::PDELab::Legacy::Basis{q1x3_gfs, con} );
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  Dune::MPIHelper::instance(argc, argv);
  return RUN_ALL_TESTS();
}
