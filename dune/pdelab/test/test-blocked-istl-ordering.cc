// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/**
 *  This test captures the bug described in #114
 *
 *  The ordering machinery wrongly decided that a blocked leaf space
 *  with auto-detected block size (descriptor block_size == 0) was not
 *  blocked, causing all manner of problems.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <random>
#include <algorithm>

#include <dune/common/float_cmp.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>

template<typename FlatGFS, typename BlockedGFS>
void check_blocked_backend(const FlatGFS& flat_gfs, const BlockedGFS& blocked_gfs)
{

  // This test assembles the L2 residual for both a blocked and a flat vector
  // with identical, random input data and checks that the entries for both
  // residual vectors are identical

  using Real = typename FlatGFS::Traits::EntitySet::ctype;

  using LOP = Dune::PDELab::L2;
  auto  lop = LOP();

  using FlatVector    = Dune::PDELab::Backend::Vector<FlatGFS,double>;
  using BlockedVector = Dune::PDELab::Backend::Vector<BlockedGFS,double>;

  auto flat_x    = FlatVector(flat_gfs,0.0);
  auto blocked_x = BlockedVector(blocked_gfs,0.0);

  auto rng  = std::mt19937_64();
  auto dist = std::uniform_real_distribution<Real>(0.0,1.0);

  std::generate(flat_x.begin(),flat_x.end(),[&]() { return dist(rng); });
  std::copy(flat_x.begin(),flat_x.end(),blocked_x.begin());

  using MatrixBackend  = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  auto  matrix_backend = MatrixBackend(1);
  using FlatGO         = Dune::PDELab::GridOperator<FlatGFS,FlatGFS,LOP,MatrixBackend,Real,Real,Real>;
  using BlockedGO      = Dune::PDELab::GridOperator<BlockedGFS,BlockedGFS,LOP,MatrixBackend,Real,Real,Real>;

  auto flat_go    = FlatGO(flat_gfs,flat_gfs,lop,matrix_backend);
  auto blocked_go = BlockedGO(blocked_gfs,blocked_gfs,lop,matrix_backend);

  auto flat_r    = FlatVector(flat_gfs,0.0);
  auto blocked_r = BlockedVector(blocked_gfs,0.0);

  flat_go.residual(flat_x,flat_r);
  blocked_go.residual(blocked_x,blocked_r);

  auto r = std::mismatch(
    flat_r.begin(),flat_r.end(),
    blocked_r.begin(),
    [](auto x, auto y) { return Dune::FloatCmp::eq(x,y); }
    );

  if (r.first != flat_r.end())
    DUNE_THROW(Dune::Exception,"Found mismatch between blocked and non-blocked version!");
}


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  using Grid = Dune::YaspGrid<3>;
  using Real = Grid::ctype;
  Grid grid({{1.0,1.0,1.0}},{{4,4,4}});

  constexpr int dim = Grid::dimension;

  using ES = Dune::PDELab::AllEntitySet<Grid::LeafGridView>;
  auto  es = ES(grid.leafGridView());

  using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<Real,Real,2,dim>;
  auto  fem = FEM();

  using Constraints    = Dune::PDELab::NoConstraints;
  using Ordering       = Dune::PDELab::EntityBlockedOrderingTag;
  using FlatBackend    = Dune::PDELab::ISTL::VectorBackend<>;
  using BlockedBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;

  using LeafGFS        = Dune::PDELab::GridFunctionSpace<ES,FEM,Constraints,FlatBackend>;
  using FlatGFS        = LeafGFS;
  using BlockedGFS     = Dune::PDELab::GridFunctionSpace<ES,FEM,Constraints,BlockedBackend>;

  {
    // test with a scalar function space

    auto flat_gfs    = FlatGFS(es,fem);
    auto blocked_gfs = BlockedGFS(es,fem);

    check_blocked_backend(flat_gfs,blocked_gfs);
  }

  {
    // test with a nested function space

    auto leaf_gfs    = FlatGFS(es,fem);

    using FlatGFS    = Dune::PDELab::PowerGridFunctionSpace<LeafGFS,2,FlatBackend,Ordering>;
    using BlockedGFS = Dune::PDELab::PowerGridFunctionSpace<LeafGFS,2,BlockedBackend,Ordering>;

    auto flat_gfs    = FlatGFS(leaf_gfs);
    auto blocked_gfs = BlockedGFS(leaf_gfs);

    check_blocked_backend(flat_gfs,blocked_gfs);
  }

  return 0;
}
