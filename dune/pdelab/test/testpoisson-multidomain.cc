// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/**
 * @brief Multiple-Domain Grid Function Space test for Poisson
 *
 * * This test tests a multidomain setup realized by having different entity
 *   sets in different parts of the GFS tree.
 * * This is only a conceptual test and implements no coupling between subdomians.
 *   For this reason it doesn't lead to the correct solution of the poisson
 *   problem on the whole domain. Instead, wach sub domain solves poisson individually.
 * * This test does not check against a known solution, we just exect that it runs
 *   without major failures.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/multidomaingrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>

#include "gridexamples.hh"

/**
 * @brief Simple version for a multidomain of ConvectionDiffusionFEM
 * @details Here, we know that child spaces may live in different
 *    entity sets. Thus, this class checks if the received entity is in one
 *    (out of two) child spaces. If check succeeds, we forward assembly to the
 *    base class with the respective child function space.
 *    Notice that interior intersection of domains is not treated at all (if
 *    treated, it should be done in skeleton method), therefore, those interior
 *    boundaries result in Neumann BC.
 */
template <class Problem, class FEM>
struct MultiDomainConvectionDiffusionFEM
    : public Dune::PDELab::ConvectionDiffusionFEM<Problem, FEM> {

  using Base = Dune::PDELab::ConvectionDiffusionFEM<Problem, FEM>;

  MultiDomainConvectionDiffusionFEM(Problem &problem) : Base{problem} {}

  template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG &eg, const LFSU &lfsu, const X &x,
                    const LFSV &lfsv, R &r) const {
    if (lfsu.child(0).gridFunctionSpace().entitySet().contains(eg.entity()))
      Base::alpha_volume(eg, lfsu.child(0), x, lfsv.child(0), r);
    if (lfsu.child(1).gridFunctionSpace().entitySet().contains(eg.entity()))
      Base::alpha_volume(eg, lfsu.child(1), x, lfsv.child(1), r);
  }

  template <typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume(const EG &eg, const LFSU &lfsu, const X &x,
                       const LFSV &lfsv, M &mat) const {
    if (lfsu.child(0).gridFunctionSpace().entitySet().contains(eg.entity()))
      Base::jacobian_volume(eg, lfsu.child(0), x, lfsv.child(0), mat);
    if (lfsu.child(1).gridFunctionSpace().entitySet().contains(eg.entity()))
      Base::jacobian_volume(eg, lfsu.child(1), x, lfsv.child(1), mat);
  }

  template <typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary(const IG &ig, const LFSU &lfsu_s, const X &x_s,
                      const LFSV &lfsv_s, R &r_s) const {
    if (lfsu_s.child(0).gridFunctionSpace().entitySet().contains(ig.inside()))
      Base::alpha_boundary(ig, lfsu_s.child(0), x_s, lfsv_s.child(0), r_s);
    if (lfsu_s.child(1).gridFunctionSpace().entitySet().contains(ig.inside()))
      Base::alpha_boundary(ig, lfsu_s.child(1), x_s, lfsv_s.child(1), r_s);
  }

  // jacobian contribution from boundary
  template <typename IG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_boundary(const IG &ig, const LFSU &lfsu_s, const X &x_s,
                         const LFSV &lfsv_s, M &mat_s) const {
    if (lfsu_s.child(0).gridFunctionSpace().entitySet().contains(ig.inside()))
      Base::jacobian_boundary(ig, lfsu_s.child(0), x_s, lfsv_s.child(0), mat_s);
    if (lfsu_s.child(1).gridFunctionSpace().entitySet().contains(ig.inside()))
      Base::jacobian_boundary(ig, lfsu_s.child(1), x_s, lfsv_s.child(1), mat_s);
  }
};

// ===============================================================
// ===============================================================
// Solve the Poisson equation
//           - \Delta u = f in \Omega,
//                    u = g on \partial\Omega_D
//  -\nabla u \cdot \nu = j on \partial\Omega_N
// ===============================================================
// ===============================================================

// ===============================================================
// Define parameter functions f,g,j and \partial\Omega_D/N
// ===============================================================

template <typename GV, typename RF> class PoissonModelProblem {
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV, RF> Traits;

  //! tensor diffusion constant per cell? return false if you want more than one
  //! evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell() { return true; }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A(const typename Traits::ElementType &e,
    const typename Traits::DomainType &x) const {
    typename Traits::PermTensorType I;
    for (std::size_t i = 0; i < Traits::dimDomain; i++)
      for (std::size_t j = 0; j < Traits::dimDomain; j++)
        I[i][j] = (i == j) ? 1 : 0;
    return I;
  }

  //! velocity field
  typename Traits::RangeType b(const typename Traits::ElementType &e,
                               const typename Traits::DomainType &x) const {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c(const typename Traits::ElementType &e,
    const typename Traits::DomainType &x) const {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f(const typename Traits::ElementType &e,
    const typename Traits::DomainType &x) const {
    const auto &xglobal = e.geometry().global(x);
    if (xglobal[0] > 0.25 && xglobal[0] < 0.375 && xglobal[1] > 0.25 &&
        xglobal[1] < 0.375)
      return 50.0;
    else
      return 0.0;
  }

  //! boundary condition type function
  BCType bctype(const typename Traits::IntersectionType &is,
                const typename Traits::IntersectionDomainType &x) const {
    typename Traits::DomainType xglobal = is.geometry().global(x);

    if (xglobal[1] < 1E-6 || xglobal[1] > 1.0 - 1E-6) {
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    }
    if (xglobal[0] > 1.0 - 1E-6 && xglobal[1] > 0.5 + 1E-6) {
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    }
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g(const typename Traits::ElementType &e,
    const typename Traits::DomainType &x) const {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    xglobal -= 0.5;
    return exp(-xglobal.two_norm2());
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  j(const typename Traits::IntersectionType &is,
    const typename Traits::IntersectionDomainType &x) const {
    typename Traits::DomainType xglobal = is.geometry().global(x);

    if (xglobal[0] > 1.0 - 1E-6 && xglobal[1] > 0.5 + 1E-6) {
      return -5.0;
    } else {
      return 0.0;
    }
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o(const typename Traits::IntersectionType &is,
    const typename Traits::IntersectionDomainType &x) const {
    return 0.0;
  }
};

//===============================================================
// Problem setup and solution
//===============================================================

template <typename Grid, typename FEM, typename CON>
void poisson(const Grid &grid, const FEM &fem_1, const FEM &fem_2,
             std::string filename, int q) {
  using Dune::PDELab::Backend::native;

  // constants and types
  using R = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType::
      Traits::RangeFieldType;

  using GV = typename Grid::SubDomainGrid::LeafGridView;
  using EntitySet = Dune::PDELab::AllEntitySet<GV>;
  // create entity sets for each subdomain
  EntitySet es_1{grid.subDomain(1).leafGridView()};
  EntitySet es_2{grid.subDomain(2).leafGridView()};

  // make function space
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using DomainGFS = Dune::PDELab::GridFunctionSpace<EntitySet, FEM, CON, VBE>;
  // create two grid function spaces with different entity sets!
  DomainGFS gfs_1(es_1, fem_1);
  DomainGFS gfs_2(es_2, fem_2);
  gfs_1.name("domain_1");
  gfs_2.name("domain_2");

  // combine them with a power GFS
  using GFS = Dune::PDELab::PowerGridFunctionSpace<DomainGFS, 2,VBE>;
  GFS gfs{gfs_1, gfs_2};

  // the root node is in charge of the assembly loop, thus we need another
  // entity set that contains entities needed in every child space
  EntitySet assembly_es{grid.subDomain(0).leafGridView()};
  gfs.setEntitySet(assembly_es);

  // make model problem
  typedef PoissonModelProblem<GV, R> Problem;
  Problem problem;

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<R>::Type C;
  C cg;
  cg.clear();
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> bctype(
      problem);
  Dune::PDELab::constraints(bctype, gfs, cg);

  // make local operator
  using LOP = MultiDomainConvectionDiffusionFEM<Problem, FEM>;
  LOP lop{problem};

  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(
      27); // 27 is too large / correct for all test cases, so should work fine

  // make grid operator
  typedef Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, double, double, double,
                                     C, C>
      GridOperator;
  GridOperator gridoperator(gfs, cg, gfs, cg, lop, mbe);

  // make coefficent Vector and initialize it from a function
  // There is some weird shuffling around here - please leave it in,
  // it's there to test the copy constructor and assignment operator of the
  // matrix wrapper
  typedef typename GridOperator::Traits::Domain DV;
  DV x0(gfs, Dune::PDELab::Backend::unattached_container());
  {
    DV x1(gfs);
    DV x2(x1);
    x2 = 0.0;
    x0 = x1;
    x0 = x2;
  }

  // initialize DOFs from Dirichlet extension
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g_1(gfs_1.gridView(), problem);
  G g_2(gfs_2.gridView(), problem);
  auto g = Dune::PDELab::PowerGridFunction<G, 2>{g_1, g_2};
  Dune::PDELab::interpolate(g, gfs, x0);
  Dune::PDELab::set_nonconstrained_dofs(cg, 0.0, x0);

  // represent operator as a matrix
  // There is some weird shuffling around here - please leave it in,
  // it's there to test the copy constructor and assignment operator of the
  // matrix wrapper
  typedef typename GridOperator::Traits::Jacobian M;
  M m;
  {
    Dune::Timer patternTimer;
    M m1(gridoperator);
    std::cout << "pattern creation:" << patternTimer.elapsed() << std::endl;
#ifndef OLD_BACKEND
    std::cout << m1.patternStatistics() << std::endl;
#endif
    M m2(m1);
    m2 = 0.0;
    m = m1;
    m = m2;
  }
  gridoperator.jacobian(x0, m);

  // evaluate residual w.r.t initial guess
  typedef typename GridOperator::Traits::Range RV;
  RV r(gfs);
  r = 0.0;
  gridoperator.residual(x0, r);

  // make ISTL solver
  Dune::MatrixAdapter<typename M::Container, typename DV::Container,
                      typename RV::Container>
      opa(native(m));
  // ISTLOnTheFlyOperator opb(gridoperator);
  Dune::SeqSSOR<typename M::Container, typename DV::Container,
                typename RV::Container>
      ssor(native(m), 1, 1.0);
  Dune::SeqILU<typename M::Container, typename DV::Container,
               typename RV::Container>
      ilu0(native(m), 1.0);
  Dune::Richardson<typename DV::Container, typename RV::Container> richardson(
      1.0);

  Dune::CGSolver<typename DV::Container> solvera(opa, ilu0, 1E-10, 5000, 2);
  // FIXME: Use ISTLOnTheFlyOperator in the second solver again
  Dune::CGSolver<typename DV::Container> solverb(opa, richardson, 1E-10, 5000,
                                                 2);
  Dune::InverseOperatorResult stat;

  // solve the jacobian system
  r *= -1.0; // need -residual
  DV x(gfs, 0.0);
  solvera.apply(native(x), native(r), stat);
  x += x0;

  // output grid function with VTKWriter
  //   using the whole space would be fine, but it would return many zero where
  //   the sub-domains is not defined. So it is more convinient to use
  //   sub-spaces (notice that using gfs_1 and gfs_2 will not work as they do
  //   not contain any ordering)
  {
    using Path = Dune::TypeTree::HybridTreePath<Dune::index_constant<0>>;
    Dune::PDELab::GridFunctionSubSpace<GFS,Path> sub_gfs{gfs};
    Dune::VTKWriter<GV> vtkwriter(grid.subDomain(1).leafGridView());
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter, sub_gfs, x);
    vtkwriter.write(filename + "-1", Dune::VTK::ascii);
  }

  {
    using Path = Dune::TypeTree::HybridTreePath<Dune::index_constant<1>>;
    Dune::PDELab::GridFunctionSubSpace<GFS,Path> sub_gfs{gfs};
    Dune::VTKWriter<GV> vtkwriter(grid.subDomain(2).leafGridView());
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter, sub_gfs, x);
    vtkwriter.write(filename + "-2", Dune::VTK::ascii);
  }
}

// creates a multidomain grid with 3 domain.
//   domain 0: assembly domain, union of all other domains
//   domain 1: all entities where x <  0.5
//   domain 2: all entities where x >= 0.5
template <class HostGrid> auto splitDomain(HostGrid &host_grid) {
  // lets have only 3 sub-domains
  using MDTraits = Dune::mdgrid::FewSubDomainsTraits<HostGrid::dimension, 3>;
  using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDTraits>;
  auto grid = std::make_shared<Grid>(host_grid);

  grid->startSubDomainMarking();
  for (const auto &cell :
       elements(grid->leafGridView(), Dune::Partitions::interior)) {
    grid->addToSubDomain(0, cell);
    if (cell.geometry().center()[0] < 0.5)
      grid->addToSubDomain(1, cell);
    else
      grid->addToSubDomain(2, cell);
  }

  grid->preUpdateSubDomains();
  grid->updateSubDomains();
  grid->postUpdateSubDomains();
  return grid;
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char **argv) {
  try {
    // Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // YaspGrid Q1 2D test
    {
      // make grid
      Dune::FieldVector<double, 2> L(1.0);
      std::array<int, 2> N(Dune::filledArray<2, int>(2));
      Dune::YaspGrid<2> grid(L, N);
      grid.globalRefine(3);

      auto md_grid = splitDomain(grid);
      using Grid = std::decay_t<decltype(*md_grid)>;

      // get view
      using GV = typename Grid::SubDomainGrid::LeafGridView;

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV, DF, double, 1> FEM;
      FEM fem_1(md_grid->subDomain(1).leafGridView());
      FEM fem_2(md_grid->subDomain(2).leafGridView());

      // solve problem
      poisson<Grid, FEM, Dune::PDELab::ConformingDirichletConstraints>(
          *md_grid, fem_1, fem_2, "poisson_yasp_Q1_2d-multidomain", 2);
    }

    // test passed
    return 0;
  } catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
