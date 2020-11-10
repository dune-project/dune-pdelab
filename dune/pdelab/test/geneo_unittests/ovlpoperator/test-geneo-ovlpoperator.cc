#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// #include <iostream>
#include <memory>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/power.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>

// includes specific to the GenEO unit tests
#include <dune/pdelab/test/geneo_unittests/utility.hh>


constexpr int DIM{ 2 }; // the dimension of the domain
constexpr int ORD{ 1 }; // the polynomial order used in the Lagrange basis


int main(int argc, char** argv)
{
  // initialize MPI
  Dune::MPIHelper& mpihelper{ Dune::MPIHelper::instance(argc, argv) };

  if (Dune::MPIHelper::isFake)
    DUNE_THROW(Dune::Exception, "GenEO needs to be run in parallel.");

  Dune::ParameterTree ptree;

  // either a single file is passed, which is then used, or no file is passed,
  // in which case a default ini file is read
  if (argc == 2)
    Dune::ParameterTreeParser::readINITree(argv[1], ptree);
  else if (argc == 1)
    Dune::ParameterTreeParser::readINITree("test-geneo-ovlpoperator.ini",
                                           ptree);
  else
    DUNE_THROW(Dune::InvalidStateException, "Error parsing .ini file. "
               "Either provide a single .ini file or no .ini file (default "
               "is used).");

  using Real = double;
  using DF = Real;  // domain field
  using RF = Real;  // range field

  // read the grid parameters form the ini file
  Utility::GridParameters<DF, DIM> gp(ptree);

  // partitioning of the grid
  using Partitioner = Dune::YaspFixedSizePartitioner<DIM>;
  auto partitioner{ std::make_unique<Partitioner>(gp.subdomLayout) };

  // instanitate grid
  using Grid = Dune::YaspGrid<DIM>;
  auto grid{ std::make_shared<Grid>(
               gp.upperRight, gp.nCells, gp.isPeriodic, gp.ovlp,
               Dune::MPIHelper::getCollectiveCommunication(),
               partitioner.get()) };

  // create GridView to finest level (the only level for now)
  using GV = Grid::LevelGridView;
  auto gv{ grid->levelGridView(grid->maxLevel()) };

  // Use Qk Lagrange Basis
  using FEM =
    Dune::PDELab::QkLocalFiniteElementMap<GV, DF, Real, ORD>;
  FEM fem(gv);

  // define the relevant backends
  using Dune::PDELab::Backend::native;
  using VBE =
    Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed, 1>;
  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;

  // constraints at domain boundary, but none at the subdomain boundaries
  using CON = Dune::PDELab::ConformingDirichletConstraints;

  using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE>;
  GFS gfs(gv, fem);

  // set up the constraints containers
  using CC = GFS::template ConstraintsContainer<RF>::Type;
  auto cc{ CC() };

  // the type of boundary condition for the problems
  Utility::BC bcType{ Utility::BC::DirichletAndNeumann };

  // create Poisson problem with heterogeneous coefficient
  using Problem = Utility::PoissonHeterogeneous<GV, DF, RF>;
  DF layerwidth{ 0.1 };
  RF contrast{ 1e2 };
  Problem problem(gv, bcType, layerwidth, contrast);

  // map problem boundary conditions to corresponding entities on the grid
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem>
  bcadapter(gv, problem);

  Dune::PDELab::constraints(bcadapter, gfs, cc, 0);

  // create a local operator and wrap the overlap operator around it
  using LOP = Dune::PDELab::ConvectionDiffusionFEM<Problem, FEM>;
  LOP lop(problem);
  using OvlpLOP = Dune::PDELab::LocalOperatorOvlpRegion<LOP, GFS>;
  OvlpLOP ovlplop(lop, gfs);

  // define corresponding GridOperator
  constexpr int maxNonzeros{ Dune::StaticPower<2*ORD+1, DIM>::power };
  using GO =
    Dune::PDELab::GridOperator<GFS, GFS, OvlpLOP, MBE, DF, RF, RF, CC, CC>;
  auto ovlpgo{ GO(gfs, cc, gfs, cc, ovlplop, MBE(maxNonzeros)) };

  // interpolate the boundary conditions
  Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>
  g(gv, problem);

  using V = Dune::PDELab::Backend::Vector<GFS, RF>;
  V x(gfs, 0.0);

  Dune::PDELab::interpolate(g, gfs, x);

  // assemble the matrix
  GO::Jacobian ovlpmat(ovlpgo);
  ovlpgo.jacobian(x, ovlpmat);
  const auto& natovlpmat{ native(ovlpmat) };

  // load the reference matrix for this process from file
  int mpirank{ mpihelper.rank() };
  std::string filename{ "referencematrix_" + std::to_string(mpirank) };

  auto refmat{ native(GO::Jacobian(ovlpgo)) };
  Dune::loadMatrixMarket(refmat, filename);

  // compare overlap matrix to reference matrix
  if (!Utility::matrices_equal(natovlpmat, refmat))
    DUNE_THROW(Dune::Exception,
               "Matrices on rank " << mpirank << " are not the same");

  return 0;
}
