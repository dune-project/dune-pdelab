#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <memory>
#include <algorithm>

#include <dune/pdelab.hh>
#include <dune/pdelab/backend/istl/geneo/geneo.hh>
#include <dune/pdelab/test/geneo_unittests/utility.hh>


constexpr int DIM{ 2 };
constexpr int ORD{ 1 };


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
    Dune::ParameterTreeParser::readINITree("test-geneo-twolevelschwarz.ini",
                                           ptree);
  else
    DUNE_THROW(Dune::InvalidStateException, "Error parsing .ini file. "
               "Either provide a single .ini file or no .ini file (default "
               "is used).");

  using Real = double;
  using DF = Real;  // domain field
  using RF = Real;  // range field

  // read parameter for the problem coefficient
  DF layerwidth{ ptree.get<Real>("coefficient.layerwidth") };
  RF contrast{ ptree.get<Real>("coefficient.contrast") };

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
  using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV, DF, Real, ORD>;
  FEM fem(gv);

  // define the relevant backends
  using VBE =
    Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed, 1>;
  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  using Dune::PDELab::Backend::native;

  // constraints at the domain boundary as well as subdomain boundaries, needed
  // for the partitions of unity and the additive Schwarz method.
  using CON_INT = Dune::PDELab::OverlappingConformingDirichletConstraints;

  // constraints at the domain boundary, but none at the subdomain boundaries.
  // Needed for the setup of the GenEO eigenproblem.
  using CON_EXT = Dune::PDELab::ConformingDirichletConstraints;

  // set up corresponding function spaces
  using GFS_INT = Dune::PDELab::GridFunctionSpace<GV, FEM, CON_INT, VBE>;
  GFS_INT gfsint(gv, fem);
  using GFS_EXT = Dune::PDELab::GridFunctionSpace<GV, FEM, CON_EXT, VBE>;
  GFS_EXT gfsext(gv, fem);

  // create Poisson problem with heterogeneous coefficient
  using Problem = Utility::PoissonHeterogeneous<GV, DF, RF>;
  Utility::BC bctype{ Utility::BC::DirichletAndNeumann };
  Problem problem(gv, bctype, layerwidth, contrast);

  // instantiate the constraints containers
  using CC_INT = typename GFS_INT::template ConstraintsContainer<RF>::Type;
  auto cctls{ CC_INT() };
  auto ccpou{ CC_INT() };
  using CC_EXT = typename GFS_EXT::template ConstraintsContainer<RF>::Type;
  auto ccext{ CC_EXT() };

  // exclude constraints at domain boundary for partition of unity
  Dune::PDELab::NoDirichletConstraintsParameters poubcadapter;

  // map constraints from the problem to grid for GenEO and two level Schwarz
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem>
  bcadapter(gv, problem);

  // set up constraints in the constraints containers
  Dune::PDELab::constraints(bcadapter, gfsint, cctls, 0);
  Dune::PDELab::constraints(bcadapter, gfsext, ccext, 0);
  Dune::PDELab::constraints(poubcadapter, gfsint, ccpou, 0);

  // create the partition of unity
  using V_INT = Dune::PDELab::Backend::Vector<GFS_INT, RF>;
  auto pou{ std::make_shared<V_INT>
              (standardPartitionOfUnity<V_INT>(gfsint, ccpou)) };

  // define the local operator for the problem
  using LOP = Dune::PDELab::ConvectionDiffusionFEM<Problem, FEM>;
  LOP lop(problem);

  // set up the overlap operator
  using LOP_OVLP = Dune::PDELab::LocalOperatorOvlpRegion<LOP, GFS_EXT>;
  LOP_OVLP lopovlp(lop, gfsext);

  // define corresponding GridOperators
  constexpr int maxnzs{ Dune::StaticPower<2*ORD+1, DIM>::power };
  using GO_INT =
    Dune::PDELab::GridOperator<GFS_INT, GFS_INT, LOP, MBE, DF, RF, RF, CC_INT,
        CC_INT>;
  GO_INT goint(gfsint, cctls, gfsint, cctls, lop, MBE(maxnzs));
  using GO_OVLP =
    Dune::PDELab::GridOperator<GFS_EXT, GFS_EXT, LOP_OVLP, MBE, DF, RF, RF,
        CC_EXT, CC_EXT>;
  GO_OVLP goovlp(gfsext, ccext , gfsext, ccext, lopovlp, MBE(maxnzs));
  using GO_EXT =
    Dune::PDELab::GridOperator<GFS_EXT, GFS_EXT, LOP, MBE, DF, RF, RF, CC_EXT,
        CC_EXT>;
  GO_EXT goext(gfsext, ccext, gfsext, ccext, lop, MBE(maxnzs));

  // create degrees of freedom vector on the fine grid
  V_INT x(gfsint, 0.0);

  // use the same storage, but different constraints for the degrees of freedom
  // vector without subdomain boundary constraints
  using V_EXT = Dune::PDELab::Backend::Vector<GFS_EXT, RF>;
  V_EXT xext(gfsext, Dune::PDELab::Backend::unattached_container());
  xext.attach(x.storage());

  // interpolate the boundary conditions
  Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>
  g(gv, problem);
  Dune::PDELab::interpolate(g, gfsint, x);

  // matrix types
  using M_INT = GO_INT::Jacobian;
  using M_EXT = GO_EXT::Jacobian;

  // assemble the fine-grid stiffness matrix with subdmain constraints
  M_INT tlsmat(goint);
  goint.jacobian(x, tlsmat);

  // assemble the fine grid stiffness matrix with homogeneous Neumann boundary
  // conditions for the GenEO eigenproblem
  M_EXT epmat(goext);
  goext.jacobian(xext, epmat);

  // assemble the overlap matrix
  M_EXT ovlpmat(goovlp);
  goovlp.jacobian(xext, ovlpmat);

  // eigenvalue threshold
  int H{ gp.nCells[0] / gp.subdomLayout[0] + gp.ovlp };
  int d{ 2 * gp.ovlp };
  double evt{ double(d) / H };

  // number of eigenvalues computed, will be set in local basis
  int nev{ 20 };

  // compute the local basis
  using LB = Dune::PDELab::GenEOBasis<GFS_INT, M_EXT, V_INT, 1>;
  auto localbasis{ std::make_shared<LB>(gfsint, epmat, ovlpmat, evt, *pou,
                                        nev, 25) };

  // communicate local bases in order to construct the global coarse space
  using PH = Dune::PDELab::ISTL::ParallelHelper<GFS_INT>;
  PH phelper(gfsint);
  using CS =
    Dune::PDELab::SubdomainProjectedCoarseSpace<GFS_INT, M_EXT, V_INT, PH>;
  auto coarsespace{ std::make_shared<CS>(gfsint, epmat, localbasis, phelper) };

  // construct the preconditioner
  using POP = Dune::PDELab::OverlappingOperator<CC_INT, M_INT, V_INT, V_INT>;
  auto pop{ std::make_shared<POP>(cctls, tlsmat) };
  using PSP = Dune::PDELab::OverlappingScalarProduct<GFS_INT, V_INT>;
  PSP psp(gfsint, phelper);
  using TLS =
    Dune::PDELab::ISTL::TwoLevelOverlappingAdditiveSchwarz<GFS_INT, M_INT,
        V_INT, V_INT>;
  auto prec{ std::make_shared<TLS>(gfsint, tlsmat, coarsespace, true, false) };

  // assemble the residual
  V_INT r(gfsint, 0.0);
  goint.residual(x, r);

  // solve the defect equation Av = r, using preconditioned CG
  V_INT v(gfsint, 0.0);
  Dune::InverseOperatorResult solverinfo;
  using Solver = Dune::CGSolver<V_INT>;
  auto solver{ std::make_shared<Solver>(*pop, psp, *prec, 1e-6, 1e3, 0, 1) };
  solver->apply(v, r, solverinfo);
  x -= v;

  int mpirank{ mpihelper.rank() };

  std::cout << "condition estimate from rank " << mpirank << ": "
            << solverinfo.condition_estimate << std::endl;

  // compute theoretical bound
  double k0{ 2.0 };
  double tbound{ (1.0 + k0) * (2.0 + k0 * (2.0 * k0 + 1.0)
                               * (1.0 + 1.0 / evt)) };
  std::cout << "theoretical bound for rank " << mpirank << ": " << tbound
            << std::endl;

  if (solverinfo.condition_estimate > tbound)
    DUNE_THROW(Dune::Exception, "condition estimate larger than theoretical"
               " bound");

  return 0;
}
