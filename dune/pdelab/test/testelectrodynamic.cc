/** @file
 *  @brief Unit test for `Electrodymic_S` and `Electrodynamic_T` operators
 *
 * This was adapted from Andreas Buhr's circ_in_rect_2d test program.
 *
 * Despite the name, it does not actually do a time-domain simulation; instead
 * it does a frequency-domain simulation.  But is does make use of the
 * operators that were originally intended for time-domain simulations.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/vertexorderfactory.hh>

#include <dune/alugrid/grid.hh>

#include <dune/istl/matrixmarket.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/edges0.5fem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/electrodynamic.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/numericalresidual.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

// store physical entitiy information from GmshReader
struct MeshInfo
{
  // note: we rely on grids where we can use the boundary segment index to
  // subscript the array of physical entities for the boundary intersections.
  std::vector<int> boundary_id_to_physical_entity;
  std::vector<int> element_index_to_physical_entity;
};

// no M_PI in C++ standard, and acos(-1.0) ain't constexpr (may set errno)
constexpr double pi = 3.14159265358979323846264338328;

template <typename Eps, typename Mu>
class Electrodynamic_Full :
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::JacobianBasedAlphaVolume<Electrodynamic_Full<Eps, Mu> >
{
public:
  static constexpr bool doPatternVolume  = true;
  static constexpr bool doAlphaVolume    = true;
  static constexpr bool doLambdaBoundary = true;

  template<class Eps_, class Mu_>
  Electrodynamic_Full(Eps_&& eps, Mu_&& mu, double frequency,
                      const MeshInfo& m, int qorder = 2)
    : edyn_s_lop_(std::forward<Mu_>(mu),   qorder)
    , edyn_t_lop_(std::forward<Eps_>(eps), qorder)
    , meshinfo_(m)
    , frequency_(frequency)
  { /* empty */ }

  template <typename EG, typename LFS, typename X, typename M>
  void jacobian_volume(const EG& eg, const LFS& lfsu, const X& x,
                       const LFS& lfsv, M& mat) const
  {
    // first assemble curl curl part:
    edyn_s_lop_.jacobian_volume(eg, lfsu, x, lfsv, mat);

    // now calculate \omega^2:
    double omega_squared = 2. * pi * frequency_;
    omega_squared *= omega_squared;

    // get a weighted accumulation view:
    auto weightedview = mat.weightedAccumulationView(-omega_squared);

    // and assemble mass part:
    edyn_t_lop_.jacobian_volume(eg, lfsu, x, lfsv, weightedview);
  }

  template <typename IG, typename LFSV, typename R>
  void lambda_boundary(const IG& ig, const LFSV& lfsv_s, R& r_s) const {
    using BasisTraits =
      typename LFSV::Traits::FiniteElementType::Traits::Basis::Traits;
    using Range = typename BasisTraits::Range;

    // see note in MeshInfo
    int boundarySegmentIndex = ig.intersection().boundarySegmentIndex();
    int physical_id =
        meshinfo_.boundary_id_to_physical_entity[boundarySegmentIndex];

    // Non-zero Neumann BC
    if (physical_id == 3) // Port1 (left side) in circ_in_rect.geo
    {
      // tangent in global coordinates
      // assumes 2D and the boundary is aligned with the y-axis
      Range tangent = { 0.0, 1.0 };

      std::vector<Range> phi(lfsv_s.size());
      const auto &geo = ig.geometryInInside();
      for (auto quadraturepoint : quadratureRule(ig.geometry(), 2))
      {
        // map quadraturepoint to tet, then evaluate
        lfsv_s.finiteElement().basis().evaluateFunction(
            geo.global(quadraturepoint.position()), phi);

        auto factor = frequency_ * quadraturepoint.weight();
        for (unsigned i = 0; i < lfsv_s.size(); ++i) {
          r_s.accumulate(lfsv_s, i, factor * (phi[i] * tangent));
        }
      }
    }
  }

private:
  Dune::PDELab::Electrodynamic_S<Mu>  edyn_s_lop_;
  Dune::PDELab::Electrodynamic_T<Eps> edyn_t_lop_;

  const MeshInfo& meshinfo_;
  double frequency_;
};

class BCTypeParam :
  public Dune::PDELab::DirichletConstraintsParameters
{
public:
  BCTypeParam(const MeshInfo& m) : meshinfo_(m) {}

  template <typename IG>
  bool isDirichlet(const IG& intersection,
                   const Dune::FieldVector<typename IG::ctype,
                                           IG::dimension - 1>&) const
  {
    // see note in MeshInfo
    int boundarySegmentIndex =
        intersection.intersection().boundarySegmentIndex();
    int physical_id =
        meshinfo_.boundary_id_to_physical_entity[boundarySegmentIndex];

    // return true;
    // physical ID 2 is dirichlet
    return physical_id == 2;
  }
  const MeshInfo& meshinfo_;
};

int main(int argc, char** argv) {
  // Maybe initialize Mpi
  Dune::MPIHelper::instance(argc, argv);

  // here we start something real:

  std::string nameprefix = "testelectrodynamic";
  std::string meshfilename = "mesh.msh";
  if(argc >= 2)
    meshfilename = argv[1];

  // create the grid
  const int dim = 2;
  using Grid = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;

  MeshInfo meshinfo;
  auto mygrid(Dune::GmshReader<Grid>::read(
      meshfilename, meshinfo.boundary_id_to_physical_entity,
      meshinfo.element_index_to_physical_entity));
  Dune::gridinfo(*mygrid);

  BCTypeParam bctype(meshinfo);

  // now we need a vertex order
  using Geometry = Grid::Codim<0>::Geometry;
  using VertexOrderFactory = Dune::VertexOrderByIdFactory<Grid::GlobalIdSet>;
  VertexOrderFactory vertexOrderFactory(mygrid->globalIdSet());
  using RF = double;
  using FiniteElementMap =
    Dune::PDELab::EdgeS0_5FiniteElementMap<Geometry, VertexOrderFactory, RF>;
  FiniteElementMap finiteElementMap(vertexOrderFactory);

  using GV = Grid::LeafGridView;
  using CON = Dune::PDELab::ConformingDirichletConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using GFS =
    Dune::PDELab::GridFunctionSpace<GV, FiniteElementMap, CON, VBE>;
  GFS gfs(mygrid->leafView(), finiteElementMap);
  gfs.name("solution"); // name for vtkwriter

  using CC = GFS::ConstraintsContainer<double>::Type;
  CC cc;
  Dune::PDELab::constraints(bctype, gfs, cc);
  std::cout << "constrained dofs: " << cc.size() << " of "
            << gfs.globalSize() << std::endl;

  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  MBE mbe(5); // the edge itself + two adjacent elements * two other edges

  // parameter functions
  auto mu = [](const auto &elem, const auto &xl) {
    // this value is the vaccuum permability mu_0 in N/A^2
    return 4e-7 * pi;
  };
  auto eps = [](const auto &elem, const auto &xl) {
    // this value is the vaccuum permittivity eps_0 in F/m
    return 8.854187817e-12;
  };

  // build matrices and write them out for analysis
  auto writeMatrix = [&](auto &&lo, const std::string &fname) {
    using LO = std::decay_t<decltype(lo)>;

    using GO = Dune::PDELab::GridOperator<GFS, GFS, LO, MBE,
                                          double, double, double, CC, CC>;
    GO go(gfs, cc, gfs, cc, lo, mbe);
    using U = typename GO::Traits::Domain;
    U u(gfs, 0);

    typename GO::Jacobian j(go);
    j = 0;
    go.jacobian(u, j);

    // get rhs
    // (strictly speaking, this has no effect besides ensuring that the
    // residual can be computed, we don't check the result in any way
    // currently)
    typename GO::Range rhs(gfs);
    go.residual(u, rhs);

    // TODO: native() isn't found via ADL, despite what CHANGELOG.md says
    Dune::storeMatrixMarket(Dune::PDELab::Backend::native(j), fname);
  };

  using Dune::PDELab::makeLocalOperatorEdynS;
  writeMatrix(makeLocalOperatorEdynS(mu),  nameprefix + "_matrixS.mtx");

  using Dune::PDELab::makeLocalOperatorEdynT;
  writeMatrix(makeLocalOperatorEdynT(eps), nameprefix + "_matrixT.mtx");

  using LocalOperatorFull = Electrodynamic_Full<decltype(eps), decltype(mu)>;
  writeMatrix(LocalOperatorFull(eps, mu, 3e7, meshinfo),
              nameprefix + "_matrix_" + std::to_string((int)3e7) + ".mtx");

  using GOFull = Dune::PDELab::GridOperator<GFS, GFS, LocalOperatorFull, MBE,
                                            double, double, double, CC, CC>;

  // now solve a problem
  using LS = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
  LS ls;
  using U = GOFull::Traits::Domain;
  using SLP = Dune::PDELab::StationaryLinearProblemSolver<GOFull, LS, U>;

  // generate some vtk output
  // TODO: the free leafGridView() function seems broken atm (returning a
  // gridview of a weird type)
  const auto& gv = mygrid->leafView();

  auto dofreq = [&](double frequency) {
    LocalOperatorFull localOperatorFull(eps, mu, frequency, meshinfo);
    GOFull goFull(gfs, cc, gfs, cc, localOperatorFull, mbe);
    U u(gfs, 0);
    SLP slp(goFull, ls, u, 1e-10);
    slp.apply();
    Dune::VTKWriter<GV> myvtkwriter(gv, Dune::VTK::nonconforming);
    Dune::PDELab::addSolutionToVTKWriter(myvtkwriter, gfs, u);
    myvtkwriter.write(nameprefix + "_" + std::to_string((int)frequency),
                      Dune::VTK::appendedraw);
  };

  // scan through frequencies from 1MHz to 300MHz:
  // for (double frequency = 1e6; frequency < 3e8; frequency *= 1.05)
  //   dofreq(frequency);

  // for the unit test:
  // do a single frequency (the one we have a screenshot for)
  dofreq(30426425);

  // Don't handle exceptions, that would prevent typical C++ libraries from
  // printing a meaningful backtrace
}
