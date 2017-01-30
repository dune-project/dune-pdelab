/** @file
 *  @brief Unit test for `Electrodymic_S` and `Electrodynamic_T` operators
 *
 * This was adapted from Andreas Buhr's circ_in_rect_2d test program.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>         // We use exceptions
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/vertexorderfactory.hh>

#include <dune/alugrid/grid.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/io.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/constraints/noconstraints.hh>
#include <dune/pdelab/finiteelementmap/edges0.5fem.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/electrodynamic.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

// store physical entitiy information from GmshReader
class MeshInfo {
public:
  std::vector<int> boundary_id_to_physical_entity;
  std::vector<int> element_index_to_physical_entity;
};

const double pi = 3.14159265358979323846264338328;

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

  Electrodynamic_Full(const Eps& eps, const Mu& mu, double frequency,
                      const MeshInfo& m, int qorder = 2)
    : edyn_s_lop_(mu, qorder)
    , edyn_t_lop_(eps, qorder)
    , meshinfo_(m)
    , frequency_(frequency) {
    // empty
  }

  template <typename EG, typename LFS, typename X, typename M>
  void jacobian_volume(const EG& eg, const LFS& lfsu, const X& x,
                       const LFS& lfsv, M& mat) const {
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

        for (unsigned i = 0; i < lfsv_s.size(); ++i) {
          r_s.accumulate(lfsv_s, i, frequency_ *
                                        (phi[i] * tangent) *
                                        quadraturepoint.weight());
        }
      }
    }
  }

private:
  Dune::PDELab::Electrodynamic_S<Mu> edyn_s_lop_;
  Dune::PDELab::Electrodynamic_T<Eps> edyn_t_lop_;

  const MeshInfo& meshinfo_;
  double frequency_;
};

template <typename cdomaintype>
struct Mu {
  template <typename T1, typename T2, typename T3>
  void evaluate(T1&, T2&, T3& result) const
  {
    // this value is the vaccuum permability mu_0 in N/A^2
    result = 4e-7 * pi;
  }

  struct Traits {
    typedef cdomaintype RangeType;
  };
};

template <typename cdomaintype>
struct Eps {
  template <typename T1, typename T2, typename T3>
  void evaluate(T1&, T2&, T3& result) const
  {
    // this value is the vaccuum permittivity eps_0 in F/m
    result = 8.854187817e-12;
  }

  struct Traits {
    typedef cdomaintype RangeType;
  };
};

class BCTypeParam : public Dune::PDELab::DirichletConstraintsParameters {
public:
  BCTypeParam(const MeshInfo& m) : meshinfo_(m) {}

  template <typename IG>
  bool isDirichlet(const IG& intersection,
                   const Dune::FieldVector<typename IG::ctype,
                                           IG::dimension - 1>&) const {
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
  try {
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
    Dune::shared_ptr<Grid> mygrid(Dune::GmshReader<Grid>::read(
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
    // using CON = Dune::PDELab::NoConstraints;

    using VBE = Dune::PDELab::ISTL::VectorBackend<>;
    using GFS =
      Dune::PDELab::GridFunctionSpace<GV, FiniteElementMap, CON, VBE>;

    GFS gfs(mygrid->leafView(), finiteElementMap);
    using CC = GFS::ConstraintsContainer<double>::Type;
    CC cc;

    Dune::PDELab::constraints(bctype, gfs, cc);

    std::cout << "constrained dofs: " << cc.size() << " of "
              << gfs.globalSize() << std::endl;

    using LocaloperatorS = Dune::PDELab::Electrodynamic_S<Mu<double> >;
    using LocaloperatorT = Dune::PDELab::Electrodynamic_T<Eps<double> >;
    using LocaloperatorFull = Electrodynamic_Full<Eps<double>, Mu<double> >;
    Mu<double> mu;
    Eps<double> eps;
    LocaloperatorS localoperatorS(mu);
    LocaloperatorT localoperatorT(eps);
    LocaloperatorFull localoperatorFull(eps, mu, 3e7, meshinfo);
    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    using GOS = Dune::PDELab::GridOperator<GFS, GFS, LocaloperatorS, MBE,
                                           double, double, double, CC, CC>;
    using GOT = Dune::PDELab::GridOperator<GFS, GFS, LocaloperatorT, MBE,
                                           double, double, double, CC, CC>;
    using GOFull = Dune::PDELab::GridOperator<GFS, GFS, LocaloperatorFull, MBE,
                                              double, double, double, CC, CC>;
    MBE mbe(5); // the edge itself + two adjacent elements * two other edges
    GOS goS(gfs, cc, gfs, cc, localoperatorS, mbe);
    GOT goT(gfs, cc, gfs, cc, localoperatorT, mbe);
    GOFull goFull(gfs, cc, gfs, cc, localoperatorFull, mbe);

    using US = GOS::Traits::Domain;
    using UT = GOT::Traits::Domain;
    using UFull = GOFull::Traits::Domain;
    US uS(gfs, 0);
    UT uT(gfs, 0);
    UFull uFull(gfs, 0);

    GOS::Jacobian jS(goS);
    GOT::Jacobian jT(goT);
    GOFull::Jacobian jFull(goFull);
    jS = 0;
    jT = 0;
    jFull = 0;

    goS.jacobian(uS, jS);
    goT.jacobian(uT, jT);
    goFull.jacobian(uFull, jFull);

    // get rhs:
    GOFull::Range rhs(gfs);
    goFull.residual(uFull, rhs);
    auto matrixS = Dune::PDELab::Backend::native(jS);
    auto matrixT = Dune::PDELab::Backend::native(jT);
    auto matrixFull = Dune::PDELab::Backend::native(jFull);

    Dune::storeMatrixMarket(matrixS,    nameprefix + "_matrixS.mtx");
    Dune::storeMatrixMarket(matrixT,    nameprefix + "_matrixT.mtx");
    Dune::storeMatrixMarket(matrixFull, nameprefix + "_matrixFull.mtx");

    // now solve a problem
    using LS = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
    LS ls;
    using U = GOFull::Traits::Domain;
    using SLP = Dune::PDELab::StationaryLinearProblemSolver<GOFull, LS, U>;
    using DGF = Dune::PDELab::DiscreteGridFunction<GFS, U>;

    // generate some vtk output
    using GV = Grid::LeafGridView;
    const GV& gv = mygrid->leafView();

    auto dofreq = [&](double frequency) {
      LocaloperatorFull localoperatorFull2(eps, mu, frequency, meshinfo);
      GOFull goFull2(gfs, cc, gfs, cc, localoperatorFull2, mbe);
      U u(gfs, 0);
      SLP slp(goFull2, ls, u, 1e-10);
      slp.apply();
      DGF udgf(gfs, u);
      Dune::VTKWriter<GV> myvtkwriter(gv, Dune::VTK::nonconforming);
      myvtkwriter.addVertexData
        (std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >
         (udgf, "solution"));
      myvtkwriter.write(nameprefix + "_" + std::to_string((int)frequency),
                        Dune::VTK::appendedraw);
    };

    // scan through frequencies from 1MHz to 300MHz:
    // for (double frequency = 1e6; frequency < 3e8; frequency *= 1.05)
    //   dofreq(frequency);

    // for the unit test:
    // do a single frequency (the one we have a screenshot for)
    dofreq(30426425);

    return 0;
  }
  catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
