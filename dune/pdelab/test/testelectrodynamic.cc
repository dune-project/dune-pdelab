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

// #include <dune/stuff/common/disable_warnings.hh>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh>         // We use exceptions

#include <dune/alugrid/grid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/utility/vertexorderfactory.hh>
#include <dune/pdelab/finiteelementmap/edges0.5fem.hh>

// for NoConstraints
#include <dune/pdelab/constraints/noconstraints.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>

#include <dune/istl/matrixmarket.hh>
// #include <dune/stuff/common/disable_warnings.hh>
// #include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istl.hh>
// #include <dune/pdelab/backend/istlmatrixbackend.hh>
// #include <dune/pdelab/backend/istlsolverbackend.hh>

#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/constraints/conforming.hh>
// #include <dune/pdelab/backend/seqistlsolverbackend.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/localfunctions/common/interfaceswitch.hh>

// #include <dune/stuff/common/reenable_warnings.hh>

// #include "localoperator.hh"
// #include "meshinfo.hh"

#include <boost/lexical_cast.hpp>

// meshinfo.hh

#include <vector>

class MeshInfo {
public:
  std::vector<int> boundary_id_to_physical_entity;
  std::vector<int> element_index_to_physical_entity;
};

// localoperator.hh

#include <dune/pdelab/localoperator/electrodynamic.hh>

// #include "meshinfo.hh"

const double pi = 3.14159265358979323846264338328;

namespace Dune {
namespace PDELab {

template <typename Eps, typename Mu>
class Electrodynamic_Full
    : public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public JacobianBasedAlphaVolume<Electrodynamic_Full<Eps, Mu>> {
public:
  enum {
    doPatternVolume = true
  };
  enum {
    doAlphaVolume = true
  };
  enum {
    doLambdaBoundary = true
  };

  Electrodynamic_Full(const Eps& eps, const Mu& mu, double frequency,
                      const MeshInfo& m, int qorder = 2)
    : edyn_t_lop_(eps, qorder)
    , edyn_s_lop_(mu, qorder)
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
    M weightedview = mat.weightedAccumulationView(-omega_squared);

    // and assemble mass part:
    edyn_t_lop_.jacobian_volume(eg, lfsu, x, lfsv, weightedview);
  }

  template <typename IG, typename LFSV, typename R>
  void lambda_boundary(const IG& ig, const LFSV& lfsv_s, R& r_s) const {
    typedef typename LFSV::Traits::FiniteElementType::Traits::Basis::Traits
    BasisTraits;

    typedef typename BasisTraits::DomainField DF;
    typedef typename BasisTraits::DomainLocal DomainLocal;

    typedef typename BasisTraits::RangeField RF;
    typedef typename BasisTraits::Range Range;

    int boundarySegmentIndex = ig.intersection().boundarySegmentIndex();
    int physical_id =
        meshinfo_.boundary_id_to_physical_entity[boundarySegmentIndex];

    if (physical_id == 3) {
      // get a quad rule for element:
      auto gt = ig.geometry().type();
      const auto& rule =
          Dune::QuadratureRules<DF, IG::Geometry::mydimension>::rule(gt, 2);
      RF result(0.);
      Range tangentialvalue(0);
      tangentialvalue[1] = 1.;

      typename IG::LocalGeometry geo(ig.geometryInInside());
      for (auto& quadraturepoint : rule) {
        std::vector<Range> phi(lfsv_s.size());
        // map quadraturepoint to tet, then evaluate
        lfsv_s.finiteElement().basis().evaluateFunction(
            geo.global(quadraturepoint.position()), phi);

        for (unsigned i = 0; i < lfsv_s.size(); ++i) {
          r_s.accumulate(lfsv_s, i, frequency_ *
                                        (phi[i] * tangentialvalue) *
                                        quadraturepoint.weight());
        }
      }
    }
  }

private:
  Electrodynamic_S<Mu> edyn_s_lop_;
  Electrodynamic_T<Eps> edyn_t_lop_;

  const MeshInfo& meshinfo_;
  double frequency_;
};
}
};     // end namespace

// circ_in_rect_2d.cc

template <typename cdomaintype>
struct Mu {
  template <typename T1, typename T2, typename T3>
  void evaluate(T1&, T2&, T3& result) const {
    result = 4e-7 * pi;
  }

  struct Traits {
    typedef cdomaintype RangeType;
  };
};

template <typename cdomaintype>
struct Eps {
  template <typename T1, typename T2, typename T3>
  void evaluate(T1&, T2&, T3& result) const {
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

    // create the grid
    const int dim = 2;
    typedef Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>
    GridType;

    MeshInfo meshinfo;
    Dune::shared_ptr<GridType> mygrid(Dune::GmshReader<GridType>::read(
        "mesh.msh", meshinfo.boundary_id_to_physical_entity,
        meshinfo.element_index_to_physical_entity));
    Dune::gridinfo(*mygrid);

    BCTypeParam bctype(meshinfo);

    // now we need a vertex order
    typedef GridType::Codim<0>::Geometry Geometry;
    typedef Dune::VertexOrderByIdFactory<GridType::GlobalIdSet>
    VertexOrderFactory;
    VertexOrderFactory vertexOrderFactory(mygrid->globalIdSet());
    typedef double RF;
    typedef Dune::PDELab::EdgeS0_5FiniteElementMap<
        Geometry, VertexOrderFactory, RF> FiniteElementMap;
    FiniteElementMap finiteElementMap(vertexOrderFactory);

    typedef GridType::LeafGridView GV;
    typedef Dune::PDELab::ConformingDirichletConstraints CON;
    // typedef Dune::PDELab::NoConstraints CON;

    typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
    typedef Dune::PDELab::GridFunctionSpace<GV, FiniteElementMap, CON, VBE>
    GFS;

    GFS gfs(mygrid->leafView(), finiteElementMap);
    typedef GFS::ConstraintsContainer<double>::Type CC;
    CC cc;

    Dune::PDELab::constraints(bctype, gfs, cc);

    std::cout << "constrained dofs: " << cc.size() << " of "
              << gfs.globalSize() << std::endl;

    typedef Dune::PDELab::Electrodynamic_S<Mu<double>> LocaloperatorS;
    typedef Dune::PDELab::Electrodynamic_T<Eps<double>> LocaloperatorT;
    typedef Dune::PDELab::Electrodynamic_Full<Eps<double>, Mu<double>>
    LocaloperatorFull;
    Mu<double> mu;
    Eps<double> eps;
    LocaloperatorS localoperatorS(mu);
    LocaloperatorT localoperatorT(eps);
    LocaloperatorFull localoperatorFull(eps, mu, 3e7, meshinfo);
    typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
    typedef Dune::PDELab::GridOperator<GFS, GFS, LocaloperatorS, MBE,
                                       double, double, double, CC, CC> GOS;
    typedef Dune::PDELab::GridOperator<GFS, GFS, LocaloperatorT, MBE,
                                       double, double, double, CC, CC> GOT;
    typedef Dune::PDELab::GridOperator<GFS, GFS, LocaloperatorFull, MBE,
                                       double, double, double, CC,
                                       CC> GOFull;
    MBE mbe(4);  // the diagonal plus 3 neighbor elements
    GOS goS(gfs, cc, gfs, cc, localoperatorS, mbe);
    GOT goT(gfs, cc, gfs, cc, localoperatorT, mbe);
    GOFull goFull(gfs, cc, gfs, cc, localoperatorFull, mbe);

    typedef GOS::Traits::Domain US;
    typedef GOT::Traits::Domain UT;
    typedef GOFull::Traits::Domain UFull;
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

    Dune::storeMatrixMarket(matrixS, "matrixS.mtx");
    Dune::storeMatrixMarket(matrixT, "matrixT.mtx");
    Dune::storeMatrixMarket(matrixFull, "matrixFull.mtx");

    // now solve a problem
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
    LS ls;
    typedef GOFull::Traits::Domain U;
    typedef Dune::PDELab::StationaryLinearProblemSolver<GOFull, LS, U> SLP;
    typedef Dune::PDELab::DiscreteGridFunction<GFS, U> DGF;

    // generate some vtk output
    typedef GridType::LeafGridView GV;
    const GV& gv = mygrid->leafView();

    for (double frequency = 1e6; frequency < 3e8; frequency *= 1.05) {
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
      myvtkwriter.write(
          "vtkout_" + boost::lexical_cast<std::string>((int)frequency),
          Dune::VTK::appendedraw);
    }

    return 0;
  }
  catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
