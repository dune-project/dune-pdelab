// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/pkqkfem.hh>
#include <dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/pattern.hh>

using namespace Dune;

template <typename T>
class LinearPrimalPlasticity
    : public PDELab::FullVolumePattern,
      public PDELab::LocalOperatorDefaultFlags,
      public PDELab::InstationaryLocalOperatorDefaultMethods<
          typename T::Domain> {
 public:
  typedef T ParameterType;

  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doLambdaVolume = true };
  enum { doLambdaBoundary = true };

  template <typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& x,
                       const LFSV& lfsv, M& mat) const {}

  // volume integral depending on test and ansatz functions
  template <typename EG, typename LFSU_HAT, typename X, typename LFSV,
            typename R>
  void alpha_volume(const EG& eg, const LFSU_HAT& lfsu_hat, const X& x,
                    const LFSV& lfsv, R& r) const {}

  // volume integral depending only on test functions
  template <typename EG, typename LFSV_HAT, typename R>
  void lambda_volume(const EG& eg, const LFSV_HAT& lfsv_hat, R& r) const {}

  /** \brief Boundary loads */
  template <typename IG, typename LFSV_HAT, typename R>
  void lambda_boundary(const IG& ig, const LFSV_HAT& lfsv_hat, R& r) const {}
};

/** \brief Stores the material parameters of a small-strain plastic material
 * with kinematic hardening
 */
template <class GridView>
class StVenantKirchhoffParameters {
  static const int dim = GridView::dimension;

 public:
  typedef FieldVector<double, dim> Domain;
};

int main(int argc, char* argv[]) try {
  // Init MPI, if present
  MPIHelper::instance(argc, argv);

  ///////////////////////////////
  //   Generate the grid
  ///////////////////////////////

  // The grid dimension
  const int dim = 2;

  // typedef UGGrid<dim> GridType;
  typedef YaspGrid<dim> GridType;
  shared_ptr<GridType> grid;
  FieldVector<double, dim> lower(0), upper(10);
  std::array<unsigned int, dim> elements = {2, 2};

  grid =
      StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

  // The number of values needed to store one plastic strain tensor
  const size_t nPlasticStrainComponents = (dim == 2) ? 2 : 5;

  ////////////////////////////////////////////////////////////////
  //  Construct grid function space
  ////////////////////////////////////////////////////////////////

  // set up vector-valued finite element space
  typedef PDELab::PkQkLocalFiniteElementMap<GridView::ctype, double, dim, 1>
      FEM;
  FEM fem;

  typedef PDELab::VectorGridFunctionSpace<
    GridView, FEM, dim, PDELab::ISTL::VectorBackend<>,
      PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none,1>, PDELab::ConformingDirichletConstraints,
      PDELab::EntityBlockedOrderingTag, PDELab::DefaultLeafOrderingTag>
      DisplacementGFS;

  DisplacementGFS displacementGFS(gridView, fem);
  displacementGFS.name("displacement");

  typedef PDELab::VectorGridFunctionSpace<
      GridView, FEM, nPlasticStrainComponents, PDELab::ISTL::VectorBackend<>,
      PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none,1>, PDELab::ConformingDirichletConstraints,
      PDELab::EntityBlockedOrderingTag, PDELab::DefaultLeafOrderingTag>
      PlasticStrainGFS;

  PlasticStrainGFS plasticStrainGFS(gridView, fem);
  plasticStrainGFS.name("plastic strain");

  typedef PDELab::ISTL::VectorBackend<PDELab::ISTL::Blocking::fixed> VBE;
  typedef Dune::PDELab::CompositeGridFunctionSpace<
      VBE, PDELab::EntityBlockedOrderingTag, DisplacementGFS, PlasticStrainGFS>
      CompositeGFS;

  // sets up composition of spaces
  CompositeGFS gfscomp(displacementGFS, plasticStrainGFS);

  ///////////////////////////////////////////////////////////////
  // make grid operator
  ///////////////////////////////////////////////////////////////

  LinearPrimalPlasticity<StVenantKirchhoffParameters<GridView> > lop;

  typedef PDELab::ISTL::BCRSMatrixBackend<> MatrixBackend;
  MatrixBackend mbe(
      27);  // estimate for the average number of matrix entries per row

  // set up linear operator acting on the FEM space
  typedef PDELab::GridOperator<
      CompositeGFS, CompositeGFS,
      LinearPrimalPlasticity<StVenantKirchhoffParameters<GridView> >,
      MatrixBackend, double, double, double>
      GridOperator;

  GridOperator gridOperator(gfscomp, gfscomp, lop,
                            mbe);  // Assertion failure here!

  typedef GridOperator::Traits::Domain V;
  typedef GridOperator::Jacobian M;
  using MatrixType = M::Container;
  using VectorType = V::Container;

  std::cout << "Matrix type: " << className<MatrixType>() << std::endl;
  std::cout << "Vector type: " << className<VectorType>() << std::endl;

} catch (Exception e) {
  std::cout << e << std::endl;
}
