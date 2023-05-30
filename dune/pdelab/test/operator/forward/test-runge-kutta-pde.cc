#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/operator/forward/runge_kutta.hh>
#include <dune/pdelab/operator/forward/instationary/assembler.hh>
#include <dune/pdelab/operator/local_assembly/archetype.hh>
#include <dune/pdelab/operator/adapter.hh>
#include <dune/pdelab/operator/inverse/istl_adapter.hh>

#include <dune/pdelab/pattern/basis_to_pattern.hh>
#include <dune/pdelab/pattern/sparsity_pattern.hh>
#include <dune/pdelab/pattern/pattern_to_matrix.hh>

#include <dune/pdelab/common/convergence/reason.hh>

#include <dune/pdelab/basis/basis.hh>

#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>

#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk.hh>

#include <gtest/gtest.h>

#include <filesystem>
#include <cmath>
#include <numbers>

using namespace Dune::PDELab::Experimental;

static constexpr std::size_t dim = 1;
static constexpr std::size_t degree = 1;
static constexpr double diffusion = 0.005;

using Duration = double;
using TimePoint = double;

struct FullVolumePattern {

  constexpr static auto localAssembleDoVolume() {
    return std::true_type{};
  }

  void localAssemblePatternVolume(
    const Dune::PDELab::Concept::LocalBasis            auto& ltrial,
    const Dune::PDELab::Concept::LocalBasis            auto& ltest,
                                                       auto& lpattern)
  {
    forEachLeafNode(ltrial.tree(), [&](auto& ltrial_node){
      forEachLeafNode(ltest.tree(), [&](auto& ltest_node){
        for (std::size_t i = 0; i != ltrial_node.size(); ++i)
          for (std::size_t j = 0; j != ltest_node.size(); ++j)
            lpattern.addLink(ltrial_node, i, ltrial_node, j);
      });
    });
  }
};

// define a class with the local integrals for l2
struct LocalL2 {

  // whether a(u,v) volume part needs to be computed
  constexpr static bool localAssembleDoVolume() { return true; }

  constexpr static bool localAssembleIsLinear() { return true; }

  // storage of local test function values`
  std::vector<Dune::FieldVector<double,1>> _v_hats;

  // function that evaluates the bi-linear form a(u,v)
  void localAssembleVolume(                 auto  time_point,
    const Concept::LocalBasis               auto& ltrial,
    const Concept::LocalConstContainer      auto& lcoefficients,
    const Concept::LocalBasis               auto& ltest,
    Concept::LocalMutableContainer          auto& lresidual)
  {
    // obtain geometry of the current grid entity
    const auto& geometry = ltest.element().geometry();

    // obtain the local basis of the current test finite element
    const auto& test_local_basis = ltest.tree().finiteElement().localBasis();

    // loop over quadrature points
    for (auto [position, weight] : quadratureRule(geometry, 2)) {
      // evaluate local test function on current quadrature point for each basis
      test_local_basis.evaluateFunction(position, _v_hats);

      // compute trial function: linear combination of gradient of local trial functions
      double u(.0);
      for (std::size_t dof = 0; dof != ltrial.tree().size(); ++dof) // galerkin: v_hats == u_hats
        u += _v_hats[dof] * lcoefficients(ltrial.tree(), dof);

      // compute integration factor
      auto factor = weight * geometry.integrationElement(position);

      // accumulate residual of u*v for each test degree of freedom
      for (std::size_t dof = 0; dof != ltest.tree().size(); ++dof) {
        auto v = _v_hats[dof];
        lresidual.accumulate(ltest.tree(), dof, u * v * factor);
      }
    }
  }

  void localAssembleJacobianVolumeApply(    auto  time_point,
    const Concept::LocalBasis               auto& ltrial,
    const Concept::LocalConstContainer      auto& llin_point,
    const Concept::LocalConstContainer      auto& lpoint,
    const Concept::LocalBasis               auto& ltest,
    Concept::LocalMutableContainer          auto& lresidual)
  {
    static_assert(localAssembleIsLinear());
    localAssembleVolume(time_point, ltrial, lpoint, ltest, lresidual);
  }

  void localAssembleJacobianVolume(            auto  time_point,
    const Concept::LocalBasis                  auto& ltrial,
    const Concept::LocalConstContainer         auto& llin_point,
    const Concept::LocalBasis                  auto& ltest,
    Concept::LocalMutableMatrix                auto& ljacobian)
  {
    const auto& geometry = ltest.element().geometry();

    // obtain the local basis of the current test finite element
    const auto& test_local_basis = ltest.tree().finiteElement().localBasis();

    // loop over quadrature points
    for (auto [position, weight] : quadratureRule(geometry, 2)) {
      // evaluate local test function on current quadrature point for each basis
      test_local_basis.evaluateFunction(position, _v_hats);

      // compute integration factor
      auto factor = weight * geometry.integrationElement(position);

      // integrate mass matrix
      for (std::size_t dof_i = 0; dof_i != ltest.tree().size(); ++dof_i)
        for (std::size_t dof_j = 0; dof_j != ltrial.tree().size(); ++dof_j)
          ljacobian.accumulate(ltest.tree(), dof_i, ltrial.tree(), dof_j, _v_hats[dof_i] * _v_hats[dof_j] * factor);
    }
  }
};

// define a class with the local integrals for the poisson problem
struct LocalPoisson {

  constexpr static auto localAssembleIsLinear() { return true; }

  // whether a(u,v) volume part needs to be computed
  constexpr static auto localAssembleDoVolume()   { return true; }

  // storage of jacobian and gradients of local test functions
  std::vector<Dune::FieldVector<double,1>>      _v_hats;
  std::vector<Dune::FieldMatrix<double,1,dim>>  _jac_v_hats;
  std::vector<Dune::FieldVector<double,dim>>    _grad_v_hats;

  // function that evaluates the bi-linear form a(u,v)
  void localAssembleVolume(                 auto  time_point,
    const Concept::LocalBasis               auto& ltrial,
    const Concept::LocalConstContainer      auto& lcoefficients,
    const Concept::LocalBasis               auto& ltest,
    Concept::LocalMutableContainer          auto& lresidual)
  {
    // obtain geometry of the current grid entity
    const auto& geometry = ltest.element().geometry();

    // obtain the local basis of the current test finite element
    const auto& test_local_basis = ltest.tree().finiteElement().localBasis();

    // loop over quadrature points
    for (auto [position, weight] : quadratureRule(geometry, 2)) {
      // evaluate jacobian of local test function on current quadrature point for each basis
      test_local_basis.evaluateJacobian(position, _jac_v_hats);

      // compute gradient of local test function for each basis
      _grad_v_hats.resize(ltest.tree().size());
      for (std::size_t dof = 0; dof != ltest.tree().size(); ++dof)
        _grad_v_hats[dof] = (_jac_v_hats[dof] * geometry.jacobianInverse(position))[0];

      // compute gradient of trial function: linear combination of gradient of local trial functions
      Dune::FieldVector<double, dim> grad_u(.0);
      for (std::size_t dof = 0; dof != ltrial.tree().size(); ++dof) // galerkin: _grad_v_hats == grad_u_hats
        grad_u += _grad_v_hats[dof] * lcoefficients(ltrial.tree(), dof);

      // compute integration factor
      auto factor = weight * geometry.integrationElement(position);

      // accumulate residual of bi-linear form a(u,v) on each local test function v
      for (std::size_t dof = 0; dof != ltest.tree().size(); ++dof) {
        auto grad_v = _grad_v_hats[dof];
        lresidual.accumulate(ltest.tree(), dof, diffusion * dot(grad_u, grad_v) * factor);
      }
    }
  }

  void localAssembleJacobianVolumeApply(    auto  time_point,
    const Concept::LocalBasis               auto& ltrial,
    const Concept::LocalConstContainer      auto& llin_point,
    const Concept::LocalConstContainer      auto& lpoint,
    const Concept::LocalBasis               auto& ltest,
    Concept::LocalMutableContainer          auto& lresidual)
  {
    static_assert(localAssembleIsLinear());
    localAssembleVolume(time_point, ltrial, lpoint, ltest, lresidual);
  }

  void localAssembleJacobianVolume(            auto  time_point,
    const Concept::LocalBasis                  auto& ltrial,
    const Concept::LocalConstContainer         auto& llin_point,
    const Concept::LocalBasis                  auto& ltest,
    Concept::LocalMutableMatrix                auto& ljacobian)
  {
    const auto& geometry = ltest.element().geometry();

    // obtain the local basis of the current test finite element
    const auto& test_local_basis = ltest.tree().finiteElement().localBasis();

    // loop over quadrature points
    for (auto [position, weight] : quadratureRule(geometry, 2)) {
      // evaluate local test function on current quadrature point for each basis
      test_local_basis.evaluateJacobian(position, _jac_v_hats);
      // compute gradient of local test function for each basis
      _grad_v_hats.resize(ltest.tree().size());
      for (std::size_t dof = 0; dof != ltest.tree().size(); ++dof)
        _grad_v_hats[dof] = (_jac_v_hats[dof] * geometry.jacobianInverse(position))[0];

      // compute integration factor
      auto factor = weight * geometry.integrationElement(position);

      // integrate mass matrix
      for (std::size_t dof_i = 0; dof_i != ltest.tree().size(); ++dof_i)
        for (std::size_t dof_j = 0; dof_j != ltrial.tree().size(); ++dof_j)
          ljacobian.accumulate(ltest.tree(), dof_i, ltrial.tree(), dof_j, diffusion * dot(_grad_v_hats[dof_i], _grad_v_hats[dof_j]) * factor);
    }
  }
};

void interpolate(const Concept::Basis auto& basis, auto& coefficients, const auto& func) {
  // local coefficients
  std::vector<double> lcoefficients;
  // loop once over the grid and interpolate
  auto lbasis = basis.localView();
  for (const auto& entity : elements(basis.entitySet())) {
    // bind local function basis to element
    lbasis.bind(entity);
    // clean up local container
    lcoefficients.assign(lbasis.size(), 0.);
    // define function f in terms of global coordinate
    auto lfunc = [&](auto xlocal){
      return func(entity.geometry().global(xlocal));
    };
    // interpolate local function into local coefficients
    lbasis.tree().finiteElement().localInterpolation().interpolate(lfunc, lcoefficients);
    // assign local coefficients to global coefficients
    for (std::size_t dof = 0; dof != lbasis.tree().size(); ++dof)
      Dune::PDELab::containerEntry(coefficients, lbasis.tree().index(dof)) = lcoefficients[dof];
    lbasis.unbind();
  }
}


template<Concept::Basis Basis, class Container>
class ScalarDiscreteFunction {

  static_assert(Concept::LeafTreeNode<typename Basis::LocalView::Tree>);

  class LocalFunction
  {
  public:
    LocalFunction(const Basis& space, const Container& container)
      : _lspace{space.localView()}
      , _lcontainer{space, std::cref(container)}
    {}


    void bind(const auto& element)
    {
      _lspace.bind(element);
      _lcontainer.load(_lspace);
    }

    auto operator()(const auto& xlocal) const
    {
      const auto& node = _lspace.tree();
      node.finiteElement().localBasis().evaluateFunction(xlocal, _range);
      auto value = _lcontainer(node, 0) * _range[0];
      for (std::size_t dof = 1; dof < node.size(); ++dof)
        value += _lcontainer(node, dof) * _range[dof];
      return value;
    }

    void unbind()
    {
      _lcontainer.clear(_lspace);
      _lspace.unbind();
    }

  private:
    using FiniteElement = typename Basis::LocalView::Tree::FiniteElement;
    using Range = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;

    mutable typename Basis::LocalView _lspace;
    mutable LocalContainerBuffer<Basis,const Container> _lcontainer;
    mutable std::vector<Range> _range;
  };

public:
  friend LocalFunction localFunction(const ScalarDiscreteFunction& discrete_function) {
    return LocalFunction{discrete_function._space, discrete_function._container.get()};
  }


  ScalarDiscreteFunction(const Basis& space, const Container& container)
    : _space{space}
    , _container{container}
  {}

private:
  Basis _space;
  std::reference_wrapper<const Container> _container;
};




TEST(TestRungeKutta, TestRungeKuttaExplicitPDE) {

  // define the parameters of the domain
  Dune::FieldVector<double,dim> grid_extent(1.);     // extent of the domain per each direction
  std::array<int,dim> grid_cells;
  std::fill(begin(grid_cells), end(grid_cells), 1);  // number of grid cells per direction

  // create a grid using yasp grid
  using Grid = Dune::YaspGrid<dim>;
  auto grid = std::make_unique<Grid>(grid_extent, grid_cells);

  // refine the square n times
  grid->globalRefine(7);

  // use leaf grid view of the grid for our entity set of our function spaces
  using GridView = typename Grid::LeafGridView;
  using EntitySet = GridView;
  EntitySet entity_set{grid->leafGridView()};

  // select finite element map: grid entity -> Q1 finite element
  using FEM = Dune::PDELab::QkLocalFiniteElementMap<EntitySet, double, double, degree>;
  auto q1_fem = std::make_shared<FEM>(entity_set);

  // select a merging strategy betwee degrees of freedom (DoFs)
  // flatByEntity: DoFs attached to a grid entity are groupped together and are layed out one after the other (flat merging)
  auto strategy = Strategy::flatByEntity(entity_set);

  // construct tree of (unordered) function spaces, in this case, we only have one node
  auto leaf_pre_space = makePreBasis(strategy, q1_fem);

  // convert the tree of (unordered) function spaces to an ordered function space
  auto trial = makeBasis(entity_set, leaf_pre_space);
  using Basis = decltype(trial);

  // since we have a galerkin method, test and trial are the same space
  Basis test = trial;

  // make container for storing the coefficients of the finite element problem  (e.g. std::vector<double>)
  using Domain = Dune::BlockVector<double>;
  using Range = Dune::BlockVector<double>;

  using RungeKuttaDomain = Dune::BlockVector<Domain>;
  using RungeKuttaRange = Dune::BlockVector<Range>;

  using RungeKuttaJacobian = Dune::DynamicMatrix<Dune::BCRSMatrix<double>>;

  auto jac_resize = [&](const auto& op, RungeKuttaJacobian& jac){
    LeafSparsePattern<Basis, Basis> pattern{test, {}, trial, {}, 5};
    jac.resize(1,1); // runge kutta stages/sub-steps

    basisToPattern(FullVolumePattern{}, pattern);
    pattern.sort();
    patternToMatrix(pattern, jac[0][0]);
    jac[0][0] = 0.;
  };

  std::shared_ptr<Operator<RungeKuttaRange,RungeKuttaDomain>> instationary;

  bool matrix_free = true;

  // jacobian type
  if (matrix_free)
    instationary = makeInstationaryMatrixFreeAssembler<RungeKuttaRange,RungeKuttaDomain>(test, trial, LocalL2{}, LocalPoisson{});
  else
    instationary = makeInstationaryMatrixBasedAssembler<RungeKuttaRange,RungeKuttaDomain,RungeKuttaJacobian>(test, trial, LocalL2{}, LocalPoisson{}, jac_resize);

  // set up solver the linear problem
  // for explicit or semi-implicit methods, we need to solve the diagonal block in the runge-kutta system
  // for fully-impliciy methids, we need to solve the whole system
  auto linear_apply = [&](Operator<RungeKuttaRange, RungeKuttaDomain>& inverse, RungeKuttaRange& b, RungeKuttaDomain& x) -> ErrorCondition {
    static_assert(std::is_same_v<RungeKuttaDomain,RungeKuttaRange>);
    // obtine the forward operator as a differentiable fuction (notice that the foward operator may have a partial set of coefficients from the shu-osher tableau)
    auto& forward = inverse.get<Operator<RungeKuttaDomain, RungeKuttaRange>>("forward");

    // define generic scalar product and a simple preconditioner
    auto scalar_product_op = std::make_shared<Dune::ScalarProduct<RungeKuttaDomain>>();
    auto pre_op = std::make_shared<Dune::Richardson<RungeKuttaDomain,RungeKuttaRange>>(0.5);

    // set derivative (jacobian) of the mass-stiffness problem as the matrix vector operation
    auto dx = forward.derivative(x);
    assert(dx);
    auto istl_derivative_op = std::make_shared<ISTL::LinearAdapter<RungeKuttaDomain,RungeKuttaRange>>(dx);

    // calculate the right hand side of the mass-stiffness problem
    forward.apply(x,b).or_throw();

    // use the conjugate gradient method
    auto solver = Dune::CGSolver<RungeKuttaDomain>{istl_derivative_op, scalar_product_op, pre_op, /*reduction*/ 1e-15, /*max_it*/ 400, /*verbosity*/ 1};

    // compute correction
    auto z = x;
    Dune::InverseOperatorResult res;
    solver.apply(z, b, res);

    // update current solution with correction
    Dune::PDELab::axpy(std::execution::par_unseq,x,-1.0,z);

    if (res.converged)
      return ErrorCondition{};
    else
      return make_error_condition(Convergence::Reason::DivergedNull);
  };

  // the runge kutta operator needs an operator that can invert the mass-stiffness problem
  std::shared_ptr<Operator<RungeKuttaRange,RungeKuttaDomain>> inverse = std::make_shared<OperatorAdapter<RungeKuttaRange, RungeKuttaDomain>>(linear_apply);

  RungeKutta<RungeKuttaDomain,RungeKuttaRange> runge_kutta;
  TimePoint time{0.5};
  runge_kutta["initial_residual"] = (Range(test.dimension()) = 0.);
  runge_kutta["inverse"] = inverse;
  runge_kutta["inverse.forward"] = std::weak_ptr(instationary);
  runge_kutta["inverse.forward.time_point"] = std::ref(time);
  runge_kutta["inverse.forward.duration"] = Duration{0.02};
  runge_kutta["instationary_coefficients"] = InstationaryCoefficients(Dune::PDELab::Alexander2Parameter<double>{});

  Domain x0(trial.dimension()), x1;
  x0 = 0.;

  auto ref = [&time](auto pos){
    double dist2{0.}; // distance to the center of the domain
    for (auto x : pos) dist2 += (x-0.5)*(x-0.5);
    return std::exp(-dist2/(4*time*diffusion)) / (2*std::sqrt(std::numbers::pi*time*diffusion));
  };
  auto uref = Dune::Functions::makeAnalyticGridViewFunction(ref, trial.entitySet());

  // set initial condition
  interpolate(trial, x0, ref);

  // create a discrete function out of function space and coefficietns
  ScalarDiscreteFunction u{trial, x0};
  // set up VTK writer
  std::filesystem::create_directory("vtk");
  Dune::VTKSequenceWriter vtkwriter{entity_set, "poisson", "vtk", "vtk"};
  vtkwriter.vtkWriter()->addVertexData(u, {"u", Dune::VTK::FieldInfo::Type::scalar, 1});
  vtkwriter.vtkWriter()->addVertexData(uref, {"u-ref", Dune::VTK::FieldInfo::Type::scalar, 1});

  // write initial condition
  vtkwriter.write(time, Dune::VTK::appendedraw);

  while(Dune::FloatCmp::lt(time, TimePoint{0.6})) {
    runge_kutta.apply(x1 = x0, x0).or_throw();
    // write results at this time step
    vtkwriter.write(time, Dune::VTK::appendedraw);

    // auto lu = localFunction(u);
    // for (auto&& entity : elements(test.entitySet())) {
    //   lu.bind(entity);
    //   for (auto [position, weight] : quadratureRule(entity.geometry(), 2))
    //     EXPECT_NEAR(lu(position), ref(position), 1e-4);
    // }

  }
  std::cout << runge_kutta << std::endl;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  Dune::MPIHelper::instance(argc, argv);
  return RUN_ALL_TESTS();
}
