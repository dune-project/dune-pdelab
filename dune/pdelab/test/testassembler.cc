// -*- tab-width: 4; indent-tabs-mode: nil -*-

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <random>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/solvers.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include <dune/pdelab/assembler/assembler.hh>
#include <dune/pdelab/assembler/residualengine.hh>
#include <dune/pdelab/assembler/patternengine.hh>
#include <dune/pdelab/assembler/jacobianengine.hh>
#include <dune/pdelab/assembler/applyjacobianengine.hh>

#include <dune/pdelab/localoperator/adapter.hh>
#include <dune/pdelab/localoperator/l2.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondgnew.hh>

struct LOP
{

  template<typename Context>
  void volumeIntegral(Context& ctx) const
  {
    std::cout << "Volume integral on cell " << ctx.entityIndex()
              << ", #DOFS = " << ctx.test().functionSpace().size()
              << ", volume = " << ctx.domain().embedding().global().volume()
              << ", centroid = " << ctx.centroid()
              << std::endl;
    std::cout << "Cell DOFs:";
    for (auto dof : ctx.argument())
      std::cout << " " << dof;
    std::cout << std::endl;
  }


  template<typename Context>
  void boundaryIntegral(Context& ctx) const
  {
    std::cout << "Boundary integral on intersection " << ctx.domain().index() << " of cell " << ctx.inside().entityIndex()
              << ", volume = " << ctx.domain().embedding().global().volume()
              << ", inside volume = " << ctx.domain().embedding().inside().volume()
              << ", centroid = " << ctx.domain().centroid()
              << std::endl;
  }


};





    template<typename GV, typename RF>
    class ConvectionDiffusionModelProblemSample3
    {
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
      const double lambda = 2.0;

    public:
      typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

      static constexpr bool permeabilityIsConstantPerCell()
      {
        return false;
      }

      //Dimension
      int dimension() {return 2;}

      //! tensor diffusion coefficient
      typename Traits::PermTensorType
      A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        typename Traits::PermTensorType I;
        typename Traits::DomainType xglobal = e.geometry().global(x);
        I[0][0] = exp(lambda*cos(13*xglobal[1]));
        I[0][1] = 0;
        I[1][0] = 0;
        I[1][1] = exp(lambda*cos(17*xglobal[0]));
        //for (std::size_t i=0; i<Traits::dimDomain; i++)
        //  for (std::size_t j=0; j<Traits::dimDomain; j++)
        //  I[i][j] = (i==j) ? 1 : 0;
        return I;
      }

      //! velocity field
      typename Traits::RangeType
      b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        typename Traits::DomainType xglobal = e.geometry().global(x);
        typename Traits::RangeType v;
        v[0] = 1.0*cos(11.0*xglobal[0]);
        v[1] = 1.0*sin(7.0*xglobal[1]);
        return v;
      }

      //! sink term
      typename Traits::RangeFieldType
      c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        typename Traits::DomainType xglobal = e.geometry().global(x);
        return 3.0*cos(23.0*xglobal[0])*sin(29*xglobal[1]);
      }

      //! source term
      typename Traits::RangeFieldType
      f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        typename Traits::DomainType xglobal = e.geometry().global(x);
        return  -2.0*exp(lambda*cos(13.0*xglobal[1])) - 2.0*exp(lambda*cos(17.0*xglobal[0])) + 2.0*xglobal[0]*cos(11.0*xglobal[0]) + 2.0*xglobal[1]*sin(7.0*xglobal[1]) + (xglobal[0]*xglobal[0] + xglobal[1]*xglobal[1])*(-11.0*sin(11.0*xglobal[0]) + 7.0*cos(7.0*xglobal[1])) + 3.0*cos(23.0*xglobal[0])*sin(29.0*xglobal[1]) * (xglobal[0]*xglobal[0] + xglobal[1]*xglobal[1]);
      }

      //Boundary condition type
      BCType
      bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;

      }

      //! automatic distinction between inflow and outflow boundary condition
      BCType
      Inflow_or_Outflow (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        typename Traits::RangeType evalb = b(is.inside(), is.geometryInInside().global(x));
        typename Traits::RangeFieldType flux = evalb * is.unitOuterNormal(x);
        if (flux > 0) {return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow;}
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      }


      //! Dirichlet boundary condition value
      typename Traits::RangeFieldType
      g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        typename Traits::DomainType xglobal = e.geometry().global(x);
        return xglobal[0]*xglobal[0] + xglobal[1]*xglobal[1];

      }

      //! Neumann boundary condition
      typename Traits::RangeFieldType
      j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        typename Traits::DomainType xglobal = is.geometry().global(x);
        typename Traits::RangeType gradient; //This will be A applied to the gradient of u: gradient = (A(x)\nabla u).
        gradient[0] = 2.0*xglobal[0]*exp(lambda*cos(13.0*xglobal[1]));
        gradient[1] = 2.0*xglobal[1]*exp(lambda*cos(17.0*xglobal[0]));
        auto normal = is.unitOuterNormal(x);
        return - normal[0]*gradient[0] - normal[1]*gradient[1] + (1.0*cos(11.0*xglobal[0])*normal[0] + 1.0*sin(7.0*xglobal[1])*normal[1]) * (xglobal[0]*xglobal[0] + xglobal[1]*xglobal[1]);
      }

      //! outflow boundary condition
      typename Traits::RangeFieldType
      o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        typename Traits::DomainType xglobal = is.geometry().global(x);
        typename Traits::RangeType gradient; //This will be A applied to the gradient of u: gradient = (A(x)\nabla u).
        gradient[0] = 2.0*xglobal[0]*exp(lambda*cos(13.0*xglobal[1]));
        gradient[1] = 2.0*xglobal[1]*exp(lambda*cos(17.0*xglobal[0]));
        auto normal = is.unitOuterNormal(x);
        return - normal[0]*gradient[0] - normal[1]*gradient[1] + (1.0*cos(11.0*xglobal[0])*normal[0] + 1.0*sin(7.0*xglobal[1])*normal[1]) * (xglobal[0]*xglobal[0] + xglobal[1]*xglobal[1]);
      }
    };





namespace DG {

  namespace LOP = Dune::PDELab::LocalOperator;

  template<typename RF>
  struct Problem
  {

    using Real = RF;

    template<typename Context>
    class ProblemContext
      : public Context
    {

    public:

      ProblemContext(Context&& ctx)
        : Context(std::move(ctx))
      {}

      using Traits = Dune::PDELab::NewConvectionDiffusionParameterTraits<Context>;

      //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
      static constexpr bool permeabilityIsConstantPerCell()
      {
        return false;
      }

      //! tensor diffusion coefficient
      template<typename P>
      typename Traits::PermeabilityTensor A(const P& p)
      {
        typename Traits::PermeabilityTensor I;
        auto x = Context::global(p);
        using std::exp;
        using std::cos;
        I[0][0] = exp(problem().lambda*cos(13*x[1]));
        I[0][1] = 0;
        I[1][0] = 0;
        I[1][1] = exp(problem().lambda*cos(17*x[0]));
        return I;
      }

      //! velocity field
      template<typename P>
      typename Traits::Velocity b(const P& p)
      {
        using std::sin;
        using std::cos;
        auto x = Context::global(p);
        return {
          cos(11.0*x[0]),
          sin(7.0*x[1])
        };
      }

      //! sink term
      template<typename P>
      LOP::RangeField<Context> c(const P& p)
      {
        auto x = Context::global(p);
        return 3.0*cos(23.0*x[0])*sin(29*x[1]);
      }

      //! source term
      template<typename P>
      LOP::RangeField<Context> f(const P& p)
      {
        auto x = Context::global(p);
        return -2.0*exp(problem().lambda*cos(13.0*x[1])) - 2.0*exp(problem().lambda*cos(17.0*x[0])) + 2.0*x[0]*cos(11.0*x[0]) + 2.0*x[1]*sin(7.0*x[1]) + (x[0]*x[0] + x[1]*x[1])*(-11.0*sin(11.0*x[0]) + 7.0*cos(7.0*x[1])) + 3.0*cos(23.0*x[0])*sin(29.0*x[1]) * (x[0]*x[0] + x[1]*x[1]);
      }

      //! boundary condition type function
      template<typename P>
      typename Traits::BoundaryCondition bctype(const P& p)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
      }

      //! Dirichlet boundary condition value
      template<typename P>
      LOP::RangeField<Context> g(const P& p)
      {
        auto x = Context::global(p);
        return x[0]*x[0] + x[1]*x[1];
      }

      //! Neumann boundary condition
      template<typename P>
      LOP::RangeField<Context> j(const P& p)
      {
        auto x = Context::global(p);
        LOP::Gradient<Context> gradient; //This will be A applied to the gradient of u: gradient = (A(x)\nabla u).
        gradient[0] = 2.0*x[0]*exp(problem().lambda*cos(13.0*x[1]));
        gradient[1] = 2.0*x[1]*exp(problem().lambda*cos(17.0*x[0]));
        auto normal = Context::domain().unitOuterNormal(p);
        //auto evalg = g(p);
        return - normal[0]*gradient[0] - normal[1]*gradient[1] + (1.0*cos(11.0*x[0])*normal[0] + 1.0*sin(7.0*x[1])*normal[1]) * (x[0]*x[0] + x[1]*x[1]);
      }

      //! outflow boundary condition
      template<typename P>
      LOP::RangeField<Context> o(const P& p)
      {
        return 0.0;
      }

      const Problem& problem()
      {
        return Context::engine().localOperator().problem();
      }

    };

    RF lambda;

  };

}

int main(int argc, char** argv)
{
  auto& helper = Dune::MPIHelper::instance(argc, argv);

  if(Dune::MPIHelper::isFake)
    std::cout<< "This is a sequential program." << std::endl;
  else{
    if(helper.rank()==0)
      std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
  }

  using Real = double;

  // Create 2D yasp grid
  constexpr int dim = 2;
  Dune::FieldVector<Real,dim> L(1.0);
  std::bitset<dim> P(false);
  using Grid = Dune::YaspGrid<dim>;

  Grid grid(L,{16,16},P,0);

  // Get GridView
  using GV = Grid::LeafGridView;
  auto gv = grid.leafGridView();

  using ES = Dune::PDELab::AllEntitySet<GV>;
  ES es(gv);

  // Create DG space
  constexpr int degree=1;
  using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,Real,degree,dim>;
  FEM fem;

  using CON = Dune::PDELab::NoConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;
  using GFS = Dune::PDELab::GridFunctionSpace<ES,FEM,CON,VBE>;
  GFS gfs(es,fem);
  gfs.name("u");

  Dune::PDELab::Assembler<ES> assembler(es);

  // LOP lop;

  using Vector = Dune::PDELab::Backend::Vector<GFS,Real>;
  using Matrix = Dune::PDELab::Backend::Matrix<Dune::PDELab::ISTL::BCRSMatrixBackend<>,Vector,Vector>;

  Vector solution(gfs,42.0), residual(gfs,0.0);

  auto rng = std::mt19937_64();
  auto dist = std::uniform_real_distribution<Real>(0.0,1.0);

  for (auto& x : solution)
    x = dist(rng);

  /*
  {
    auto engine = Dune::PDELab::ResidualEngine(solution,residual,lop);
    assembler.assemble(engine);
  }

  {
    Dune::PDELab::L2 l2;
    Dune::PDELab::LocalOperatorAdapter adapter(l2);

    auto engine = Dune::PDELab::ResidualEngine(solution,residual,adapter);
    assembler.assemble(engine);
  }

  for (auto& x : residual)
    std::cout << " " << x;
  std::cout << std::endl;

  residual = 0.0;

  {
    Dune::PDELab::NewL2 l2(2);
    auto engine = Dune::PDELab::ResidualEngine(solution,residual,l2);
    assembler.assemble(engine);
  }

  for (auto& x : residual)
    std::cout << " " << x;
  std::cout << std::endl;

  residual = 0.0;

  {
    using Problem = ConvectionDiffusionModelProblemSample3<ES,Real>;
    auto  problem = Problem{};

    using LOP = Dune::PDELab::ConvectionDiffusionDG<Problem,FEM>;
    auto  lop = LOP{
      problem,
      Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,
      Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,
      1.0,
      1
    };

    Dune::PDELab::LocalOperatorAdapter adapter(lop);

    auto engine = Dune::PDELab::ResidualEngine(solution,residual,adapter);
    assembler.assemble(engine);
  }

  for (auto& x : residual)
    std::cout << " " << x;
  std::cout << std::endl;

  residual = 0.0;

  {
    using Problem = DG::Problem<Real>;
    auto  problem = Problem{2.0};

    using LOP = Dune::PDELab::ConvectionDiffusionDGNew<Problem>;
    auto  lop = LOP{
      problem,
      Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,
      Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,
      1.0,
      1
    };

    auto engine = Dune::PDELab::ResidualEngine(solution,residual,lop);
    assembler.assemble(engine);
  }
  */
  {
    using Problem = DG::Problem<Real>;
    auto  problem = Problem{2.0};

    using LOP = Dune::PDELab::ConvectionDiffusionDGNew<Problem>;
    auto  lop = LOP{
      problem,
      Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,
      Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,
      1.0,
      1
    };

    auto mbe = Dune::PDELab::ISTL::BCRSMatrixBackend(5ul);
    using MBE = decltype(mbe);

    using Dune::PDELab::Backend::native;
    using Dune::PDELab::Backend::Native;

    auto pattern_engine = Dune::PDELab::PatternEngine(solution,residual,lop,mbe);
    auto pattern = assembler.assemble(pattern_engine);
    auto jac = Matrix(Dune::PDELab::Backend::attached_container());
    auto stats = std::vector<MBE::Statistics>();
    Dune::PDELab::ISTL::Impl::allocate_bcrs_matrix(gfs.ordering(),gfs.ordering(),*pattern,native(jac),stats);
    jac.attach(gfs,gfs);

    auto jac_engine = Dune::PDELab::JacobianEngine(solution,jac,lop);
    assembler.assemble(jac_engine);

    residual = 0.0;
    auto apply_engine = Dune::PDELab::ApplyJacobianEngine(solution,solution,residual,lop);
    assembler.assemble(apply_engine);

    residual = 0.0;
    auto res_engine = Dune::PDELab::ResidualEngine(solution,residual,lop);
    assembler.assemble(res_engine);

    auto matrix_operator = Dune::MatrixAdapter<Native<decltype(jac)>,Native<Vector>,Native<Vector>>(native(jac));

    using Preconditioner = Dune::SeqILU<Native<decltype(jac)>,Native<Vector>,Native<Vector>>;
    Preconditioner preconditioner(native(jac),1,.9,false);

    auto solver = Dune::CGSolver<Native<Vector>>(matrix_operator,preconditioner,1e-9,5000,2);

    Vector correction(gfs,0.0);

    Dune::InverseOperatorResult res;
    solver.apply(native(correction),native(residual),res);

    solution -= correction;

    Dune::VTKWriter<GV> vtk_writer(gv);

    Dune::PDELab::addSolutionToVTKWriter(vtk_writer,gfs,solution);
    vtk_writer.write("test");

  }

  return 0;
}
