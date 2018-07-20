// -*- tab-width: 4; indent-tabs-mode: nil -*-

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <random>

#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/common/benchmarkhelper.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>

#include <dune/pdelab/assembler/assembler.hh>
#include <dune/pdelab/assembler/residualengine.hh>
#include <dune/pdelab/assembler/patternengine.hh>
#include <dune/pdelab/assembler/jacobianengine.hh>

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


template<typename GV, typename FEM>
void benchmark(const Dune::ParameterTree& params, const GV& gv, const FEM& fem)
{

  constexpr int dim = GV::dimension;
  using Real = typename GV::ctype;

  using ES = Dune::PDELab::AllEntitySet<GV>;
  ES es(gv);

  using CON = Dune::PDELab::NoConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;
  using GFS = Dune::PDELab::GridFunctionSpace<ES,FEM,CON,VBE>;
  GFS gfs(es,fem);
  gfs.ordering();

  auto mbe = Dune::PDELab::ISTL::BCRSMatrixBackend(5ul);
  using MBE = decltype(mbe);

  using Vector = Dune::PDELab::Backend::Vector<GFS,Real>;
  using Matrix = Dune::PDELab::Backend::Matrix<MBE,Vector,Vector>;

  Vector solution(gfs,0.0), residual(gfs,0.0);

  auto rng = std::mt19937_64();
  auto dist = std::uniform_real_distribution<Real>(0.0,1.0);

  for (auto& x : solution)
    x = dist(rng);

  using OldProblem   = ConvectionDiffusionModelProblemSample3<ES,Real>;
  auto  old_problem = OldProblem{};

  using Problem = DG::Problem<Real>;
  auto  problem = Problem{2.0};

  using OldLocalOperator   = Dune::PDELab::ConvectionDiffusionDG<OldProblem,FEM>;
  auto  old_local_operator = OldLocalOperator{
    old_problem,
    Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,
    Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,
    1.0,
    1
  };

  auto wrapped_local_operator = Dune::PDELab::LocalOperatorAdapter(old_local_operator);

  using LocalOperator  = Dune::PDELab::ConvectionDiffusionDGNew<Problem>;
  auto  local_operator = LocalOperator{
    problem,
    Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,
    Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,
    1.0,
    1
  };

  Dune::PDELab::GridOperator<GFS,GFS,OldLocalOperator,MBE,Real,Real,Real> grid_operator(gfs,gfs,old_local_operator,mbe);

  Dune::PDELab::Assembler<ES> assembler(es);

  Vector r_old(gfs,0.0), r_wrapped(gfs,0.0), r_new(gfs,0.0);
  Matrix j_old(grid_operator,0.0), j_wrapped(grid_operator,0.0), j_new(grid_operator,0.0);

  auto runs = params.get<std::size_t>("benchmark.runs");
  {

    grid_operator.residual(solution,r_old);
    grid_operator.jacobian(solution,j_old);

    auto helper = Dune::PDELab::BenchmarkHelper("old operator",runs);

    for (std::size_t i = 0 ; i < runs ; ++i)
      {
        helper.start_run(std::cout);
        helper.start("reset residual",std::cout);
        residual = 0.0;
        helper.end("reset residual",std::cout);

        helper.start("residual",std::cout);
        grid_operator.residual(solution,residual);
        helper.end("residual",std::cout);

        helper.start("matrix setup",std::cout);
        auto jacobian = Matrix(grid_operator,0.0);
        helper.end("matrix setup",std::cout);

        helper.start("jacobian",std::cout);
        grid_operator.jacobian(solution,jacobian);
        helper.end("jacobian",std::cout);

        helper.end_run(std::cout);
      }

    helper.print(std::cout);
  }

  {

    {
      auto residual_engine = Dune::PDELab::ResidualEngine(solution,r_wrapped,wrapped_local_operator);
      assembler.assemble(residual_engine);
      auto jacobian_engine = Dune::PDELab::JacobianEngine(solution,j_wrapped,wrapped_local_operator);
      assembler.assemble(jacobian_engine);
    }

    auto helper = Dune::PDELab::BenchmarkHelper("wrapped operator",runs);

    for (std::size_t i = 0 ; i < runs ; ++i)
      {
        helper.start_run(std::cout);
        helper.start("reset residual",std::cout);
        residual = 0.0;
        helper.end("reset residual",std::cout);

        helper.start("residual",std::cout);
        {
          auto engine = Dune::PDELab::ResidualEngine(solution,residual,wrapped_local_operator);
          assembler.assemble(engine);
        }
        helper.end("residual",std::cout);

        helper.start("matrix setup",std::cout);
        auto jacobian = Matrix(Dune::PDELab::Backend::attached_container());
        {
          using Dune::PDELab::Backend::native;
          auto engine = Dune::PDELab::PatternEngine(solution,residual,wrapped_local_operator,mbe);
          auto pattern = assembler.assemble(engine);
          auto stats = std::vector<typename MBE::Statistics>();
          Dune::PDELab::ISTL::Impl::allocate_bcrs_matrix(gfs.ordering(),gfs.ordering(),*pattern,native(jacobian),stats);
          jacobian.attach(gfs,gfs);
          jacobian = 0.0;
        }
        helper.end("matrix setup",std::cout);

        helper.start("jacobian",std::cout);
        {
          auto engine = Dune::PDELab::JacobianEngine(solution,jacobian,wrapped_local_operator);
          assembler.assemble(engine);
        }
        helper.end("jacobian",std::cout);

        helper.end_run(std::cout);
      }

    helper.print(std::cout);
  }

  {

    {
      auto residual_engine = Dune::PDELab::ResidualEngine(solution,r_new,local_operator);
      assembler.assemble(residual_engine);
      auto jacobian_engine = Dune::PDELab::JacobianEngine(solution,j_new,local_operator);
      assembler.assemble(jacobian_engine);
    }

    auto helper = Dune::PDELab::BenchmarkHelper("new operator",runs);

    for (std::size_t i = 0 ; i < runs ; ++i)
      {
        helper.start_run(std::cout);
        helper.start("reset residual",std::cout);
        residual = 0.0;
        helper.end("reset residual",std::cout);

        helper.start("residual",std::cout);
        {
          auto engine = Dune::PDELab::ResidualEngine(solution,residual,local_operator);
          assembler.assemble(engine);
        }
        helper.end("residual",std::cout);

        helper.start("matrix setup",std::cout);
        auto jacobian = Matrix(Dune::PDELab::Backend::attached_container());
        {
          using Dune::PDELab::Backend::native;
          auto engine = Dune::PDELab::PatternEngine(solution,residual,local_operator,mbe);
          auto pattern = assembler.assemble(engine);
          auto stats = std::vector<typename MBE::Statistics>();
          Dune::PDELab::ISTL::Impl::allocate_bcrs_matrix(gfs.ordering(),gfs.ordering(),*pattern,native(jacobian),stats);
          jacobian.attach(gfs,gfs);
          jacobian = 0.0;
        }
        helper.end("matrix setup",std::cout);

        helper.start("jacobian",std::cout);
        {
          auto engine = Dune::PDELab::JacobianEngine(solution,jacobian,local_operator);
          assembler.assemble(engine);
        }
        helper.end("jacobian",std::cout);

        helper.end_run(std::cout);
      }

    helper.print(std::cout);
  }

  {
    Vector diff(gfs,0.0);
    diff = r_old;
    diff -= r_wrapped;
    std::cout << "residual: old vs wrapped: " << diff.two_norm() << std::endl;
    diff = r_old;
    diff -= r_new;
    std::cout << "residual: old vs new   : " << diff.two_norm() << std::endl;
  }

  {
    using Dune::PDELab::Backend::native;
    using std::sqrt;
    Real norm(0.0);
    for (auto rit1 = native(j_old).begin(), rend = native(j_old).end(), rit2 = native(j_wrapped).begin() ;
         rit1 != rend ;
         ++rit1, ++rit2
         )
      for (auto cit1 = rit1->begin(), cend = rit1->end(), cit2 = rit2->begin() ;
           cit1 != cend ;
           ++cit1, ++cit2
           )
        {
          auto diff = *cit1;
          diff -= *cit2;
          norm += diff.frobenius_norm2();
        }
    std::cout << "jacobian: old vs wrapped: " << sqrt(norm) << std::endl;

    norm = 0.0;
    for (auto rit1 = native(j_old).begin(), rend = native(j_old).end(), rit2 = native(j_new).begin() ;
         rit1 != rend ;
         ++rit1, ++rit2
         )
      for (auto cit1 = rit1->begin(), cend = rit1->end(), cit2 = rit2->begin() ;
           cit1 != cend ;
           ++cit1, ++cit2
           )
        {
          auto diff = *cit1;
          diff -= *cit2;
          norm += diff.frobenius_norm2();
        }
    std::cout << "jacobian: old vs new: " << sqrt(norm) << std::endl;
  }

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

  Dune::ParameterTree params;
  Dune::ParameterTreeParser::readOptions(argc,argv,params);

  using Real = double;

  // Create 2D yasp grid
  constexpr int dim = 2;
  Dune::FieldVector<Real,dim> L(1.0);
  std::bitset<dim> P(false);
  using Grid = Dune::YaspGrid<dim>;

  auto cells_per_dir = params.get<int>("grid.cells_per_dir");
  Grid grid(L,{cells_per_dir,cells_per_dir},P,0);

  // Get GridView
  auto gv = grid.leafGridView();

  constexpr int degree = 1;
  using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,Real,degree,dim>;
  FEM fem;

  benchmark(params,gv,fem);

  return 0;
}
