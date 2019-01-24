// -*- tab-width: 4; indent-tabs-mode: nil -*-

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
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

#include <dune/pdelab/assembler/gridoperator.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondgnew.hh>
#include <dune/pdelab/localoperator/adapter.hh>


template<typename GV, typename RF>
class ParameterA
{
  const GV gv;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef RF RangeFieldType;
  using Real = RF;
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  ParameterA( const GV gv_ ) : gv(gv_)
  {
  }

  std::string name() const {return "A";};

  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static bool permeabilityIsConstantPerCell()
  {
    return false;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1 : 0;
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    typename Traits::RangeFieldType norm = xglobal.two_norm2();
    return (2.0*GV::dimension-4.0*norm)*exp(-norm);
  }

  //! boundary condition type function
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    typename Traits::RangeFieldType norm = xglobal.two_norm2();
    return exp(-norm);
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
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
        for (std::size_t i=0; i < LOP::dimension<Context>; i++)
          for (std::size_t j=0; j < LOP::dimension<Context>; j++)
            I[i][j] = (i==j) ? 1 : 0;
        return I;
      }

      //! velocity field
      template<typename P>
      typename Traits::Velocity b(const P& p)
      {
        return typename Traits::Velocity(0.0);
      }

      //! sink term
      template<typename P>
      Real c(const P& p)
      {
        return 0.0;
      }

      //! source term
      template<typename P>
      LOP::RangeField<Context> f(const P& p)
      {
        auto x = Context::global(p);
        auto norm = x.two_norm2();
        using std::exp;
        return (2.0 * LOP::dimension<Context> - 4.0 * norm) * exp(-norm);
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
        using std::exp;
        return exp(-Context::global(p).two_norm2());
      }

      //! Neumann boundary condition
      template<typename P>
      LOP::RangeField<Context> j(const P& p)
      {
        return 0.0;
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



//! solve problem with DG method
template<class GV, class FEM, class Problem, class OldProblem, int degree>
bool runDG(const GV& gv, const FEM& fem, Problem& problem, OldProblem& old_problem)
{
  // Coordinate and result type
  typedef typename Problem::Real Real;
  const int dim = GV::Grid::dimension;

  // Make grid function space
  typedef Dune::PDELab::NoConstraints CON;
  const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
  typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,blocksize> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("u_h");

  // Make local operator
  Dune::PDELab::ConvectionDiffusionDGMethod::Type m = Dune::PDELab::ConvectionDiffusionDGMethod::SIPG;
  Dune::PDELab::ConvectionDiffusionDGWeights::Type w = Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn;

  typedef Dune::PDELab::ConvectionDiffusionDG<OldProblem,FEM> OldLOP;
  auto old_lop = OldLOP{old_problem, m, w, 2.0};
  using LOPWrapper = Dune::PDELab::LocalOperatorAdapter<OldLOP>;
  auto lop_wrapper = LOPWrapper{old_lop};

  using LOP = Dune::PDELab::ConvectionDiffusionDGNew<Problem>;
  auto lop = LOP{problem, m, w, 2.0};


  // Constraints
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(9); // number of nonzeroes per row can be cross-checked by patternStatistics().
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;

  // GridOperator
  // typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  typedef Dune::PDELab::NewGridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  using OldGO = Dune::PDELab::NewGridOperator<GFS,GFS,LOPWrapper,MBE,Real,Real,Real,CC,CC>;
  OldGO old_go(gfs,cc,gfs,cc,lop_wrapper,mbe);


  // Make a vector of degree of freedom vectors and initialize it with Dirichlet extension
  typedef typename GO::Traits::Domain U;
  U u(gfs,0.0);
  //typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  //G g(gv,problem);

  {
    U r(gfs,0.0);
    U r_old(gfs,0.0);

    go.residual(u,r);
    old_go.residual(u,r_old);

    for (auto it = r.begin(), it2 = r_old.begin(), end = r.end();
         it != end;
         ++it, ++it2
      )
      if (*it != *it2)
        return true;
  }


  {
    typename GO::Traits::Jacobian j(go,0.0);
    typename GO::Traits::Jacobian j_old(old_go,0.0);

    go.jacobian(u,j);
    old_go.jacobian(u,j_old);
  }

  {
    // Make linear solver
    int ls_verbosity = 2;
    // CG only works for symmetric operator (SIPG). In case of (NIPG) use e.g.:
    //
    // typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0 LS;
    // LS ls(10000,ls_verbosity);
    typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
    LS ls(10000,ls_verbosity);

    // Solve problem
    typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
    SLP slp(go,ls,u,1e-12);
    slp.apply();

    // compute L2 error if analytical solution is available
    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> UDGF;
    UDGF udgf(gfs,u);
    /*
      typedef Dune::PDELab::DifferenceSquaredAdapter<G,UDGF> DifferenceSquared;
      DifferenceSquared differencesquared(g,udgf);
      typename DifferenceSquared::Traits::RangeType l2errorsquared(0.0);
      Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,12);
      std::cout << "l2 error squared: " << l2errorsquared << std::endl;
    */
    // write vtk file
    Dune::SubsamplingVTKWriter<typename GV::GridView> vtkwriter(gv.gridView(),Dune::refinementIntervals(degree));
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,u);
    //vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<UDGF>>(udgf,"u_h"));
    vtkwriter.write("testnewassembler",Dune::VTK::appendedraw);
  }

  {
    // Make linear solver
    int ls_verbosity = 2;
    // CG only works for symmetric operator (SIPG). In case of (NIPG) use e.g.:
    //
    // typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0 LS;
    // LS ls(10000,ls_verbosity);
    typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
    LS ls(10000,ls_verbosity);

    // Solve problem
    typedef Dune::PDELab::StationaryLinearProblemSolver<OldGO,LS,U> SLP;
    SLP slp(old_go,ls,u,1e-12);
    slp.apply();

    // compute L2 error if analytical solution is available
    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> UDGF;
    UDGF udgf(gfs,u);
    /*
      typedef Dune::PDELab::DifferenceSquaredAdapter<G,UDGF> DifferenceSquared;
      DifferenceSquared differencesquared(g,udgf);
      typename DifferenceSquared::Traits::RangeType l2errorsquared(0.0);
      Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,12);
      std::cout << "l2 error squared: " << l2errorsquared << std::endl;
    */
    // write vtk file
    Dune::SubsamplingVTKWriter<typename GV::GridView> vtkwriter(gv.gridView(),Dune::refinementIntervals(degree));
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,u);
    //vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<UDGF>>(udgf,"u_h"));
    vtkwriter.write("testnewassembler_old",Dune::VTK::appendedraw);
  }

  /*
  bool test_fail = false;
  if (l2errorsquared>1e-08)
    test_fail = true;
  return test_fail;
  */
  return false;
}


int main(int argc, char** argv)
{
  //Maybe initialize Mpi
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  if(Dune::MPIHelper::isFake)
    std::cout<< "This is a sequential program." << std::endl;
  else{
    if(helper.rank()==0)
      std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
  }

  bool test_fail = false;

  typedef double Real;

  // Create 2D yasp grid
  const int dim = 2;
  Dune::FieldVector<Real,dim> L(1.0);
  std::bitset<dim> P(false);
  typedef Dune::YaspGrid<dim> Grid;
  Grid grid(L,Dune::filledArray<dim, int>(1),P,0);

  // Refine grid
  grid.globalRefine(1);

  // Get GridView
  typedef Grid::LeafGridView GV;
  const GV gv = grid.leafGridView();

  using ES = Dune::PDELab::AllEntitySet<GV>;
  ES es(gv);

  // Create problem

  using OldProblem = ParameterA<ES,Real>;
  OldProblem old_problem(es);
  using Problem = DG::Problem<Real>;
  Problem problem;


  // Create DG space
  const int degree=1;
  typedef Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,Real,degree,dim> FEMDG;
  FEMDG femdg;

  // Solve problem
  return runDG <ES, FEMDG, Problem, OldProblem, degree>(es, femdg, problem, old_problem);
}
