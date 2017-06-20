// -*- tab-width: 4; indent-tabs-mode: nil -*-

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/solvers.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/pdelab/finiteelementmap/qkdg.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include<dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>

#include<dune/pdelab/gridoperator/fastdg.hh>
#include"convectiondiffusionfastdg.hh"

template<typename GV, typename RF>
class ParameterA
{
  const GV gv;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef RF RangeFieldType;
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  ParameterA( const GV gv_ ) : gv(gv_)
  {
  }

  std::string name() const {return "A";};

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


//! solve problem with DG method
template<class GV, class FEM, class Problem, int degree>
bool runDG(const GV& gv, const FEM& fem, Problem& problem)
{
  // Coordinate and result type
  typedef typename Problem::RangeFieldType Real;
  const int dim = GV::Grid::dimension;

  // Make grid function space
  typedef Dune::PDELab::NoConstraints CON;
  const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
  typedef Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed,blocksize> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  // Make local operator
  Dune::PDELab::ConvectionDiffusionDGMethod::Type m = Dune::PDELab::ConvectionDiffusionDGMethod::SIPG;
  Dune::PDELab::ConvectionDiffusionDGWeights::Type w = Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn;
  // typedef Dune::PDELab::ConvectionDiffusionDG<Problem,FEM> LOP;
  typedef Dune::PDELab::ConvectionDiffusionFastDG<Problem,FEM> LOP;
  LOP lop(problem,m,w,2.0);

  // Constraints
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(9); // number of nonzeroes per row can be cross-checked by patternStatistics().
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;

  // GridOperator
  // typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  typedef Dune::PDELab::FastDGGridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // Make a vector of degree of freedom vectors and initialize it with Dirichlet extension
  typedef typename GO::Traits::Domain U;
  U u(gfs,0.0);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(gv,problem);

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
  typedef Dune::PDELab::DifferenceSquaredAdapter<G,UDGF> DifferenceSquared;
  DifferenceSquared differencesquared(g,udgf);
  typename DifferenceSquared::Traits::RangeType l2errorsquared(0.0);
  Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,12);
  std::cout << "l2 error squared: " << l2errorsquared << std::endl;

  // write vtk file
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,degree-1);
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<UDGF>>(udgf,"u_h"));
  vtkwriter.write("testfastdgassembler",Dune::VTK::appendedraw);

  bool test_fail = false;
  if (l2errorsquared>1e-08)
    test_fail = true;
  return test_fail;
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

  try{
    typedef double Real;

    // Create 2D yasp grid
    const int dim = 2;
    Dune::FieldVector<Real,dim> L(1.0);
    Dune::array<int,dim> N(Dune::fill_array<int,dim>(1));
    std::bitset<dim> P(false);
    typedef Dune::YaspGrid<dim> Grid;
    Grid grid(L,N,P,0);

    // Refine grid
    grid.globalRefine(5);

    // Get GridView
    typedef Grid::LeafGridView GV;
    const GV gv = grid.leafGridView();

    // Create problem
    typedef ParameterA<GV,Real> Problem;
    Problem problem(gv);

    // Create DG space
    const int degree=1;
    typedef Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,Real,degree,dim> FEMDG;
    FEMDG femdg;

    // Solve problem
    test_fail = runDG <GV, FEMDG, Problem, degree>(gv, femdg, problem);
  }

  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }

  return test_fail;
}
