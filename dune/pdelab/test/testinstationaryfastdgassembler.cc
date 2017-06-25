// -*- tab-width: 4; indent-tabs-mode: nil -*-

// NOTE: Right now the test compares with the exact solution afer one
// time step. It would be better to compare residuals with the non
// fast dg case.


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
#include <dune/pdelab/localoperator/l2.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/boilerplate/pdelab.hh>

#include<dune/pdelab/gridoperator/fastdg.hh>
#include"convectiondiffusionfastdg.hh"
#include"l2fastdg.hh"

#define FAST_DG

template<typename GV, typename RF>
class ParameterA
{
  const GV gv;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef RF RangeFieldType;
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  ParameterA( const GV gv_ ) : gv(gv_), time(0.0)
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

  //! set time for subsequent evaluation
  void setTime (typename Traits::RangeFieldType t)
  {
    time = t;
  }

private:
  RF time;
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
#ifdef FAST_DG
  typedef Dune::PDELab::ConvectionDiffusionFastDG<Problem,FEM> LOP;
#else
  typedef Dune::PDELab::ConvectionDiffusionDG<Problem,FEM> LOP;
#endif
  LOP lop(problem,m,w,2.0);

  // Local operator for mass matrix
#ifdef FAST_DG
  typedef Dune::PDELab::L2FastDG MLOP;
#else
  typedef Dune::PDELab::L2 MLOP;
#endif
  MLOP mlop(2*degree);

  // Constraints
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(9); // number of nonzeroes per row can be cross-checked by patternStatistics().
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;

  // GridOperator
#ifdef FAST_DG
  typedef Dune::PDELab::FastDGGridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO0;
#else
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO0;
#endif
  GO0 go0(gfs,cc,gfs,cc,lop,mbe);
#ifdef FAST_DG
  typedef Dune::PDELab::FastDGGridOperator<GFS,GFS,MLOP,MBE,Real,Real,Real,CC,CC> GO1;
#else
  typedef Dune::PDELab::GridOperator<GFS,GFS,MLOP,MBE,Real,Real,Real,CC,CC> GO1;
#endif
  GO1 go1(gfs,cc,gfs,cc,mlop,mbe);
  typedef Dune::PDELab::OneStepGridOperator<GO0, GO1> IGO;
  IGO igo(go0, go1);

  // Make a vector of degree of freedom vectors and initialize it with Dirichlet extension
  typedef typename IGO::Traits::Domain V;
  V x(gfs,0.0);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(gv,problem);
  Dune::PDELab::interpolate(g,gfs,x);

  // Make linear solver
  int ls_verbosity = 1;
  // CG only works for symmetric operator (SIPG). In case of (NIPG) use e.g.:
  //
  // typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0 LS;
  // LS ls(10000,ls_verbosity);
  typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
  LS ls(10000,ls_verbosity);
  typedef Dune::PDELab::StationaryLinearProblemSolver<IGO, LS,V> PDESOLVER;
  PDESOLVER pdesolver(igo,ls,1e-10);

  Dune::PDELab::Alexander2Parameter<Real> method;
  // Dune::PDELab::ImplicitEulerParameter<Real> method;
  Dune::PDELab::OneStepMethod<Real,IGO,PDESOLVER,V,V> osm(method,igo,pdesolver);
  osm.setVerbosityLevel(1);

  auto stationaryVTKWriter = std::make_shared<Dune::SubsamplingVTKWriter<GV> >(gv,degree-1);
  Dune::VTKSequenceWriter<GV> vtkwriter(stationaryVTKWriter,"testinstationaryfastdgassembler","","");
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(xdgf,"x_h"));

   Real time = 0.0;
   Real dt = 0.1;
   Real T = 0.1;
   vtkwriter.write(time,Dune::VTK::appendedraw);
   while (time<T-1e-10){
     // No time dependent parameters in problem
     // problem.setTime(time+dt);

     // do time step
     V xnew(gfs,0.0);
     osm.apply(time,dt,x,xnew);

     // accept time step
     x = xnew;
     time += dt;

     // VTK output
     vtkwriter.write(time,Dune::VTK::appendedraw);
   }

  // compute L2 error if analytical solution is available
  typedef Dune::PDELab::DifferenceSquaredAdapter<G,DGF> DifferenceSquared;
  DifferenceSquared differencesquared(g,xdgf);
  typename DifferenceSquared::Traits::RangeType l2errorsquared(0.0);
  Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,12);
  std::cout << "l2 error squared: " << l2errorsquared << std::endl;

  bool test_fail = false;
  if (l2errorsquared>5e-6)
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
    grid.globalRefine(3);

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
