// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/localoperator/l2.hh>

#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>

//***********************************************************************
//***********************************************************************
// diffusion problem with time dependent coefficients
//***********************************************************************
//***********************************************************************

/** \brief Parameter class for solving the linear convection-diffusion equation
 *
 * A parameter class for the linear convection-diffusion equation
 * \f{align*}{
 *   -\nabla\cdot(A(x) \nabla u) + b(x)\cdot \nabla u + c(x)u &=& f \mbox{ in } \Omega,  \\
 *                                                          u &=& g \mbox{ on } \partial\Omega_D \\
 *                            (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
 *                                    -(A(x)\nabla u) \cdot n &=& o \mbox{ on } \partial\Omega_O
 * \f}
 * Note:
 *  - This formulation is valid for velocity fields which are non-divergence free.
 *  - Outflow boundary conditions should only be set on the outflow boundary
 *
 * \tparam T a traits class defining the necessary types
 */
const double kx = 3.0, ky = 3.0;
template<typename GV, typename RF>
class ConvectionDiffusionModelProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

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
    return 0.0;
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
    return std::exp(-(kx*kx+ky*ky)*M_PI*M_PI*time) * sin(kx*M_PI*xglobal[0]) * sin(ky*M_PI*xglobal[1]);
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
  void setTime (RF t)
  {
    time = t;
  }

private:
  RF time;
};

//***********************************************************************
//***********************************************************************
// a function that does the simulation on a given grid
//***********************************************************************
//***********************************************************************

template<typename GM, unsigned int degree, Dune::GeometryType::BasicType elemtype,
         Dune::PDELab::MeshType meshtype, Dune::SolverCategory::Category solvertype>
void do_simulation (double T, double dt, GM& grid, std::string basename)
{
  // define parameters
  typedef double NumberType;

  typedef typename GM::LeafGridView GV;
  typedef typename GV::Grid::ctype Coord;
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,NumberType,degree> FEM;
  FEM fem(grid.leafGridView());
  typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON; // ovlp
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> FS;
  FS fs(grid.leafGridView(),fem);

  // define problem parameters
  typedef ConvectionDiffusionModelProblem<GV,NumberType> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(grid.leafGridView(),problem);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(grid.leafGridView(),problem);

  typedef typename FS::template ConstraintsContainer<NumberType>::Type C;
  C cg;
  Dune::PDELab::constraints(bctype, fs, cg);

  // make grid operator space for time-dependent problem
  typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
  LOP lop(problem,1);
  typedef Dune::PDELab::L2 MLOP;
  MLOP mlop(2*degree + 2);
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
  //Dune::PDELab::FractionalStepParameter<Real> method;
  Dune::PDELab::Alexander2Parameter<NumberType> method;
  typedef Dune::PDELab::GridOperator<FS,FS,LOP,MBE,NumberType,NumberType,NumberType,C,C> GO0;
  GO0 go0(fs,cg,fs,cg,lop,mbe);
  typedef Dune::PDELab::GridOperator<FS,FS,MLOP,MBE,NumberType,NumberType,NumberType,C,C> GO1;
  GO1 go1(fs,cg,fs,cg,mlop,mbe);
  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);
  typedef typename IGO::Traits::Domain V;

  // make a degree of freedom vector and set initial value
  V x(fs,0.0);
  problem.setTime(0.0);
  Dune::PDELab::interpolate(g,fs,x);

  // linear problem solver
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILU0<FS, C> LS;
  LS ls(fs, cg);

  typedef Dune::PDELab::StationaryLinearProblemSolver<IGO,LS,V> PDESOLVER;
  PDESOLVER pdesolver(igo,ls,1e-10);

  Dune::PDELab::OneStepMethod<NumberType,IGO,PDESOLVER,V,V> osm(method,igo,pdesolver);
  osm.setVerbosityLevel(2);

  auto stationaryVTKWriter = std::make_shared<Dune::SubsamplingVTKWriter<typename GM::LeafGridView> >(grid.leafGridView(),degree-1);
  Dune::VTKSequenceWriter<typename GM::LeafGridView> vtkwriter(stationaryVTKWriter,basename,"","");
  typedef Dune::PDELab::DiscreteGridFunction<FS,V> DGF;
  DGF xdgf(fs,x);
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(xdgf,"x_h"));

  // time loop
  NumberType time = 0.0;
  vtkwriter.write(time,Dune::VTK::appendedraw);
  while (time<T-1e-10)
    {
      problem.setTime(time+dt);

      // do time step
      V xnew(fs,0.0);
      osm.apply(time,dt,x,xnew);

      // accept time step
      x = xnew;
      time += dt;

      // VTK output
      vtkwriter.write(time,Dune::VTK::appendedraw);
    }
}

//***********************************************************************
//***********************************************************************
// the main function
//***********************************************************************
//***********************************************************************

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  double T = 0.075;
  double dt = 0.0075;
  int cells = 8;

  // start try/catch block to get error messages from dune
  try {

    const int dim=2;
    const int degree=1;
    const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::overlapping;
    const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
    const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::conforming;

    typedef Dune::YaspGrid<dim> GM;
    Dune::FieldVector<double,dim> L(1.0);
    std::array<int,dim> N(Dune::fill_array<int,dim>(cells));

    std::bitset<dim> periodic (false);
    periodic[0] = true;
    periodic[1] = true;
    int overlap = 1;

    GM grid(L,N,periodic, overlap);
    grid.refineOptions (true);
    grid.globalRefine (1);
    grid.loadBalance();

    std::stringstream basename;
    basename << "heat_instationary_periodic" << "_dim" << dim << "_degree" << degree;
    do_simulation<GM,degree,elemtype,meshtype,solvertype>(T,dt,grid,basename.str());
  }
  catch (std::exception & e) {
    std::cout << "Exception: " << e.what() << std::endl;
    throw;
  }
  catch (...) {
    std::cout << "Unknown ERROR" << std::endl;
    return 1;
  }

  // done
  return 0;
}
