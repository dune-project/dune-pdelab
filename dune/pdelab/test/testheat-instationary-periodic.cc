// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/localoperator/l2.hh>

#include <dune/pdelab/localoperator/convectiondiffusion.hh>

//***********************************************************************
//***********************************************************************
// diffusion problem with time dependent coefficients
//***********************************************************************
//***********************************************************************

const double kx = 3.0, ky = 3.0;
//! base class for parameter class
template<typename GV, typename RF>
class ConvectionDiffusionProblem :
  public Dune::PDELab::ConvectionDiffusionParameterInterface<
  Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF>,
  ConvectionDiffusionProblem<GV,RF>
  >
{
public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! source/reaction term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
     typename Traits::RangeFieldType u) const
  {
    return 0.0;
  }

  //! nonlinearity under gradient
  typename Traits::RangeFieldType
  w (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
     typename Traits::RangeFieldType u) const
  {
    return u;
  }

  //! nonlinear scaling of diffusion tensor
  typename Traits::RangeFieldType
  v (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
     typename Traits::RangeFieldType u) const
  {
    return 1.0;
  }

  //! tensor permeability
  typename Traits::PermTensorType
  D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType kabs;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        kabs[i][j] = (i==j) ? 1 : 0;
    return kabs;
  }

  //! nonlinear flux vector
  typename Traits::RangeType
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
     typename Traits::RangeFieldType u) const
  {
    typename Traits::RangeType flux;
    flux[0] = 0.0;
    flux[1] = 0.0;
    return flux;
  }

  //! boundary condition type function
  template<typename I>
  bool isDirichlet(
                   const I & intersection,               /*@\label{bcp:name}@*/
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
    return true;  // Dirichlet b.c. on all boundaries
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType x = e.geometry().global(xlocal);

    return std::exp(-(kx*kx+ky*ky)*M_PI*M_PI*time) * sin(kx*M_PI*x[0]) * sin(ky*M_PI*x[1]);
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
     typename Traits::RangeFieldType u) const
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
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> FS;
  FS fs(grid.leafGridView(),fem);

  // define problem parameters
  typedef ConvectionDiffusionProblem<GV,NumberType> Problem;
  Problem problem;
  Dune::PDELab::BCTypeParam_CD<Problem> bctype(grid.leafGridView(),problem);
  typedef Dune::PDELab::DirichletBoundaryCondition_CD<Problem> G;
  G g(grid.leafGridView(),problem);

  typedef typename FS::template ConstraintsContainer<NumberType>::Type C;
  C cg;
  Dune::PDELab::constraints(bctype, fs, cg);

  // make grid operator space for time-dependent problem
  typedef Dune::PDELab::ConvectionDiffusion<Problem> LOP;
  LOP lop(problem,4);
  typedef Dune::PDELab::L2 MLOP;
  MLOP mlop(4);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
  //Dune::PDELab::FractionalStepParameter<Real> method;
  Dune::PDELab::Alexander3Parameter<NumberType> method;
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

  typedef Dune::PDELab::Newton<IGO,LS,V> PDESOLVER;
  PDESOLVER pdesolver(igo,ls);
  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(0);
  pdesolver.setReduction(0.9);
  pdesolver.setMinLinearReduction(1e-9);

  Dune::PDELab::OneStepMethod<NumberType,IGO,PDESOLVER,V,V> osm(method,igo,pdesolver);
  osm.setVerbosityLevel(2);

  // graphics for initial guess
  Dune::PDELab::FilenameHelper fn(basename);
  { // start a new block to automatically delete the VTKWriter object
    typedef Dune::PDELab::DiscreteGridFunction<FS,V> DGF;
    Dune::SubsamplingVTKWriter<typename GM::LeafGridView> vtkwriter(grid.leafGridView(),degree-1);
    DGF xdgf(fs,x);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"x_h"));
    vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
    fn.increment();
  }

  // time loop
  NumberType time = 0.0;
  while (time<T-1e-10)
    {
      problem.setTime(time+dt);

      // do time step
      V xnew(fs,0.0);
      osm.apply(time,dt,x,xnew);

      // output to VTK file
      {
        typedef Dune::PDELab::DiscreteGridFunction<FS,V> DGF;
        Dune::SubsamplingVTKWriter<typename GM::LeafGridView> vtkwriter(grid.leafGridView(),degree-1);
        DGF xdgf(fs,xnew);
        vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"x_h"));
        vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
        fn.increment();
      }

      // accept time step
      x = xnew;
      time += dt;
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
  Dune::MPIHelper& mpihelper = Dune::MPIHelper::instance(argc,argv);

  double T = 0.075;
  double dt = 0.0075;
  int cells = 32;

  // start try/catch block to get error messages from dune
  try {

    const int dim=2;
    const int degree=1;
    const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::overlapping;
    const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
    const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::conforming;

    typedef Dune::YaspGrid<dim> GM;
    Dune::FieldVector<double,dim> L(1.0);
    Dune::array<int,dim> N(Dune::fill_array<int,dim>(8));

    std::bitset<dim> periodic (false);
    periodic[0] = true;
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
    std::cout << "STL ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (Dune::Exception & e) {
    std::cout << "DUNE ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cout << "Unknown ERROR" << std::endl;
    return 1;
  }

  // done
  return 0;
}
