// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/localoperator/l2.hh>
#include "dune/pdelab/gridfunctionspace/gridfunctionadapter.hh"

// Poisson problem
template<typename GV, typename RF>
class Poisson
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  Poisson () : time(0.0) {}

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1.0 : 0.0;
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
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType x = e.geometry().global(xlocal);
    Dune::FieldVector<double,2> c(0.5);
    c-= x;
    return 4.*(1.-c.two_norm2())*std::exp(-1.*c.two_norm2());
    // return 0.0;
  }

  //! boundary condition type function
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType x = e.geometry().global(xlocal);
    Dune::FieldVector<double,2> c(0.5);
    c-= x;
    return std::exp(-1.*c.two_norm2());
    // return std::exp(-(kx*kx+ky*ky)*M_PI*M_PI*time) * sin(kx*M_PI*x[0]) * sin(ky*M_PI*x[1]);
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
    //std::cout << "setting time to " << time << std::endl;
  }

private:
  RF time;
};


// Solve poisson problem
template<typename GM, unsigned int degree, Dune::GeometryType::BasicType elemtype,
         Dune::PDELab::MeshType meshtype, Dune::SolverCategory::Category solvertype>
bool do_simulation (double T, double dt, GM& grid, std::string basename)
{
  // define parameters
  typedef double NumberType;

  // make problem parameters
  typedef Poisson<typename GM::LeafGridView,NumberType> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(grid.leafGridView(),problem);

  // make a finite element space
  typedef Dune::PDELab::CGSpace<GM,NumberType,degree,BCType,elemtype,meshtype,solvertype> FS;
  FS fs(grid,bctype);

  // assemblers for finite element problem
  typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,typename FS::FEM> LOP;
  LOP lop(problem,4);
  typedef Dune::PDELab::GalerkinGlobalAssembler<FS,LOP,solvertype> SASS;
  SASS sass(fs,lop,6);
  typedef Dune::PDELab::L2 MLOP;
  MLOP mlop(2*degree);
  typedef Dune::PDELab::GalerkinGlobalAssembler<FS,MLOP,solvertype> TASS;
  TASS tass(fs,mlop,6);
  typedef Dune::PDELab::OneStepGlobalAssembler<SASS,TASS> ASSEMBLER_IMPLICIT;
  ASSEMBLER_IMPLICIT assembler_implicit(sass,tass);
  typedef Dune::PDELab::OneStepGlobalAssembler<SASS,TASS,false> ASSEMBLER_EXPLICIT;
  ASSEMBLER_EXPLICIT assembler_explicit(sass,tass);

  // make a degree of freedom vector and set initial value
  typedef typename FS::DOF V;
  V x_implicit(fs.getGFS(),0.0);
  V x_explicit(fs.getGFS(),0.0);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(grid.leafGridView(),problem);
  problem.setTime(0.0);
  Dune::PDELab::interpolate(g,fs.getGFS(),x_implicit);
  Dune::PDELab::interpolate(g,fs.getGFS(),x_explicit);

  // linear solver backend
  typedef Dune::PDELab::ISTLSolverBackend_CG_AMG_SSOR<FS,ASSEMBLER_IMPLICIT,solvertype> SBE_IMPLICIT;
  SBE_IMPLICIT sbe_implicit(fs,assembler_implicit,5000,1);
  typedef Dune::PDELab::ISTLSolverBackend_CG_AMG_SSOR<FS,ASSEMBLER_EXPLICIT,solvertype> SBE_EXPLICIT;
  SBE_EXPLICIT sbe_explicit(fs,assembler_explicit,5000,1);

  // linear problem solver
  typedef Dune::PDELab::StationaryLinearProblemSolver<typename ASSEMBLER_IMPLICIT::GO,typename SBE_IMPLICIT::LS,V>
    PDESOLVER_IMPLICIT;
  PDESOLVER_IMPLICIT pdesolver_implicit(*assembler_implicit,*sbe_implicit,1e-6);
  typedef Dune::PDELab::StationaryLinearProblemSolver<typename ASSEMBLER_EXPLICIT::GO,typename SBE_EXPLICIT::LS,V>
    PDESOLVER_EXPLICIT;
  PDESOLVER_EXPLICIT pdesolver_explicit(*assembler_explicit,*sbe_explicit,1e-6);

  // time-stepper
  Dune::PDELab::OneStepThetaParameter<NumberType> method_implicit(1.0);
  typedef Dune::PDELab::OneStepMethod<NumberType,typename ASSEMBLER_IMPLICIT::GO,PDESOLVER_IMPLICIT,V>
    OSM_IMPLICIT;
  OSM_IMPLICIT osm_implicit(method_implicit,*assembler_implicit,pdesolver_implicit);
  osm_implicit.setVerbosityLevel(1);
  Dune::PDELab::ExplicitEulerParameter<NumberType> method_explicit;
  typedef Dune::PDELab::ExplicitOneStepMethod<NumberType,typename ASSEMBLER_EXPLICIT::GO,typename SBE_EXPLICIT::LS,V,V>
    OSM_EXPLICIT;
  OSM_EXPLICIT osm_explicit(method_explicit,*assembler_explicit,*sbe_explicit);
  osm_explicit.setVerbosityLevel(1);

  // graphics for initial guess
  auto stationaryVTKWriter = std::make_shared<Dune::SubsamplingVTKWriter<typename GM::LeafGridView> >(grid.leafGridView(),degree-1);
  Dune::VTKSequenceWriter<typename GM::LeafGridView> vtkwriter(stationaryVTKWriter,basename,"","");
  typedef Dune::PDELab::DiscreteGridFunction<typename FS::GFS,V> DGF;
  DGF xdgf_implicit(fs.getGFS(),x_implicit);
  DGF xdgf_explicit(fs.getGFS(),x_explicit);
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(xdgf_implicit,"u_h_implicit"));
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(xdgf_explicit,"u_h_explicit"));

  // time loop
  NumberType time = 0.0;
  while (time<T-1e-10){
    // assemble constraints for new time step (assumed to be constant for all substeps)
    problem.setTime(time+dt);
    fs.assembleConstraints(bctype);

    // do time step
    V xnew_implicit(fs.getGFS(),0.0);
    V xnew_explicit(fs.getGFS(),0.0);
    osm_implicit.apply(time,dt,x_implicit,g,xnew_implicit);
    osm_explicit.apply(time,dt,x_explicit,g,xnew_explicit);

    // output to VTK file
    vtkwriter.write(time,Dune::VTK::appendedraw);

    // accept time step
    x_implicit = xnew_implicit;
    x_explicit = xnew_explicit;
    time += dt;
  }

  // Compare implicit and explicit solutions with interpolation of exact solution
  DGF xdgf_implicit_after(fs.getGFS(),x_implicit);
  DGF xdgf_explicit_after(fs.getGFS(),x_explicit);
  typedef Dune::PDELab::DifferenceSquaredAdapter<G,DGF> DifferenceSquared;
  DifferenceSquared differencesquared_implicit(g, xdgf_implicit_after);
  DifferenceSquared differencesquared_explicit(g, xdgf_explicit_after);
  typename DifferenceSquared::Traits::RangeType l2error_implicit(0.0);
  typename DifferenceSquared::Traits::RangeType l2error_explicit(0.0);
  Dune::PDELab::integrateGridFunction(differencesquared_implicit, l2error_implicit, 10);
  Dune::PDELab::integrateGridFunction(differencesquared_explicit, l2error_explicit, 10);
  std::cout << "Accumlated L2 error for implicit time stepping: " << l2error_implicit << std::endl;
  std::cout << "Accumlated L2 error for explicit time stepping: " << l2error_explicit << std::endl;

  // Compare implicit and explicit version
  typedef Dune::PDELab::DifferenceSquaredAdapter<DGF,DGF> DifferenceSquaredCompare;
  DifferenceSquaredCompare differencesquared_compare(xdgf_implicit_after, xdgf_explicit_after);
  typename DifferenceSquaredCompare::Traits::RangeType l2error_compare(0.0);
  Dune::PDELab::integrateGridFunction(differencesquared_compare, l2error_compare, 10);
  std::cout << "Accumlated L2 norm of difference between implicit and explicit: " << l2error_compare << std::endl;

  // Decide wether the errors are too big
  bool fail = false;
  if (l2error_implicit>1e-3 || l2error_explicit>1e-3 || l2error_compare>1e-10)
     fail = true;

  return fail;
}


int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  double T = 0.01;
  double dt = 0.001;
  int cells = 4;

  // start try/catch block to get error messages from dune
  try {

    const int dim=2;
    const int degree=1;
    const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::sequential;
    const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
    const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::conforming;

    typedef Dune::YaspGrid<dim> GM;
    typedef Dune::PDELab::StructuredGrid<GM> Grid;
    Grid grid(elemtype,cells);
    grid->loadBalance();

    std::stringstream basename;
    basename << "test_instationary_with_boundary_constraints";
    bool fail = do_simulation<GM,degree,elemtype,meshtype,solvertype>(T,dt,*grid,basename.str());
    if (fail)
      return 1;
  }
  catch (Dune::Exception & e) {
    std::cout << "DUNE ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (std::exception & e) {
    std::cout << "STL ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cout << "Unknown ERROR" << std::endl;
    return 1;
  }

  // done
  return 0;
}
