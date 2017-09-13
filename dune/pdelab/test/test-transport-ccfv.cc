// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Test case for "convectiondiffusionccfv.hh". Solve linear transport problem with cell-centered finite volumes in stationary and instationary (explicit and implicit). Newton is being tested here only for the linear case.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/typetraits.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/backend/istl.hh>

#include<dune/pdelab/localoperator/convectiondiffusionccfv.hh>

#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>


//==============================================================================
// Parameter class for the convection diffusion problem
//==============================================================================

//! base class for parameter class

template<typename GV, typename RF>
class TransportProblem //:
// public Dune::PDELab::TransportSpatialParameterInterface<Dune::PDELab::TransportParameterTraits<GV,RF>,
//                                                          TransportProblem<GV,RF> >
{

private:
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
  RF time;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;
  //typedef Dune::PDELab::TransportParameterTraits<GV,RF> Traits;

  //! tensor permeability
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    RF k = 1E-10;

    for (int i=0; i<GV::dimension; i++)
      for (int j=0; j<GV::dimension; j++)
        if (i==j)
          I[i][i] = k;
        else
          I[i][j] = 0.0;

    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType flux;
    flux[0] = 1.0;
    flux[1] = 0.5;
    return flux;
  }

  //! reaction term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& ) const
  {
    return 1.0;
  }

  //! capacity term
  typename Traits::RangeFieldType
  d (const typename Traits::ElementType& e, const typename Traits::DomainType& ) const
  {
    return 1.0;
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
    typename Traits::DomainType xglobal = is.geometry().global(x);
    if (xglobal[0]<1E-6 || xglobal[1]<1E-6)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow;
  }


  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType global = e.geometry().global(x);
    typename Traits::RangeFieldType k = std::abs( std::cos( 2.0 * time ) );
    if (global[0]>1E-12 && global[0]<1.0-1E-12 &&           // this is needed for parallel
        global[1]>1E-12 && global[1]<1.0-1E-12) return 0.0; // overlapping version !
    if (global[0]<1E-6 && global[1]>0.25 && global[1]<0.5+k*0.25)
      return k; // 1.0;
    else
      return 0.0;
  }

  //! Neumann boundary condition
  // Good: The dependence on u allows us to implement Robin type boundary conditions.
  // Bad: This interface cannot be used for mixed finite elements where the flux is the essential b.c.
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary function
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
};


//===============================================================
// Some variants to solve the linear convection diffusion problem
//===============================================================

// example solving stationary transport problem
template<class GV>
void stationary (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  auto gt = Dune::GeometryTypes::cube(dim);
  typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem(gt); // works only for cubes
  typedef Dune::PDELab::P0ParallelConstraints CON;
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> define problem parameters
  typedef TransportProblem<GV,Real> Param;

  Param param;
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Param> bctype(gv,param);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Param> G;
  G g(gv,param);

  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(bctype,gfs,cc);
  std::cout << "constrained dofs=" << cc.size()
            << " of " << gfs.globalSize() << std::endl;

  // <<<4>>> Make grid operator
  typedef Dune::PDELab::ConvectionDiffusionCCFV<Param> LOP;
  LOP lop(param);
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // <<<5>>> Compute affine shift
  typedef typename GO::Traits::Domain V;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cc,0.0,x);

  // <<<6>>> Make a linear solver
#if HAVE_SUPERLU
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<GO> LS;
  LS ls;
#else
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);
#endif

  // <<<7>>> solve nonlinear problem
  Dune::PDELab::Newton<GO,LS,V> newton(go,x,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.setReduction(0.9);
  newton.setMinLinearReduction(1e-9);
  newton.apply();

  // <<<8>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(xdgf,"solution"));
  vtkwriter.write("test_transport_ccfv_stationary",Dune::VTK::appendedraw);

}


// example using implicit time-stepping
template<class GV>
void implicit_scheme (const GV& gv, double Tend, double timestep)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  auto gt = Dune::GeometryTypes::cube(dim);
  typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem(gt); // works only for cubes
  typedef Dune::PDELab::P0ParallelConstraints CON;
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> define problem parameters
  typedef TransportProblem<GV,Real> Param;
  Param param;
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Param> bctype(gv,param);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Param> G;
  G g(gv,param);

  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;

  Dune::PDELab::constraints(bctype,gfs,cc);
  std::cout << "constrained dofs=" << cc.size()
            << " of " << gfs.globalSize() << std::endl;

  // <<<4>>> Make grid operator
  typedef Dune::PDELab::ConvectionDiffusionCCFV<Param> LOP;
  LOP lop(param);
  typedef Dune::PDELab::ConvectionDiffusionCCFVTemporalOperator<Param> SLOP;
  SLOP slop(param);
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // number of nonzero entries per row can be cross-checked by patternStatistics().
  Dune::PDELab::FractionalStepParameter<Real> method;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO0;
  GO0 go0(gfs,cc,gfs,cc,lop,mbe);
  typedef Dune::PDELab::GridOperator<GFS,GFS,SLOP,MBE,Real,Real,Real,CC,CC> GO1;
  GO1 go1(gfs,cc,gfs,cc,slop,mbe);
  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);

  // <<<5>>> Compute affine shift
  typedef typename IGO::Traits::Domain V;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cc,0.0,x);

  // <<<6>>> Make a linear solver
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS;
  LS ls(gfs,cc,5000,1,1);

  // <<<7>>> make Newton for time-dependent problem
  typedef Dune::PDELab::Newton<IGO,LS,V> PDESOLVER;
  PDESOLVER tnewton(igo,ls);
  tnewton.setReassembleThreshold(0.0);
  tnewton.setVerbosityLevel(0);
  tnewton.setReduction(0.9);
  tnewton.setMinLinearReduction(1e-7);

  // <<<8>>> time-stepper
  Dune::PDELab::OneStepMethod<Real,IGO,PDESOLVER,V,V> osm(method,igo,tnewton);
  osm.setVerbosityLevel(2);

  // <<<9>>> initial value and initial value for first time step with b.c. set
  V xold(gfs,0.0);
  xold = 0.0;

  // <<<10>>> graphics for initial guess
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,xold);
  Dune::VTKSequenceWriter<GV> vtkwriter(gv,"test_transport_ccfv_implicit_fs","vtk","");
  vtkwriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(xdgf,"solution"));
  vtkwriter.write(0.0,Dune::VTK::appendedraw);


  // <<<11>>> time loop
  Real time = 0.0;
  Real dt = timestep;
  x = 0;
  while (time < Tend)
    {
      // do time step
      osm.apply(time,dt,xold,x);

      // graphics
      xold = x;
      time += dt;
      vtkwriter.write(time,Dune::VTK::appendedraw);

    }
}


// example using explicit time-stepping
template<class GV>
void explicit_scheme (const GV& gv, double Tend, double timestep)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  auto gt = Dune::GeometryTypes::cube(dim);
  typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem(gt); // works only for cubes
  typedef Dune::PDELab::P0ParallelConstraints CON;
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> define problem parameters
  typedef TransportProblem<GV,Real> Param;
  Param param;
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Param> bctype(gv,param);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Param> G;
  G g(gv,param);

  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(bctype,gfs,cc);
  std::cout << "constrained dofs=" << cc.size()
            << " of " << gfs.globalSize() << std::endl;


  // <<<5>>> Make grid operator space
  typedef Dune::PDELab::ConvectionDiffusionCCFV<Param> LOP;
  LOP lop(param);
  typedef Dune::PDELab::ConvectionDiffusionCCFVTemporalOperator<Param> SLOP;
  SLOP slop(param);
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // number of nonzero entries per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO0;
  GO0 go0(gfs,cc,gfs,cc,lop,mbe);
  Dune::PDELab::ExplicitEulerParameter<Real> method;
  typedef Dune::PDELab::GridOperator<GFS,GFS,SLOP,MBE,Real,Real,Real,CC,CC> GO1;
  GO1 go1(gfs,cc,gfs,cc,slop,mbe);
  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1,false> IGO;
  IGO igo(go0,go1);
  igo.setMethod(method);

  // <<<4>>> Compute affine shift
  typedef typename IGO::Traits::Domain V;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cc,0.0,x);

  // <<<6>>> Make a linear solver backend
  typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS> LS;
  LS ls(gfs);

  // <<<8>>> time-stepper
  typedef Dune::PDELab::CFLTimeController<Real,IGO> TC;
  TC tc(0.999,igo);
  Dune::PDELab::ExplicitOneStepMethod<Real,IGO,LS,V,V,TC> osm(method,igo,ls,tc);
  osm.setVerbosityLevel(2);

  // <<<9>>> initial value and initial value for first time step with b.c. set
  V xold(gfs,0.0);
  xold = 0.0;

  // <<<10>>> graphics for initial guess
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,xold);
  Dune::VTKSequenceWriter<GV> vtkwriter(gv,"test_transport_ccfv_explicit_euler","vtk","");
  vtkwriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(xdgf,"solution"));
  vtkwriter.write(0.0,Dune::VTK::appendedraw);


  // <<<11>>> time loop
  Real time = 0.0;
  Real dt = timestep;
  x = 0;
  while (time < Tend)
    {
      // do time step
      dt = osm.apply(time,dt,xold,x);

      // graphics
      xold = x;
      time += dt;
      vtkwriter.write(time,Dune::VTK::appendedraw);

    }
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else {
      if(helper.rank()==0)
        std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
    }

    double Tend;
    double timestep;
    int n;
    int o;
    int choice;

    if (argc==1){
      Tend = 4.8;
      timestep = 0.06;
      n = 64;           // number of cells at one boundary = 64
      o = 1;            // overlap = 1
      choice = 2;       // Implicit Euler
    }
    else if( argc == 6 ){
      sscanf(argv[1],"%lg",&Tend);
      sscanf(argv[2],"%lg",&timestep);
      sscanf(argv[3],"%d",&n);
      sscanf(argv[4],"%d",&o);
      sscanf(argv[5],"%d",&choice);
    }
    else {
      if(helper.rank()==0){
        std::cout << "usage: " << std::endl;
        std::cout << "./test-transport-ccfv <end time> <time step> <elements on a side> <overlap> <method>" << std::endl;
        std::cout << "  <method> = <1|2|3>, 1=stationary, 2=implicit, 3=explicit " << std::endl;
        std::cout << std::endl;
        std::cout << "examples: " << std::endl;
        std::cout << "mpirun -np 4 ./test-transport-ccfv 6.0 1.0 64 1 1" << std::endl;
        std::cout << "mpirun -np 4 ./test-transport-ccfv 6.0 0.06 64 1 2" << std::endl;
        std::cout << "mpirun -np 4 ./test-transport-ccfv 6.0 0.0006 64 1 3" << std::endl;
        std::cout << std::endl;
      }
      return 1;
    }


    // parallel overlapping version
    if (true) {
      const int dim = 2;
      Dune::FieldVector<double,dim> L(1.0);
      std::array<int,dim> N(Dune::fill_array<int,dim>(n));
      std::bitset<dim> periodic(false);
      int overlap=o;
      Dune::YaspGrid<dim> grid(L,N,periodic,overlap,helper.getCommunicator());
      typedef Dune::YaspGrid<dim>::LeafGridView GV;
      GV gv = grid.leafGridView();

      if( choice==1 ){
        std::cout << "\n stationary" << std::endl;
        stationary(gv);
      }

      if( choice==2 ){
        std::cout << "\n implicit_scheme" << std::endl;
        implicit_scheme(gv,Tend,timestep);
      }

      if( choice==3 ){
        std::cout << "\n explicit_scheme" << std::endl;
        explicit_scheme(gv,Tend,timestep);
      }
    }

  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (std::exception &e){
    std::cerr << "STL error: " << e.what() << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
