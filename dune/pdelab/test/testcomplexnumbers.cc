// -*- tab-width: 2; indent-tabs-mode: nil -*-
/** \file

    \brief Solve complex-valued Helmholtz equation.
    \author Philipp Stekl, Marian Piatkowski

    * This test checks the support for computations with complex field types
    * by solving a complex-valued PDE. The problem arises from the time-harmonic
    * treatment of the wave equation which has been studied in the master thesis
    * by Philipp Stekl.
    * This test uses a slim version of one of the test problems used in this work.
    * It solves the stationary complex-valued Helmholtz equation
    * \f[
    * -\Delta u - \omega^2 u = f
    *  \f]
    * with Robin boundary conditions
    * \f[
    * \nabla u \cdot n - \mathrm{i}\omega u = 0
    * \f]
    * as a first order approximation to the Sommerfeld radiation condition.
    * Thus we also have an example on imposing Robin boundary conditions which are
    * treated as Neumann boundary conditions by allowing the flux also to depend
    * on the solution u.
    *
    * In this test we employ a broad collection of linear solvers which are called
    * either directly from dune-istl or from the backends in dune-pdelab to solve
    * the same problem again and again.
    */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<string>
#include<complex>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/float_cmp.hh>
#include<dune/common/typetraits.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/io/file/vtk.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>
// set number type for SuperLU via preprocessor macro because SuperLU is in C where there are no templates
// here is the summary
// -------------------
// SUPERLU_NTYPE==0 float
// SUPERLU_NTYPE==1 double
// SUPERLU_NTYPE==2 std::complex<float>
// SUPERLU_NTYPE==3 std::complex<double>
#include<dune/istl/umfpack.hh>

#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include"testcomplexnumbers-problem.hh"
#include"helmholtzoperator.hh"

//===============================================================
// set up Helmholtz problem and solve it
//===============================================================

template<int k, class GV, class PARAM>
void helmholtz_Qk (const GV& gv, PARAM& param)
{
  using Dune::PDELab::Backend::native;
  using Dune::PDELab::Backend::Native;

  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  const int dim = GV::dimension;
  typedef double Real;
  typedef std::complex<Real> RF;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,RF,k> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ConformingDirichletConstraints CON; // constraints class
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");
  //BCTypeParam bctype; // boundary condition type

  typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
  CC cc;
  Dune::PDELab::constraints( param, gfs, cc ); // assemble constraints
  gfs.update();
  std::cout << "constrained dofs=" << cc.size() << " of " << gfs.globalSize() << std::endl;


  // <<<3>>> Make grid operator
  typedef HelmholtzLocalOperator<PARAM > LOP;

  LOP lop(param, 2*k);
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop,MBE(std::pow(2*k+1,dim)));

  typedef typename GO::Traits::Domain U;
  RF zero(0.);
  U u(gfs,zero);                                      // initial value

#ifdef USE_ISTL_BACKEND
  // <<<4>>> Select a linear solver backend
  // <<<4.1>>> BCGS SSOR
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR BCGS_SSOR;
  BCGS_SSOR bcgs_ssor;
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,BCGS_SSOR,U> BCGS_SSOR_SLP;
  BCGS_SSOR_SLP bcgs_ssor_slp(go,bcgs_ssor,u,1e-7);
  bcgs_ssor_slp.apply();

  // <<<4.2>>> BCGS ILUn
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILUn BCGS_ILUn;
  BCGS_ILUn bcgs_ilun(1);
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,BCGS_ILUn,U> BCGS_ILUn_SLP;
  BCGS_ILUn_SLP bcgs_ilun_slp(go,bcgs_ilun,1e-7);
  u = 0.0;
  bcgs_ilun_slp.apply(u);

  // <<<4.3>>> BCGS ILU0
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0 BCGS_ILU0;
  BCGS_ILU0 bcgs_ilu0;
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,BCGS_ILU0,U> BCGS_ILU0_SLP;
  BCGS_ILU0_SLP bcgs_ilu0_slp(go,bcgs_ilu0,1e-7);
  u = 0.0;
  bcgs_ilu0_slp.apply(u);

  // <<<4.4>>> GMRes ILU0
  typedef Dune::PDELab::ISTLBackend_SEQ_GMRES_ILU0 GMRES_ILU0;
  GMRES_ILU0 gmres_ilu0(5000,5000);
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,GMRES_ILU0,U> GMRES_ILU0_SLP;
  GMRES_ILU0_SLP gmres_ilu0_slp(go,gmres_ilu0,1e-7);
  u = 0.0;
  gmres_ilu0_slp.apply(u);
#endif // USE_ISTL_BACKEND

#ifdef USE_ISTL

  // do it the manual way by creating all the objects explicitly
  // and calling the linear solver directly instead
  using V = Dune::PDELab::Backend::Vector<GFS,RF>;
  typedef typename GO::Jacobian M;
  using ISTLM = Native<M>;
  using ISTLV = Native<V>;

  M m(go,0.);
  // How well did we estimate the number of entries per matrix row?
  // => print Jacobian pattern statistics
  // std::cout << m.patternStatistics() << std::endl;
  // return;
  go.jacobian(u,m);
  V r(gfs);
  r = 0.0;
  go.residual(u,r);
  r *= -1.0; // need -residual
  V b(r);

  // <<<4>>> Select linear solver
  Dune::InverseOperatorResult stat;
  Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(native(m));

#if HAVE_SUITESPARSE_UMFPACK
  // <<<4.1>>> UMFPack
  std::cout << "=== Using UMFPack as a direct solver" << std::endl;
  Dune::UMFPack<ISTLM> umfpack(native(m), 1);
  umfpack.apply(native(u),native(b),stat);
#endif

#if HAVE_SUPERLU
  // <<<4.2>> SuperLU
  std::cout << "=== Using SuperLU as a direct solver" << std::endl;
  Dune::SuperLU<ISTLM> superlu(native(m), 1);
  u = 0;
  b = r;
  superlu.apply(native(u),native(b),stat);
#endif

  // <<<4.3>>> GMRes ILU0
  std::cout << "=== Using GMRes ILU0" << std::endl;
  Dune::SeqILU0<ISTLM,ISTLV,ISTLV> ilu0(native(m),1.0);
  Dune::RestartedGMResSolver<ISTLV> gmres_ilu0(opa, ilu0, 1E-7, 5000, 5000, 1);
  u = 0;
  b = r;
  gmres_ilu0.apply(native(u),native(b),stat);

  // <<<4.4>>> GMRes ILUn
  std::cout << "=== Using GMRes ILUn" << std::endl;
  Dune::SeqILUn<ISTLM,ISTLV,ISTLV> ilun(native(m), 1, 1.0);
  Dune::RestartedGMResSolver<ISTLV> gmres_ilun(opa, ilun, 1E-7, 5000, 5000, 1);
  u = 0;
  b = r;
  gmres_ilun.apply(native(u),native(b),stat);

  // <<<4.5>>> GMRes SSOR
  std::cout << "=== Using GMRes SSOR" << std::endl;
  Dune::SeqSSOR<ISTLM,ISTLV,ISTLV> ssor(native(m), 3, 1.); //ssor with \omega = 1 is SGS (symmetric gauss seidel)
  Dune::RestartedGMResSolver<ISTLV> gmres_ssor(opa, ssor, 1E-7, 5000, 5000, 1);
  u = 0;
  b = r;
  gmres_ssor.apply(native(u),native(b),stat);

  // <<<4.6>>> BCGS ILU0
  std::cout << "=== Using BiCGSTAB ILU0" << std::endl;
  Dune::BiCGSTABSolver<ISTLV> bcgs_ilu0(opa,ilu0,1E-7,20000, 1);
  u = 0;
  b = r;
  bcgs_ilu0.apply(native(u),native(b),stat);

  // <<<4.7>>> BCGS ILUn
  std::cout << "=== Using BiCGSTAB ILUn" << std::endl;
  Dune::BiCGSTABSolver<ISTLV> bcgs_ilun(opa,ilun,1E-7,20000, 1);
  u = 0;
  b = r;
  bcgs_ilun.apply(native(u),native(b),stat);

  // <<<4.8>>> BCGS SSOR
  std::cout << "=== Using BiCGSTAB SSOR" << std::endl;
  Dune::BiCGSTABSolver<ISTLV> bcgs_ssor(opa,ssor,1E-7,20000, 1);
  u = 0;
  b = r;
  bcgs_ssor.apply(native(u),native(b),stat);
#endif // USE_ISTL

  // Make a real-valued grid function space
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,Real,k> FEMr;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEMr,CON,VBE> GFSr;
  FEMr femr(gv);
  GFSr gfsr(gv,femr);

  //create real-valued analytic solution vectors
  using Vr = Dune::PDELab::Backend::Vector<GFSr,Real>;
  Vr reu(gfsr, 0.0);  // real part u_h
  Vr imu(gfsr, 0.0);  // imag part u_h

  // treat real part and imaginary part separately
  auto reut = reu.begin();
  auto imut = imu.begin();

  for(const auto& x : u) {
    *reut = std::real(x);
    *imut = std::imag(x);

    ++reut;
    ++imut;
  }

  //<<<5>>> graphical output
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, k-1);

  gfsr.name("real");
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfsr,reu);

  gfsr.name("imaginary");
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfsr,imu);

#ifdef USE_ISTL
  std::string filename = "vtk/testcomplexnumbers_istl";
#endif
#ifdef USE_ISTL_BACKEND
  std::string filename = "vtk/testcomplexnumber_istlbackend";
#endif
  vtkwriter.write(filename,Dune::VTK::appendedraw);
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
    else
      {
        if(helper.rank()==0)
          std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
      }

    // sequential version
    if (helper.size()==1) {
      typedef double R;
      typedef std::complex<R> C;

      Dune::FieldVector<double,2> L(1.);
      std::array<int,2> N(Dune::fill_array<int,2>(8));
      std::bitset<2> periodic(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(L,N,periodic,overlap);
      typedef Dune::YaspGrid<2>::LeafGridView GV;

      // get leafGridView
      GV gv = grid.leafGridView();

      //define problem
      typedef ParametersSphericalWave<GV,C,R> PARAM;

      const double omega = 40.0;
      PARAM param(omega);

      helmholtz_Qk<1,GV,PARAM>(gv,param);

    }
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
