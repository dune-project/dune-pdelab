// -*- tab-width: 2; indent-tabs-mode: nil -*-
/** \file

    \brief Solve complex-valued Helmholtz equation.
    \author Philipp Stekl, Marian Piatkowski
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<string>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/float_cmp.hh>
#include<dune/common/typetraits.hh>
#include<dune/grid/yaspgrid.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#else
#error UGERR!
#endif
#include<dune/grid/io/file/vtk.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>

#include<complex>
#define SUPERLU_NTYPE 3 // needed for complex superlu, because superlu is in c where there are no templates
#include<dune/istl/superlu.hh>
// SUPERLU_NTYPE==0 float
// SUPERLU_NTYPE==1 double
// SUPERLU_NTYPE==2 std::complex<float>
// SUPERLU_NTYPE==3 std::complex<double>

#include <dune/istl/umfpack.hh>

#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/backend/istl.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>

#include"testcomplexnumbers-problem.hh"
#include"helmholtzoperator.hh"

//===============================================================
// set up Helmholtz problem and solve it
//===============================================================

template<int k, class GV, class PARAM>
void helmholtz_Qk (const GV& gv, PARAM& param, std::string& errornorm, std::string solver)
{
  using Dune::PDELab::Backend::native;
  using Dune::PDELab::Backend::Native;

  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  const int dim = GV::dimension;
  typedef double R;
  typedef std::complex<R> C;
  typedef C RF;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,C,k> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ConformingDirichletConstraints CON; // constraints class
  typedef Dune::PDELab::istl::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");
  //BCTypeParam bctype; // boundary condition type

  typedef typename GFS::template ConstraintsContainer<C>::Type CC;
  CC cc;
  Dune::PDELab::constraints( param, gfs, cc ); // assemble constraints
  gfs.update();
  std::cout << "constrained dofs=" << cc.size() << " of " << gfs.globalSize() << std::endl;


  // <<<3>>> Make grid operator
  typedef HelmholtzLocalOperator<PARAM > LOP;

  LOP lop(param, 2*k);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,C,C,C,CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop,MBE(std::pow(2*k+1,dim)));

  typedef typename GO::Traits::Domain U;
  C zero(0.);
  U u(gfs,zero);                                      // initial value




  // <<<4>>> Select a linear solver backend
  // typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  //typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS ;
  //typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS;
  // typedef Dune::PDELab::ISTLBackend_SEQ_UMFPack LS;
  // LS ls(5000,0);
  // typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
  // SLP slp(go,ls,u,1e-10);
  // slp.apply();

  //-------------------------------------------------------
  // we do not use the backend, we call solver directly instead


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





    if(solver == "UMFPACK") {
      Dune::UMFPack<ISTLM> solver(native(m), 0);
      r *= -1.0; // need -residual
      //u = r;
      // u = 0;
      Dune::InverseOperatorResult stat;
      solver.apply(native(u),native(r),stat);
    }
    if(solver == "SuperLU") {
      Dune::SuperLU<ISTLM> solver(native(m), 0);
      r *= -1.0; // need -residual
      //u = r;
      // u = 0;
      Dune::InverseOperatorResult stat;
      solver.apply(native(u),native(r),stat);
    }



    Dune::InverseOperatorResult stat;

    if (solver == "GMRESILU0") {
      Dune::SeqILU0<ISTLM,ISTLV,ISTLV> ilu0(native(m),1.0);
      Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(native(m));
      Dune::RestartedGMResSolver<ISTLV> solver(opa, ilu0, 1E-7, 5000, 5000, 0);
      r *= -1.0; // need -residual
      //u = r;
      // u = 0;
      solver.apply(native(u),native(r),stat);
      std::cout<<"Iterations: "<< stat.iterations<<std::endl;
      std::cout<<"Time: "<<  stat.elapsed<< std::endl;
    }
    if (solver == "GMRESILU1") {
      Dune::SeqILUn<ISTLM,ISTLV,ISTLV> ilun(native(m), 1, 1.0);
      Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(native(m));
      Dune::RestartedGMResSolver<ISTLV> solver(opa, ilun, 1E-7, 5000, 5000, 0);
      r *= -1.0; // need -residual
      //u = r;
      // u = 0;
      solver.apply(native(u),native(r),stat);
      std::cout<<"Iterations: "<< stat.iterations<<std::endl;
      std::cout<<"Time: "<<  stat.elapsed<< std::endl;
    }
    if (solver == "GMRESSGS") {
      Dune::SeqSSOR<ISTLM,ISTLV,ISTLV> ssor(native(m), 3, 1.); //ssor with \omega = 1 is SGS (symmetric gauss seidel)
      Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(native(m));
      Dune::RestartedGMResSolver<ISTLV> solver(opa, ssor, 1E-7, 5000, 5000, 0);
      r *= -1.0; // need -residual
      //u = r;
      // u = 0;
      solver.apply(native(u),native(r),stat);
      std::cout<<"Iterations: "<< stat.iterations<<std::endl;
      std::cout<<"Time: "<<  stat.elapsed<< std::endl;
    }
    if (solver == "BiCGSILU0") {
      Dune::SeqILU0<ISTLM,ISTLV,ISTLV> ilu0(native(m),1.0);
      Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(native(m));
      Dune::BiCGSTABSolver<ISTLV> solver(opa,ilu0,1E-7,20000, 0);
      // solve the jacobian system
      r *= -1.0; // need -residual
      //u = r;
      // u = 0;
      solver.apply(native(u),native(r),stat);
      std::cout<<"Iterations: "<< stat.iterations<<std::endl;
      std::cout<<"Time: "<<  stat.elapsed<< std::endl;
    }

    //this solver breaks down at 4000 DOFs by solving the helmholtz eq
    if (solver == "BiCGSILU1") {
      Dune::SeqILUn<ISTLM,ISTLV,ISTLV> ilun(native(m), 1, 1.0);
      Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(native(m));
      Dune::BiCGSTABSolver<ISTLV> solver(opa,ilun,1E-7,20000, 0);
      // solve the jacobian system
      r *= -1.0; // need -residual
      // //u = r;
      // u = 0;
      solver.apply(native(u),native(r),stat);
      std::cout<<"Iterations: "<< stat.iterations<<std::endl;
      std::cout<<"Time: "<<  stat.elapsed<< std::endl;
    }
    if (solver == "BiCGSSGS") {
      Dune::SeqSSOR<ISTLM,ISTLV,ISTLV> ssor(native(m), 3, 1.); //ssor with \omega = 1 is SGS (symmetric gauss seidel)
      Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(native(m));
      Dune::BiCGSTABSolver<ISTLV> solver(opa,ssor,1E-7,20000, 0);
      r *= -1.0; // need -residual
      // //u = r;
      // u = 0;
      solver.apply(native(u),native(r),stat);
      std::cout<<"Iterations: "<< stat.iterations<<std::endl;
      std::cout<<"Time: "<<  stat.elapsed<< std::endl;
    }



  // // Make a real grid function space
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,R,k> FEMr;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEMr,CON,VBE> GFSr;
  FEMr femr(gv);
  GFSr gfsr(gv,femr);

  //create real analytic solution vectors
  using Vr = Dune::PDELab::Backend::Vector<GFSr,R>;
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

  //<<<7>>> graphical output
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, k-1);

  gfsr.name("real");
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfsr,reu);

  gfsr.name("imaginary");
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfsr,imu);

  std::stringstream basename;
  basename << "solvers01_Q" << k;
  vtkwriter.write(basename.str(),Dune::VTK::appendedraw);
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

    std::string errornorm = "L2";
    std::string solver = "SuperLU";

    // sequential version
    if (helper.size()==1) {
      std::vector<double> dof;
      std::vector<double> tsolve; // number of iterations + solver time

      typedef double R;
      typedef std::complex<R> C;

      Dune::FieldVector<double,2> L(1.);
      Dune::array<int,2> N(Dune::fill_array<int,2>(128));
      std::bitset<2> periodic(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(L,N,periodic,overlap);
      typedef Dune::YaspGrid<2>::LeafGridView GV;

      // refine grid
      // grid.globalRefine(level);
      // get leafGridView
      GV gv = grid.leafGridView();

      //define problem
      typedef ParametersPlaneWave<GV, C, R> PARAM;
      //typedef ParametersSphericalWave<GV, C, R> PARAM;

      const double omega = 20.0;
      PARAM param(omega);

      helmholtz_Qk<1,GV,PARAM>(gv,param,errornorm,solver);

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
