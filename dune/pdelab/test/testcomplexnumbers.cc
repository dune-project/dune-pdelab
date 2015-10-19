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

#if 0
#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#endif
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/backend/istl.hh>
#if 0
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/instationary/onestep.hh> //Filename helper
#include <dune/common/parametertreeparser.hh>
#endif

#if 0
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/boilerplate/pdelab.hh>
#endif


#if 0
#include"parameters_planewave.hh"
#include"parameters_sphericalwave.hh"
#include"helmholtz_bcextension.hh"
#include"helmholtz_bcanalytic.hh"
#endif

#include"testcomplexnumbers-problem.hh"
#include"helmholtzoperator.hh"

//===============================================================
// set up Helmholtz problem and solve it
//===============================================================

template<int k, class GV, class PARAM>
void helmholtz_Qk (const GV& gv, PARAM& param, std::string& errornorm, std::string solver)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
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

  LOP lop( param, 2*k );
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  // Structured 2D grid, Q1 finite elements => 9-point stencil / Q2 => 25
  MBE mbe(k == 1 ? 9 : 25);

  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,C,C,C,CC,CC> GO;
  //GO go(gfs,gfs,lop,mbe);
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // How well did we estimate the number of entries per matrix row?
  // => print Jacobian pattern statistics
  //typename GO::Traits::Jacobian jac(go);
  //std::cout << jac.patternStatistics() << std::endl;



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


  typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type V;
  typedef typename GO::Jacobian M;

  typedef typename M::BaseT ISTLM;
  typedef typename V::BaseT ISTLV;


  M m(go,0.);
  //std::cout << m.patternStatistics() << std::endl;
  go.jacobian(u,m);
  V r(gfs);
  r = 0.0;
  go.residual(u,r);





    if(solver == "UMFPACK") {
      Dune::UMFPack<ISTLM> solver(m.base(), 0);
      r *= -1.0; // need -residual
      //u = r;
      // u = 0;
      Dune::InverseOperatorResult stat;
      solver.apply(u.base(),r.base(),stat);
    }
    if(solver == "SuperLU") {
      Dune::SuperLU<ISTLM> solver(m.base(), 0);
      r *= -1.0; // need -residual
      //u = r;
      // u = 0;
      Dune::InverseOperatorResult stat;
      solver.apply(u.base(),r.base(),stat);
    }



    Dune::InverseOperatorResult stat;

    if (solver == "GMRESILU0") {
      Dune::SeqILU0<ISTLM,ISTLV,ISTLV> ilu0(m.base(),1.0);
      Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(m.base());
      Dune::RestartedGMResSolver<ISTLV> solver(opa, ilu0, 1E-7, 5000, 5000, 0);
      r *= -1.0; // need -residual
      //u = r;
      // u = 0;
      solver.apply(u.base(),r.base(),stat);
      std::cout<<"Iterations: "<< stat.iterations<<std::endl;
      std::cout<<"Time: "<<  stat.elapsed<< std::endl;
    }
    if (solver == "GMRESILU1") {
      Dune::SeqILUn<ISTLM,ISTLV,ISTLV> ilun(m.base(), 1, 1.0);
      Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(m.base());
      Dune::RestartedGMResSolver<ISTLV> solver(opa, ilun, 1E-7, 5000, 5000, 0);
      r *= -1.0; // need -residual
      //u = r;
      // u = 0;
      solver.apply(u.base(),r.base(),stat);
      std::cout<<"Iterations: "<< stat.iterations<<std::endl;
      std::cout<<"Time: "<<  stat.elapsed<< std::endl;
    }
    if (solver == "GMRESSGS") {
      Dune::SeqSSOR<ISTLM,ISTLV,ISTLV> ssor(m.base(), 3, 1.); //ssor with \omega = 1 is SGS (symmetric gauss seidel)
      Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(m.base());
      Dune::RestartedGMResSolver<ISTLV> solver(opa, ssor, 1E-7, 5000, 5000, 0);
      r *= -1.0; // need -residual
      //u = r;
      // u = 0;
      solver.apply(u.base(),r.base(),stat);
      std::cout<<"Iterations: "<< stat.iterations<<std::endl;
      std::cout<<"Time: "<<  stat.elapsed<< std::endl;
    }
    if (solver == "BiCGSILU0") {
      Dune::SeqILU0<ISTLM,ISTLV,ISTLV> ilu0(m.base(),1.0);
      Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(m.base());
      Dune::BiCGSTABSolver<ISTLV> solver(opa,ilu0,1E-7,20000, 0);
      // solve the jacobian system
      r *= -1.0; // need -residual
      //u = r;
      // u = 0;
      solver.apply(u.base(),r.base(),stat);
      std::cout<<"Iterations: "<< stat.iterations<<std::endl;
      std::cout<<"Time: "<<  stat.elapsed<< std::endl;
    }

    //this solver breaks down at 4000 DOFs by solving the helmholtz eq
    if (solver == "BiCGSILU1") {
      Dune::SeqILUn<ISTLM,ISTLV,ISTLV> ilun(m.base(), 1, 1.0);
      Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(m.base());
      Dune::BiCGSTABSolver<ISTLV> solver(opa,ilun,1E-7,20000, 0);
      // solve the jacobian system
      r *= -1.0; // need -residual
      // //u = r;
      // u = 0;
      solver.apply(u.base(),r.base(),stat);
      std::cout<<"Iterations: "<< stat.iterations<<std::endl;
      std::cout<<"Time: "<<  stat.elapsed<< std::endl;
    }
    if (solver == "BiCGSSGS") {
      Dune::SeqSSOR<ISTLM,ISTLV,ISTLV> ssor(m.base(), 3, 1.); //ssor with \omega = 1 is SGS (symmetric gauss seidel)
      Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(m.base());
      Dune::BiCGSTABSolver<ISTLV> solver(opa,ssor,1E-7,20000, 0);
      r *= -1.0; // need -residual
      // //u = r;
      // u = 0;
      solver.apply(u.base(),r.base(),stat);
      std::cout<<"Iterations: "<< stat.iterations<<std::endl;
      std::cout<<"Time: "<<  stat.elapsed<< std::endl;
    }



#if 0
  // // // compare helmholtz solution to analytic helmholtz solution
  U ua(gfs, zero);                                      // initial value
  typedef BCAnalytic<PARAM> Ga;                      // boundary value = extension
  Ga ga(gv, param);
  Dune::PDELab::interpolate(ga,gfs,ua);                // interpolate coefficient vector


  // // Make a real grid function space
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,R,k> FEMr;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEMr,CON,VBE> GFSr;
  FEMr femr(gv);
  GFSr gfsr(gv,femr);

  //create real analytic solution vectors
  typename Dune::PDELab::BackendVectorSelector<GFSr,R>::Type reu(gfsr, 0.);  // real part u_h
  typename Dune::PDELab::BackendVectorSelector<GFSr,R>::Type imu(gfsr, 0.);  // imag part u_h
  typename Dune::PDELab::BackendVectorSelector<GFSr,R>::Type reua(gfsr, 0.); // real part analytic solution
  typename Dune::PDELab::BackendVectorSelector<GFSr,R>::Type imua(gfsr, 0.); // imega part analytic solution
  // typename Dune::PDELab::BackendVectorSelector<GFSr,R>::Type reue(gfsr, 0.); // real part error
  // typename Dune::PDELab::BackendVectorSelector<GFSr,R>::Type imue(gfsr, 0.); // imag part error

  typename Dune::PDELab::BackendVectorSelector<GFSr,R>::Type::iterator reut = reu.begin();
  typename Dune::PDELab::BackendVectorSelector<GFSr,R>::Type::iterator imut = imu.begin();
  typename Dune::PDELab::BackendVectorSelector<GFSr,R>::Type::iterator reuat = reua.begin();
  typename Dune::PDELab::BackendVectorSelector<GFSr,R>::Type::iterator imuat = imua.begin();
  // typename Dune::PDELab::BackendVectorSelector<GFSr,R>::Type::iterator reuet = reue.begin();
  // typename Dune::PDELab::BackendVectorSelector<GFSr,R>::Type::iterator imuet = imue.begin();


  typename U::const_iterator ut = u.begin();
  typename U::const_iterator uat = ua.begin();

  // //compute error i.e. |ua - u| separated into real and imag part
  while( ut != u.end() ) {

    *reut = std::real(*ut);
    *imut = std::imag(*ut);
    *reuat = std::real(*uat);
    *imuat = std::imag(*uat);
    // *reuet = std::abs(*reut - *reuat) ;
    // *imuet = std::abs(*imut - *imuat) ;

    // std::cout<<h<<"\t"<<*reut - (*reuat)<<std::endl;
    //   ++h;

    ++ut;
    ++uat;

    ++reut;
    ++imut;
    ++reuat;
    ++imuat;
    // ++reuet;
    // ++imuet;

  }
#endif


#if 0
  //<<<7>>> graphical output

  // output u_h
  //Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,k == 1 ? 0 : 3);
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, 0);


  gfsr.name("real");
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfsr,reu);

  gfsr.name("imaginary");
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfsr,imu);

  gfsr.name("real_analytic");
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfsr,reua);
  gfsr.name("imaginary_analytic");
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfsr,imua);


  // gfsr.name("real_error");
  // Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfsr,reue);
  // gfsr.name("imaginary_error");
  // Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfsr,imue);


  std::stringstream basename;
  basename << "solvers01_Q" << k;
  vtkwriter.write(basename.str(),Dune::VTK::appendedraw);
#endif





 // //  // compute  DISCRETIZATION L2 error
 // //  /* solution */
 //  if(errornorm == "L2") {
 //    // define discrete grid function for solution
 //    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
 //    // DGF udgf(gfs, u);

 //    typedef BCAnalytic<PARAM> ES;
 //    ES es(gv, param);
 //    U esh(gfs, zero);
 //    Dune::PDELab::interpolate(es,gfs,esh);
 //    DGF esdgf(gfs, esh);

 //    typedef DifferenceSquaredAdapter<DGF,ES> DifferenceSquared;
 //    DifferenceSquared difference(esdgf, es);

 //    typename DifferenceSquared::Traits::RangeType l2normsquared(0.0);
 //    Dune::PDELab::integrateGridFunction(difference,l2normsquared,10);

 //    std::cout<<"L2DiscrError: "<<std::sqrt(std::abs(l2normsquared[0]))<<std::endl;
 //  }



 // //  // compute  POLLUTION L2 error
 // //  /* solution */
 //  if(errornorm == "L2") {
 //    // define discrete grid function for solution
 //    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
 //    DGF udgf(gfs, u);

 //    typedef BCAnalytic<PARAM> ES;
 //    ES es(gv, param);
 //    U esh(gfs, zero);
 //    Dune::PDELab::interpolate(es,gfs,esh);
 //    DGF esdgf(gfs, esh);


 //    typedef DifferenceSquaredAdapter<DGF,DGF> DifferenceSquared;
 //    DifferenceSquared difference(udgf,esdgf);

 //    typename DifferenceSquared::Traits::RangeType l2normsquared(0.0);
 //    Dune::PDELab::integrateGridFunction(difference,l2normsquared,10);

 //    std::cout<<"L2PolError: "<<std::sqrt(std::abs(l2normsquared[0]))<<std::endl;

 //  }




#if 0
 //  // compute L2 error
 //  /* solution */
  if(errornorm == "L2") {
    // define discrete grid function for solution
    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
    DGF udgf(gfs, u);

    typedef BCAnalytic<PARAM> ES;
    ES es(gv, param);

    typedef DifferenceSquaredAdapter<DGF,ES> DifferenceSquared;
    DifferenceSquared difference(udgf,es);

    typename DifferenceSquared::Traits::RangeType l2normsquared(0.0);
    Dune::PDELab::integrateGridFunction(difference,l2normsquared,10);

    std::cout<<"L2Error: "<<std::sqrt(std::abs(l2normsquared[0]))<<std::endl;
  }



  // compute H1 semi norm error
  /* solution */
  if(errornorm == "H1") {
    //discrete function gradient
    typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS,U> DGFG;
    DGFG udgfg(gfs, u);

    //gradient of exact solution
    typedef BCAnalyticGrad<PARAM> ESG;
    ESG esg(gv, param);

    // difference of gradient
    typedef DifferenceSquaredAdapter<DGFG,ESG> DifferenceSquaredAdapterg;
    DifferenceSquaredAdapterg differenceg(udgfg,esg);

    typename DifferenceSquaredAdapterg::Traits::RangeType l2normsquared(0.0);
    typename DifferenceSquaredAdapter<DGFG,ESG>::Traits::RangeType normsquared(0.0);
    Dune::PDELab::integrateGridFunction(differenceg,normsquared,10);

    std::cout<<"H1Error: "<<std::sqrt(std::abs(l2normsquared[0]))<<std::endl;

  }
#endif



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

#if 0
    if (argc!=2)
      {
        if (helper.rank()==0)
          std::cout << "usage: ./error <cfg-file>" << std::endl;
        return 1;
      }
#endif

#if 0
    // Parse configuration file.
    std::string config_file(argv[1]);
    Dune::ParameterTree configuration;
    Dune::ParameterTreeParser parser;

    try{
      parser.readINITree( config_file, configuration );
    }
    catch(...){
      std::cerr << "Could not read config file \""
                << config_file << "\"!" << std::endl;
      exit(1);
    }

    int level = configuration.get<int>("grid.level");
    int polynomialdegree = configuration.get<int>("problem_parameter.polynomialdegree");
    double omega = configuration.get<double>("problem_parameter.omega") ;
    std::string errornorm = configuration.get<std::string>("norm.errornorm");
    std::string solver = configuration.get<std::string>("solvers.solver");
#endif
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

      const double omega = 80.0;
      PARAM param(omega);


#if 0
      if(polynomialdegree == 1) {
        helmholtz_Qk<1,GV,PARAM>(gv, param,  errornorm, solver);
      }

      if(polynomialdegree == 2) {
        helmholtz_Qk<2,GV,PARAM>(gv, param, errornorm, solver);
      }
#endif
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
