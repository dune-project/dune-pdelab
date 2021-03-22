// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <ostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/pdelab.hh>

class DummyLOP :
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::NumericalJacobianVolume<DummyLOP>,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  template<typename EG, typename LFSV, typename R>
  void lambda_volume (const EG& eg, const LFSV& lfsv,
                      R& r) const
  {}
  // jacobian of volume term
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu,
                        const X& x, const LFSV& lfsv,
                        M& mat) const
  {}
  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu,
                     const X& x, const LFSV& lfsv,
                     R& r) const
  {}
};

int main(int argc, char** argv) {
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // create classes which NewtonMethod needs for constructor
    using RF = double;
    constexpr int dim=1;
    constexpr int degree=1;
    using Grid = Dune::OneDGrid;
    using DF = Grid::ctype;
    Dune::FieldVector<DF,dim> L; L[0]={1.};
    std::vector<DF> intervals{0,0.5,1};
    Grid grid(intervals);
    using GV = Grid::LeafGridView;
    GV gv=grid.leafGridView();
    using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV,RF,DF,degree>;
    FEM fem(gv);
    using CON = Dune::PDELab::NoConstraints;
    using VBE = Dune::PDELab::ISTL::VectorBackend<>;
    using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
    GFS gfs(gv,fem);
    using LOP = DummyLOP;
    LOP lop;
    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MBE mbe(1);
    using CC = Dune::PDELab::EmptyTransformation;
    using GO = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC>;
    GO go(gfs,gfs,lop,mbe);
    using LS = Dune::PDELab::ISTLBackend_SEQ_UMFPack;
    LS ls(100,0);

    // construct NewtonMethod
    Dune::PDELab::NewtonMethod<GO,LS> newton(go,ls);
    // parameterTree for newton
    Dune::ParameterTree params;
    params["VerbosityLevel"] = "9";
    params["Reduction"]="1e-6";
    params["AbsoluteLimit"]="1e-13";
    params["KeepMatrix"]="false";
    params["UseMaxNorm"]= "true";
    params["MinLinearReduction"]="1e-4";
    params["FixedLinearReduction"]="true";
    params["ReassembleThreshold"]= "1e-5";
    params["LineSearchStrategy"]="hackbuschReusken";
    params["HangingNodeModifications"]="true";
    params["Terminate.MaxIterations"]="11";
    params["Terminate.ForceIteration"]="true";
    params["LineSearch.MaxIterations"] = "8";
    params["LineSearch.DampingFactor"] = "0.6";
    params["LineSearch.AcceptBest"]="true";
    // do not setParameters yet, first check default output

    // output with default parameters
    auto defaultText =
      std::string("NewtonMethod parameters:\n")+
      "Verbosity............... 0\n"+
      "Reduction............... 1e-08\n"+
      "AbsoluteLimit........... 1e-12\n"+
      "KeepMatrix.............. true\n"+
      "UseMaxNorm.............. false\n"+
      "MinLinearReduction...... 0.001\n"+
      "FixedLinearReduction.... false\n"+
      "ReassembleThreshold..... 0\n"+
      "HangingNodeModifications false\n"+
      "Terminate.MaxIterations. 40\n"+
      "Terminate.ForceIteration false\n"+
      "LineSearch.Type........... Hackbusch-Reusken\n"+
      "LineSearch.MaxIterations.. 10\n"+
      "LineSearch.DampingFactor.. 0.5\n"+
      "LineSearch.AcceptBest..... false\n";
    // output with parameters from ParameterTree
    auto testText =
      std::string("NewtonMethod parameters:\n")+
      "Verbosity............... 9\n"+
      "Reduction............... 1e-06\n"+
      "AbsoluteLimit........... 1e-13\n"+
      "KeepMatrix.............. false\n"+
      "UseMaxNorm.............. true\n"+
      "MinLinearReduction...... 0.0001\n"+
      "FixedLinearReduction.... true\n"+
      "ReassembleThreshold..... 1e-05\n"+
      "HangingNodeModifications true\n"+
      "Terminate.MaxIterations. 11\n"+
      "Terminate.ForceIteration true\n"+
      "LineSearch.Type........... Hackbusch-Reusken\n"+
      "LineSearch.MaxIterations.. 8\n"+
      "LineSearch.DampingFactor.. 0.6\n"+
      "LineSearch.AcceptBest..... true\n";

    // variable to store printed output
    std::stringstream buffer0;
    // Redirect std::cout to buffer0
    std::streambuf* prevcoutbuf = std::cout.rdbuf(buffer0.rdbuf());
    // fill buffer0
    newton.printParameters();
    // Restore original buffer before exiting
    std::cout.rdbuf(prevcoutbuf);
    // extract output
    std::string text0 = buffer0.str();
    // compare printed text to default (hard-coded) text, expects result0==0
    int result0 = text0.compare(defaultText);

    // Second test case with changed parameters
    newton.setParameters(params);
    std::stringstream buffer1;
    prevcoutbuf = std::cout.rdbuf(buffer1.rdbuf());
    newton.printParameters();
    std::string text1 = buffer1.str();
    int result1 = text1.compare(testText);
    // Restore original buffer before exiting
    std::cout.rdbuf(prevcoutbuf);

    // check whether default parameters are printed correctly
    if (result0!=0)
      DUNE_THROW(Dune::Exception,"NewtonMethod::printParameters() with default parameters does not match test output!");
    // check whether parameters are printed correctly after they were
    // changed with setParameters. This error hints at problems with
    // setParameters, since printParameters has correct default output.
    if (result1!=0)
      DUNE_THROW(Dune::Exception,"NewtonMethod::printParameters() output does not match test, possibly setParameters problem!");
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