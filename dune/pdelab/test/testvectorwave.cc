// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <iostream>
#include <string>

#include <dune/common/array.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/misc.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>
#include <dune/common/tuples.hh>

#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/common/vertexorder.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/finiteelement/edges0.5.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/finiteelementmap/global.hh>
#include <dune/pdelab/function/const.hh>
#include <dune/pdelab/function/memberadaptor.hh>
#include <dune/pdelab/gridfunctionspace/constraints.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/localoperator/vectorwave.hh>
#include <dune/pdelab/multistep/gridoperatorspace.hh>
#include <dune/pdelab/multistep/method.hh>
#include <dune/pdelab/multistep/parameter.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

//===============================================================
//===============================================================
// Solve 1D homogenous wave equation
//   \partial_t^2(eps*E) + curl((1/mu) * curl E) = 0 in \Omega,
//                                         n × E = 0 on \partial\Omega
//===============================================================
//===============================================================

//===============================================================
// Define parameter functions epsilon, mu and initial values
//===============================================================

template<class GV, class RF = double, class Time = double>
class Parameters :
  public Dune::PDELab::VectorWave::ConstantParameters<GV, RF, Time>
{
  typedef Dune::PDELab::VectorWave::ConstantParameters<GV, RF, Time> Base;

  Time time;

public:
  Parameters() : Base(1, 1) { }

  void setTime(Time time_) { Base::setTime(time); time = time_; }

  void initialValues(const typename Base::Element &e,
                     const typename Base::Domain &xl,
                     typename Base::Range& y) const
  {
    const typename Base::Domain &xg = e.geometry().global(xl);
    y = 0;
    y[1] = std::exp(-Dune::SQR(xg[0]-.5+time)/(2*Dune::SQR(0.05)));
    y[0] = std::exp(-Dune::SQR(xg[1]-.5+time)/(2*Dune::SQR(0.05)));
  }
};


//===============================================================
// Problem setup and solution
//===============================================================

template<class CON, class Time, class GV, class FEM>
void vectorWave(const GV& gv, const FEM& fem, Time dt, std::size_t steps,
                std::string filename)
{
  // constants and types
  typedef typename GV::ctype DF;
  typedef typename FEM::Traits::FiniteElementType::Traits::Basis Basis;
  typedef typename Basis::Traits::RangeField RF;
  static const std::size_t dimRange = Basis::Traits::dimRange;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS;
  GFS gfs(gv,fem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<RF>::Type C;
  C cg;
  cg.clear();
  typedef Dune::PDELab::ConstBoundaryGridFunction<GV, std::size_t> BType;
  BType b(gv, 1);
  Dune::PDELab::constraints(b,gfs,cg);

  typedef Parameters<GV, RF, Time> Params;
  Params params;

  // make coefficent Vector and initialize it from a function
  typedef typename GFS::template VectorContainer<RF>::Type V;
  std::deque<Dune::shared_ptr<V> > oldvalues;
  oldvalues.push_front(Dune::shared_ptr<V>(new V(gfs)));
  params.setTime(-dt);
  Dune::PDELab::interpolate
    ( Dune::PDELab::makeMemberFunctionToGridFunctionAdaptor<RF, dimRange>
      (params, &Params::initialValues, gv),
      gfs,*oldvalues.front());
  Dune::PDELab::set_constrained_dofs(cg, 0, *oldvalues.front());
  oldvalues.push_front(Dune::shared_ptr<V>(new V(gfs)));
  params.setTime(0);
  Dune::PDELab::interpolate
    ( Dune::PDELab::makeMemberFunctionToGridFunctionAdaptor<RF, dimRange>
      (params, &Params::initialValues, gv),
      gfs,*oldvalues.front());
  Dune::PDELab::set_constrained_dofs(cg, 0, *oldvalues.front());

  // make parameters
  Dune::PDELab::CentralDifferencesParameters<DF> msParams;

  typedef Dune::PDELab::VectorWave::R0<Params> R0;
  typedef Dune::PDELab::VectorWave::R1<Params> R1;
  typedef Dune::PDELab::VectorWave::R2<Params> R2;

  R0 r0(params);
  R1 r1(params);
  R2 r2(params);

  typedef Dune::tuple<R0, R1, R2> LOPs;
  typedef Dune::tuple<R0&, R1&, R2&> LOPRefs;

  typedef Dune::PDELab::MultiStepGridOperatorSpace<DF,V,GFS,GFS,
    LOPs,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1>, true> MGOS;
  MGOS mgos(msParams, gfs,cg,gfs,cg, LOPRefs(r0, r1, r2));

// #if HAVE_SUPERLU
//   typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
//   LS ls(false);
// #else
  typedef Dune::PDELab::ISTLBackend_NOVLP_CG_NOPREC<GFS> LS;
  LS ls(gfs,5000,2);
  //#endif

  // <<<7>>> make Newton for time-dependent problem
  typedef Dune::PDELab::StationaryLinearProblemSolver<MGOS, LS, V> PDESOLVER;
  PDESOLVER pdesolver(mgos, ls, 1e-9);

  // <<<8>>> time-stepper
  Dune::PDELab::MultiStepMethod<DF,MGOS,PDESOLVER,V,V>
    msMethod(msParams,mgos,pdesolver);
  msMethod.setVerbosityLevel(2);

  // output grid function with VTKWriter
  Dune::VTKSequenceWriter<GV> vtkwriter(gv,filename,"","",
                                        Dune::VTK::nonconforming);

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  typedef Dune::PDELab::DiscreteGridFunctionCurl<GFS,V> DGFCurl;
  {
    DGF dgf(gfs,*oldvalues[1]);
    DGFCurl dgfCurl(gfs,*oldvalues[1]);
    vtkwriter.addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGFCurl>(dgfCurl,"curl"));
    vtkwriter.write(-dt,Dune::VTK::appendedraw);
    vtkwriter.clear();
  }
  {
    DGF dgf(gfs,*oldvalues[0]);
    DGFCurl dgfCurl(gfs,*oldvalues[0]);
    vtkwriter.addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGFCurl>(dgfCurl,"curl"));
    vtkwriter.write(0,Dune::VTK::appendedraw);
    vtkwriter.clear();
  }

  DF time = 0;

  for(unsigned step = 0; step < steps; ++step) {
    Dune::Timer allTimer;
    Dune::Timer subTimer;

    subTimer.reset();
    std::cout << "== setup result vector" << std::endl;
    Dune::shared_ptr<V> xnew(new V(*oldvalues.front()));
    std::cout << "== setup result vector (" << subTimer.elapsed() << "s)"
              << std::endl;

    dt = msMethod.apply(time, dt, oldvalues, *xnew);
    time += dt;

    oldvalues.pop_back();
    oldvalues.push_front(xnew);

    subTimer.reset();
    std::cout << "== write output" << std::endl;
    // output grid function with VTKWriter
    DGF dgf(gfs,*xnew);
    DGFCurl dgfCurl(gfs,*xnew);
    vtkwriter.addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGFCurl>(dgfCurl,"curl"));
    vtkwriter.write(time,Dune::VTK::appendedraw);
    vtkwriter.clear();
    std::cout << "== write output (" << subTimer.elapsed() << "s)"
              << std::endl;

    std::cout << "= time step total time: " << allTimer.elapsed() << "s"
              << std::endl;
  }
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

#if HAVE_UG
    // UG Pk 2D test
    {
      static const std::size_t elems = 32;
      static const std::size_t time_strech = 4;

      // make grid
      static const std::size_t dim = 2;
      typedef Dune::UGGrid<dim> Grid;
      typedef Grid::ctype DF;
      typedef Dune::FieldVector<DF, dim> Domain;
      Dune::array<unsigned, dim> elemCount;
      std::fill(elemCount.begin(), elemCount.end(), elems);

      Dune::shared_ptr<Grid>grid
        (Dune::StructuredGridFactory<Grid>::createSimplexGrid(Domain(0),
                                                              Domain(1),
                                                              elemCount));

      // get view
      typedef Grid::LeafGridView GV;
      const GV& gv=grid->leafView();

      // make finite element map
      typedef double RF;
      typedef Dune::PDELab::EdgeS0_5FiniteElementFactory<
        Grid::Codim<0>::Geometry, RF
        > FEFactory;
      FEFactory feFactory;
      typedef Dune::PDELab::VertexOrderByIdFactory<Grid::GlobalIdSet>
        VOFactory;
      VOFactory voFactory(grid->globalIdSet());
      typedef Dune::PDELab::GeometryVertexOrderFiniteElementMap<
        FEFactory, VOFactory
        > FEM;
      FEM fem(feFactory, voFactory);

      // solve problem
      vectorWave<Dune::PDELab::ConformingDirichletConstraints>
        (gv,fem,1.0/time_strech/elems,time_strech*elems,"vectorwave_UG_EdgeS0.5_2D");
    }
#endif

    // test passed
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    throw;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    throw;
  }
}
