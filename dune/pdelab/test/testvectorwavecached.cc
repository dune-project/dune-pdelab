// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <limits>
#include <string>

#include <dune/common/configparser.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/misc.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>
#include <dune/common/tuples.hh>

#include <dune/grid/common/genericreferenceelements.hh>
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
// Easy access structure for program parameters
//===============================================================

template<class Time_, class DF, class RF, std::size_t dim>
class Config {
  typedef Dune::FieldVector<DF, dim> Domain;

public:
  typedef Time_ Time;

  RF epsilon;
  RF mu;
  Time start;
  Time end;
  Time dt;
  std::string vtkprefix;

  Config(const Dune::ParameterTree &params, DF smallest_edge,
         const std::string &vtkprefix_) :
    epsilon(params.get("epsilon", RF(1))),
    mu(params.get("mu", RF(1))),
    start(params.get("start", Time(0))),
    end(params.get("end", Time(std::sqrt(mu*epsilon)))),
    dt(params.get("dt",
                  Time(smallest_edge*std::sqrt(mu*epsilon/dim)*
                       params.get("dt_stretch", Time(0.35))))),
    vtkprefix(params.get("vtkprefix", vtkprefix_))
  { }
};

template<typename GV>
typename GV::ctype smallest_edge(const GV &gv) {
  typedef typename GV::ctype DF;
  typedef typename GV::template Codim<0>::Iterator Iterator;
  typedef typename GV::template Codim<0>::Geometry Geometry;
  static const std::size_t dim = GV::dimension;
  typedef Dune::FieldVector<DF, GV::dimensionworld> DomainW;

  typedef Dune::GenericReferenceElements<DF, dim> Refelems;
  typedef Dune::GenericReferenceElement<DF, dim> Refelem;

  DF smallest = std::numeric_limits<DF>::infinity();

  const Iterator &end = gv.template end<0>();
  for(Iterator it = gv.template begin<0>(); it != end; ++it) {
    const Refelem & refelem = Refelems::general(it->type());
    const Geometry &geo = it->geometry();
    for(int edge = 0; edge < refelem.size(dim-1); ++edge) {
      const DomainW &cornerpos0 =
        geo.corner(refelem.subEntity(edge, dim-1, 0, dim));
      const DomainW &cornerpos1 =
        geo.corner(refelem.subEntity(edge, dim-1, 1, dim));
      smallest = std::min(smallest, (cornerpos0-cornerpos1).two_norm2());
    }
  }

  return std::sqrt(smallest);
}

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
  Parameters(RF epsilon, RF mu) : Base(1, 1) { }

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

template<class CON, class Config, class GV, class FEM>
void vectorWave(const Config &config, const GV& gv, const FEM& fem)
{
  // constants and types
  typedef typename GV::ctype DF;
  typedef typename FEM::Traits::FiniteElementType::Traits::Basis Basis;
  typedef typename Basis::Traits::RangeField RF;
  static const std::size_t dimRange = Basis::Traits::dimRange;
  typedef typename Config::Time Time;

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
  Params params(config.epsilon, config.mu);

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

  typedef typename GFS::template VectorContainer<RF>::Type V;

  typedef Dune::PDELab::CachedMultiStepGridOperatorSpace<int,DF,RF,GFS,GFS,
    LOPs,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1>, true> MGOS;
  MGOS mgos(msParams, gfs,cg,gfs,cg, LOPRefs(r0, r1, r2));

  typedef Dune::PDELab::VectorWave::CachePolicy<Params, int, DF> CachePolicy;
  mgos.getCache()->setPolicy(Dune::shared_ptr<CachePolicy>
                             (new CachePolicy(params)));

  // make coefficent Vector and initialize it from a function
  {
    Dune::shared_ptr<V> x;

    x.reset(new V(gfs));
    params.setTime(config.start-config.dt);
    Dune::PDELab::interpolate
      ( Dune::PDELab::makeMemberFunctionToGridFunctionAdaptor<RF, dimRange>
        (params, &Params::initialValues, gv),
        gfs,*x);
    mgos.getCache()->setUnknowns(-1, x);

    x.reset(new V(gfs));
    params.setTime(config.start);
    Dune::PDELab::interpolate
      ( Dune::PDELab::makeMemberFunctionToGridFunctionAdaptor<RF, dimRange>
        (params, &Params::initialValues, gv),
        gfs,*x);
    mgos.getCache()->setUnknowns(0, x);
  }

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
  Dune::shared_ptr<Dune::VTKSequenceWriter<GV> > vtkwriter;
  if(config.vtkprefix != "")
    vtkwriter.reset(new Dune::VTKSequenceWriter<GV>(gv,config.vtkprefix,"","",
                                                    Dune::VTK::nonconforming));

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  typedef Dune::PDELab::DiscreteGridFunctionCurl<GFS,V> DGFCurl;
  if(vtkwriter) {
    DGF dgf(gfs,*mgos.getCache()->getUnknowns(-1));
    DGFCurl dgfCurl(gfs,*mgos.getCache()->getUnknowns(-1));
    vtkwriter->addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter->addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGFCurl>(dgfCurl,"curl"));
    vtkwriter->write(config.start-config.dt,Dune::VTK::appendedraw);
    vtkwriter->clear();
  }
  if(vtkwriter) {
    DGF dgf(gfs,*mgos.getCache()->getUnknowns(0));
    DGFCurl dgfCurl(gfs,*mgos.getCache()->getUnknowns(0));
    vtkwriter->addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter->addVertexData
      (new Dune::PDELab::VTKGridFunctionAdapter<DGFCurl>(dgfCurl,"curl"));
    vtkwriter->write(config.start,Dune::VTK::appendedraw);
    vtkwriter->clear();
  }

  DF time = config.start;

  while(time < config.end) {
    Dune::Timer allTimer;
    Dune::Timer subTimer;

    Dune::shared_ptr<const V> xnew = msMethod.apply(time, config.dt);
    time += config.dt;

    if(vtkwriter) {
      subTimer.reset();
      std::cout << "== write output" << std::endl;
      // output grid function with VTKWriter
      DGF dgf(gfs,*xnew);
      DGFCurl dgfCurl(gfs,*xnew);
      vtkwriter->addVertexData
        (new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
      vtkwriter->addVertexData
        (new Dune::PDELab::VTKGridFunctionAdapter<DGFCurl>(dgfCurl,"curl"));
      vtkwriter->write(time,Dune::VTK::appendedraw);
      vtkwriter->clear();
      std::cout << "== write output (" << subTimer.elapsed() << "s)"
                << std::endl;
    }

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

    Dune::ConfigParser paramtree;
    if(argc > 1)
      paramtree.parseFile(argv[1]);

#if HAVE_UG
    // UG Pk 2D test
    {
      static const std::size_t dim = 2;
      typedef Dune::UGGrid<dim> Grid;
      typedef Grid::ctype DF;
      typedef double RF;
      typedef DF Time;
      typedef Dune::FieldVector<DF, dim> Domain;

      const Dune::ParameterTree &myParams = paramtree.hasSub("ug.2d")
        ? paramtree.Dune::ParameterTree::sub("ug.2d")
        : Dune::ParameterTree();

      // make grid
      Dune::shared_ptr<Grid>grid
        (Dune::StructuredGridFactory<Grid>::createSimplexGrid(myParams));

      // get view
      typedef Grid::LeafGridView GV;
      const GV& gv=grid->leafView();

      std::cout << "= Number of elements: " << gv.size(0) << std::endl;

      // get configuration
      Config<Time, DF, RF, dim> config(myParams, smallest_edge(gv),
                                       "vectorwavecached_UG_EdgeS0.5_2D");

      // make finite element map
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
        (config,gv,fem);
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
