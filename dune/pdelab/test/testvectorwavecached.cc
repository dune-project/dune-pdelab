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
#include <utility>

#include <dune/common/array.hh>
#include <dune/common/configparser.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/misc.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>
#include <dune/common/tuples.hh>

#include <dune/grid/alugrid.hh>
#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/common/vertexorder.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/finiteelementmap/edges0.5fem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/function/memberadaptor.hh>
#include <dune/pdelab/gridfunctionspace/constraints.hh>
#include <dune/pdelab/gridfunctionspace/dofinfo.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/printvector.hh>
#include <dune/pdelab/linearsolver/stationarymatrix.hh>
#include <dune/pdelab/localoperator/vectorwave.hh>
#include <dune/pdelab/multistep/gridoperatorspace.hh>
#include <dune/pdelab/multistep/method.hh>
#include <dune/pdelab/multistep/parameter.hh>

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
  struct {
    std::string prefix;
    Time min_interval;
  } vtkoutput;
  struct {
    bool dof_positions;
    int coord_precision;
  } debug;

  Config(const Dune::ParameterTree &params, DF smallest_edge,
         const std::string &vtkprefix_) :
    epsilon(params.get("epsilon", RF(1))),
    mu(params.get("mu", RF(1))),
    start(params.get("start", Time(0))),
    end(params.get("end", Time(std::sqrt(mu*epsilon)))),
    dt(params.get("dt",
                  Time(smallest_edge*std::sqrt(mu*epsilon/dim)*
                       params.get("dt_stretch", Time(0.35)))))
  {
    vtkoutput.prefix = params.get("vtkoutput.prefix", vtkprefix_);
    vtkoutput.min_interval = params.get("vtkoutput.min_interval", Time(0));

    debug.dof_positions = params.get("debug.dof_positions", false);
    debug.coord_precision = params.get("debug.coord_precision", int(3));
  }
};

template<typename GV>
typename GV::ctype smallest_edge(const GV &gv) {
  typedef typename GV::ctype DF;
  typedef typename GV::template Codim<0>::
    template Partition<Dune::Interior_Partition>::Iterator Iterator;
  typedef typename GV::template Codim<0>::Geometry Geometry;
  static const std::size_t dim = GV::dimension;
  typedef Dune::FieldVector<DF, GV::dimensionworld> DomainW;

  typedef Dune::GenericReferenceElements<DF, dim> Refelems;
  typedef Dune::GenericReferenceElement<DF, dim> Refelem;

  DF smallest = std::numeric_limits<DF>::infinity();

  const Iterator &end = gv.template end<0, Dune::Interior_Partition>();
  for(Iterator it = gv.template begin<0, Dune::Interior_Partition>();
      it != end; ++it)
  {
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

  return std::sqrt(gv.comm().min(smallest));
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

  std::size_t bcType
  ( const Dune::PDELab::IntersectionGeometry<typename GV::Intersection> &is,
    const Dune::FieldVector<typename GV::ctype, GV::dimension-1> &xl) const
  {
    if(is.boundary() && ! is.neighbor())
      return 1;
    else
      return 0;
  }

};


//===============================================================
// Problem setup and solution
//===============================================================

template<class Config, class GV, class FEM>
void vectorWave(const Config &config, const GV& gv, const FEM& fem)
{
  // constants and types
  typedef typename GV::ctype DF;
  static const std::size_t dimw = GV::dimensionworld;
  typedef Dune::FieldVector<DF, dimw> DomainW;
  typedef typename FEM::Traits::FiniteElementType::Traits::Basis Basis;
  typedef typename Basis::Traits::RangeField RF;
  static const std::size_t dimRange = Basis::Traits::dimRange;
  typedef typename Config::Time Time;

  // make constraints evaluator
  typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints CE;
  CE ce;

  // make function space and constraints evaluator
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CE,
    Dune::PDELab::ISTLVectorBackend<1> > GFS;
  GFS gfs(gv, fem, ce);
  ce.compute_ghosts(gfs); // con stores indices of ghost dofs

  if(config.debug.dof_positions) {
    typedef typename GFS::template VectorContainer<std::size_t>::Type IntV;
    IntV dummy(gfs);
    for(std::size_t i = 0; i < gfs.size(); ++i)
      IntV::Backend::access(dummy, i) = i;
    Dune::PDELab::printVector(std::cout, gfs, dummy, "dof positions",
                              config.debug.coord_precision);
  }

  typedef Parameters<GV, RF, Time> Params;
  Params params(config.epsilon, config.mu);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<RF>::Type C;
  C cg;
  cg.clear();
  Dune::PDELab::constraints
    ( Dune::PDELab::make2ArgsMemberFunctionToBoundaryGridFunctionAdaptor
        <std::size_t, 1>(params, &Params::bcType, gv),
      gfs, cg);

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
  typedef Dune::PDELab::StationaryMatrixLinearSolver<MGOS, LS, RF> PDESOLVER;
  PDESOLVER pdesolver(mgos, ls, 1e-9);

  // <<<8>>> time-stepper
  Dune::PDELab::MultiStepMethod<DF,MGOS,PDESOLVER,V,V>
    msMethod(msParams,mgos,pdesolver);
  msMethod.setVerbosityLevel(2);

  // output grid function with VTKWriter
  Dune::shared_ptr<Dune::VTKSequenceWriter<GV> > vtkwriter;
  if(config.vtkoutput.prefix != "")
    vtkwriter.reset(new Dune::VTKSequenceWriter<GV>(gv,
                                                    config.vtkoutput.prefix,
                                                    "","",
                                                    Dune::VTK::nonconforming));

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  typedef Dune::PDELab::DiscreteGridFunctionCurl<GFS,V> DGFCurl;
  if(vtkwriter && config.dt >= config.vtkoutput.min_interval) {
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

  Time time = config.start;
  Time last_planned_vtkoutput = config.start;

  while(time < config.end) {
    Dune::Timer allTimer;
    Dune::Timer subTimer;

    Dune::shared_ptr<const V> xnew = msMethod.apply(time, config.dt);
    time += config.dt;

    if(vtkwriter && time >= last_planned_vtkoutput) {
      subTimer.reset();
      if(gv.comm().rank() == 0)
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

      last_planned_vtkoutput += config.vtkoutput.min_interval;

      if(gv.comm().rank() == 0)
        std::cout << "== write output (" << subTimer.elapsed() << "s)"
                  << std::endl;
    }

    if(gv.comm().rank() == 0)
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

#if HAVE_ALUGRID
    {
      static const std::size_t dim = 3;
      typedef Dune::ALUSimplexGrid<dim, dim> Grid;
      typedef Grid::ctype DF;
      typedef double RF;
      typedef DF Time;
      typedef Dune::FieldVector<DF, dim> Domain;

      const Dune::ParameterTree &myParams = paramtree.hasSub("alu.3d")
        ? paramtree.Dune::ParameterTree::sub("alu.3d")
        : Dune::ParameterTree();

      // make grid
      std::pair<Domain, Domain> bbox;
      bbox.first = myParams.get("bbox.lower", Domain(0));
      bbox.second = myParams.get("bbox.upper", Domain(1));
      Dune::array<unsigned, dim> elements;
      std::fill(elements.begin(), elements.end(), 1);
      elements = myParams.get("elements", elements);
      Dune::shared_ptr<Grid>grid
        ( Dune::StructuredGridFactory<Grid>::createSimplexGrid
          ( bbox.first, bbox.second, elements));
      grid->loadBalance();
      grid->globalRefine(myParams.get("global_refines", 0));

      // get view
      typedef Grid::LeafGridView GV;
      const GV& gv=grid->leafView();

      {
        int sizes[gv.comm().size()];
        int overlapsizes[gv.comm().size()];
        int ghostsizes[gv.comm().size()];

        int mysize = gv.size(0);
        gv.comm().gather(&mysize, sizes, 1, 0);
        mysize = std::distance(gv.begin<0, Dune::Overlap_Partition>(),
                               gv.end<0, Dune::Overlap_Partition>());
        gv.comm().gather(&mysize, overlapsizes, 1, 0);
        mysize = std::distance(gv.begin<0, Dune::Ghost_Partition>(),
                               gv.end<0, Dune::Ghost_Partition>());
        gv.comm().gather(&mysize, ghostsizes, 1, 0);

        if(gv.comm().rank() == 0) {
          int allsize = 0, alloverlapsize = 0, allghostsize = 0;
          for(int i = 0; i < gv.comm().size(); ++i) {
            allsize += sizes[i];
            alloverlapsize += overlapsizes[i];
            allghostsize += ghostsizes[i];
          }
          std::cout << "= Total number of elements: " << allsize << " "
                    << "(interior: "
                    << allsize-alloverlapsize-allghostsize << ", overlap: "
                    << alloverlapsize << ", ghost: " << allghostsize << ")"
                    << std::endl;
          for(int i = 0; i < gv.comm().size(); ++i)
            std::cout << "= Number of elements (rank " << i << "): "
                      << sizes[i] << " (interior: "
                      << sizes[i]-overlapsizes[i]-ghostsizes[i] << ", "
                      << "overlap: " << overlapsizes[i] << ", ghost: "
                      << ghostsizes[i] << ")" << std::endl;
        }
      }
      DF smallest = smallest_edge(gv);
      if(gv.comm().rank() == 0)
        std::cout << "= Smallest edge: " << smallest << std::endl;

      // get configuration
      Config<Time, DF, RF, dim> config(myParams, smallest,
                                       "vectorwavecached_ALU_EdgeS0.5_3D");

      // make finite element map
      typedef Dune::PDELab::VertexOrderByIdFactory<Grid::GlobalIdSet>
        VOFactory;
      VOFactory voFactory(grid->globalIdSet());
      typedef Dune::PDELab::EdgeS0_5FiniteElementMap<
        Grid::Codim<0>::Geometry, VOFactory, RF
        > FEM;
      FEM fem(voFactory);

      // solve problem
      vectorWave(config,gv,fem);
    }
#endif // HAVE_ALUGRID

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
