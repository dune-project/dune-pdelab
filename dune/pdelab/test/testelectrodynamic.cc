// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#define DISABLE_LOCAL_INTERFACE

#include <iostream>
#include <string>

#include <unistd.h>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/smartpointer.hh>

#include <dune/grid/albertagrid/agrid.hh>

#ifdef HAVE_ALUGRID
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#endif

#ifdef HAVE_ALBERTA
#include <dune/grid/albertagrid/dgfparser.hh>
#endif

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>

#include "../backend/istlmatrixbackend.hh"
#include "../backend/istlsolverbackend.hh"
#include "../backend/istlvectorbackend.hh"
#include "../common/function.hh"
#include "../common/geometrywrapper.hh"
#include "../finiteelementmap/conformingconstraints.hh"
#include "../finiteelementmap/edges03dfem.hh"
#include "../gridfunctionspace/constraints.hh"
#include "../gridfunctionspace/gridfunctionspace.hh"
#include "../gridfunctionspace/gridfunctionspaceutilities.hh"
#include "../gridfunctionspace/interpolate.hh"
#include "../gridoperatorspace/gridoperatorspace.hh"
#include "../localoperator/electrodynamic.hh"

#include "gnuplotgraph.hh"
#include "gridexamples.hh"
#include "l2difference.hh"
#include "physicalconstants.hh"
#include "probe.hh"
#include "resonatorsolution.hh"

//===============================================================
//===============================================================
// Solve the simple electrodynamic wave equation
//    rot(1/mu * rot E) + eps * \partial_t^2 E = 0 in \Omega, 
//                                       n x E = 0 on \partial\Omega_D
//===============================================================
//===============================================================

//
//  CONFIGURATION
//

// default limit for the total convergence
// (alberta with the triangulated unit cube uses a limit derived from this
// since albertas refinement algorithm seem to produce bad results)
const double conv_limit = 0.85;

// stop refining after the grid has more than this many elements (that means
// in 3D that the fine grid may have up to 8 times as may elements)
const unsigned maxelements = 2<<7;

// multiplier for the stepsize obtained from the FDTD criterion
const double stepadjust = 0.25;

// how long to run the simulation
const double duration = 1/c0;

// whether to run the check for all levels of refinement.  If false, will run
// the check only for the first and the last level.  Running the check for all
// levels is useful when debugging.
const bool do_all_levels = true;

// Order of quadrature rules to use
const unsigned int quadrature_order = 3;

// Whether to write computed and exact solution as vtk files after each step
const bool do_vtk_output = false;

// where to place a probe to measure the E-field
const double probe_location[] = {0.5, 0.5, 0.5};

// probe_location as FieldVector, initialized from probe_location in main()
Dune::FieldVector<double, 3> probe_location_fv;

//
//  CODE
//

//===============================================================
// Define parameter function mu
//===============================================================

// function for defining the source term
template<typename GV, typename RF>
class ConstFunc
  : public Dune::PDELab::AnalyticGridFunctionBase<
      Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
      ConstFunc<GV,RF>
    >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,ConstFunc<GV,RF> > BaseT;

  ConstFunc (const GV& gv, RF val_ = 1)
    : BaseT(gv)
    , val(val_)
  {}

  inline void
  evaluateGlobal (const typename Traits::DomainType& x, 
                  typename Traits::RangeType& y) const
  {
    y = val;
  }

private:
  RF val;
};

// boundary grid function selecting boundary conditions 
template<typename GV>
class B
  : public Dune::PDELab::BoundaryGridFunctionBase<
      Dune::PDELab::BoundaryGridFunctionTraits<
        GV,
        int,1,Dune::FieldVector<int,1>
      >,
      B<GV>
    >
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B<GV> > BaseT;

  B (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void
  evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
            const typename Traits::DomainType& x,
            typename Traits::RangeType& y) const
  {
    y = 1; // Dirichlet everywhere
  }

  //! get a reference to the GridView
  inline const GV& getGridView () const
  {
    return gv;
  }
};

// function for initialization
template<typename GV, typename RF>
class Init
  : public Dune::PDELab::AnalyticGridFunctionBase<
      Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,3>,
      Init<GV,RF>
    >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,3> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Init<GV,RF> > BaseT;

  Init (const GV& gv, typename Traits::DomainFieldType radius = .1)
    : BaseT(gv), center(0.5), radius2(radius*radius)
  { }

  inline void
  evaluateGlobal (const typename Traits::DomainType& x, 
                  typename Traits::RangeType& y) const
  {
    y = 0;
    if((x-center).two_norm2() < radius2)
      y[0] = -1;
  }

private:
  typename Traits::DomainType center;
  typename Traits::DomainFieldType radius2;
};

//! calculate length of smallest edge.  If no edge was found, return +inf
template<typename GV>
typename GV::Grid::ctype smallestEdge(const GV& gv)
{
  typedef typename GV::Grid::ctype ct;
  static const int dim = GV::dimension;
  static const int dimw = GV::dimensionworld;

  ct min = std::numeric_limits<ct>::infinity();
  typedef Dune::FieldVector<ct, dimw> CV;
  typedef typename GV::template Codim<0>::Iterator Iterator;
  const Iterator end = gv.template end<0>();

  for(Iterator it = gv.template begin<0>(); it!=end; ++it) {
    const Dune::GenericReferenceElement<ct, dim> &refElem
      = Dune::GenericReferenceElements<ct, dim>::general(it->type());

    for(int i = 0; i < refElem.size(dim-1); ++i) {
      CV diff = it->geometry().corner(refElem.subEntity(i,dim-1, 0,dim));
      diff   -= it->geometry().corner(refElem.subEntity(i,dim-1, 1,dim));
      ct distance = diff.two_norm();
      if(distance < min)
        min = distance;
    }
  }
  return min;
}


//===============================================================
// Problem setup and solution 
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename CON, typename ReferenceFactory,
         typename Probe> 
double electrodynamic (const GV& gv, const FEM& fem, unsigned integrationOrder,
                       const ReferenceFactory &referenceFactory,
                       double Delta_t, unsigned steps,
                       const std::string filename, Probe &probe)
{
  // constants and types
  typedef typename GV::Grid::ctype DF;
  typedef typename FEM::Traits::LocalFiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType RangeField;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS; 
  GFS gfs(gv,fem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<RangeField>::Type C;
  C cg;
  cg.clear();
  typedef B<GV> BType;
  BType b(gv);
  Dune::PDELab::constraints(b,gfs,cg);

  // make coefficent Vector and initialize it from a function
  typedef typename GFS::template VectorContainer<RangeField>::Type V;
  typedef Dune::PDELab::DiscreteGridFunctionGlobal<GFS,V> DGF;

  //initial solution for time-step -1
  Dune::SmartPointer<V> xprev(new V(gfs));
  *xprev = 0.0;
  Dune::PDELab::interpolateGlobal(*referenceFactory.function(gv, -Delta_t),gfs,*xprev);
  // actually, it is an error if any constrained dof is != 0, but set it here nevertheless
  Dune::PDELab::set_constrained_dofs(cg,0.0,*xprev);

  //initial solution for time-step 0
  Dune::SmartPointer<V> xcur(new V(gfs));
  *xcur = 0.0;
  Dune::PDELab::interpolateGlobal(*referenceFactory.function(gv, 0),gfs,*xcur);
  // actually, it is an error if any constrained dof is != 0, but set it here nevertheless
  Dune::PDELab::set_constrained_dofs(cg,0.0,*xcur);

  Dune::SmartPointer<V> xnext(0);

  // we're using dirichlet 0 everywhere, simply leave everything as 0
  V affineShift(gfs);
  affineShift = 0.0;
  //Dune::PDELab::interpolateGlobal(g,gfs,affineShift);
  //Dune::PDELab::set_nonconstrained_dofs(cg,0.0,affineShift);

  // make grid function operator
  typedef ConstFunc<GV,RangeField> MuType;
  MuType mu(gv,mu0);
  typedef ConstFunc<GV,RangeField> EpsType;
  EpsType eps(gv,eps0);
  typedef Dune::PDELab::Electrodynamic<EpsType,MuType,V> LOP; 
  LOP lop(eps, mu, Delta_t, integrationOrder);
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    LOP,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // represent operator as a matrix
  typedef typename GOS::template MatrixContainer<RangeField>::Type M;
  M m(gos);
  m = 0;
  // just use some values here
  lop.setEprev(*xprev);
  lop.setEcur(*xcur);
  gos.jacobian(affineShift,m);
//   Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);
  Dune::Richardson<V,V> prec(1.0);
//   Dune::SeqILU0<M,V,V> prec(m,1.0);
  Dune::MatrixAdapter<M,V,V> op(m);
  Dune::CGSolver<V> solver(op,prec,1E-10,1000,0);

//   typedef Dune::SuperLU<typename M::BaseT> Solver;
//   Solver solver(m);

  Dune::SmartPointer<Dune::VTKSequenceWriter<GV> > vtkwriter(0);
  if(do_vtk_output)
    vtkwriter = new Dune::VTKSequenceWriter<GV>(gv,filename,".","", Dune::VTKOptions::nonconforming);

  probe.measure(DGF(gfs, *xprev), -Delta_t);
//  std::cout << "u[-1]\n" << *xprev << std::endl;
  probe.measure(DGF(gfs, *xcur), 0);
//   std::cout << "u[0]\n" << *xcur << std::endl;

  double errsum = 0;

  std::cout << "Number of steps " << steps << std::endl;
  for(unsigned step = 1; step <= steps; ++step) {
    std::cout << "Doing step " << step << std::endl;


    lop.setEprev(*xprev);
    lop.setEcur(*xcur);

    // evaluate residual w.r.t initial guess
    V rhs(gfs);
    rhs = 0.0;
    gos.residual(affineShift,rhs);
    rhs *= -1;

//     std::cout << "RHS\n" << rhs << std::endl;

    Dune::InverseOperatorResult stat;

    xnext = new V(gfs);
    *xnext = *xcur; // copy values to have a better start value for the solver
    *xnext -= affineShift;
    solver.apply(*xnext,rhs,stat);
    std::cout << "Solver results:" << std::endl
              << "  iterations: " << stat.iterations << std::endl
              << "  reduction:  " << stat.reduction << std::endl
              << "  converged:  " << stat.converged << std::endl
              << "  conv_rate:  " << stat.conv_rate << std::endl
              << "  elapsed:    " << stat.elapsed << std::endl;
    if(!stat.converged)
      DUNE_THROW(Dune::Exception, "Solver did not converge");

    *xnext += affineShift;

    errsum += l2difference2(gv,DGF(gfs,*xnext),*referenceFactory.function(gv, steps*Delta_t),integrationOrder);
    probe.measure(DGF(gfs, *xnext), Delta_t*step);
//     std::cout << "u[" << step << "]\n" << *xnext << std::endl;

    xprev = xcur;
    xcur = xnext;
  }
  
  probe.measureFinal(DGF(gfs, *xcur), Delta_t*steps);

  return std::sqrt(errsum/steps);
}

template<typename P>
struct GridPtrTraits;
template<typename G>
struct GridPtrTraits<Dune::SmartPointer<G> > {
  typedef G Grid;
};
template<typename G>
struct GridPtrTraits<Dune::GridPtr<G> > {
  typedef G Grid;
};


template<typename GV, typename LPF>
void testLevel(const GV &gv, unsigned level, const std::string &prefix, std::ostream &dat,
               LPF &lpf, double &error, double &mean_h)
{
  typedef Dune::PDELab::EdgeS03DLocalFiniteElementMap<GV, double> FEM;
  std::ostringstream levelprefix;
  levelprefix << prefix << ".level" << level;

  typedef typename LPF::template Traits<GV>::Probe Probe;
  Dune::SmartPointer<Probe> probe(lpf.getProbe(gv, level));

  std::cout << "electrodynamic level " << level << std::endl;
  // time step
  double Delta_t = stepadjust*smallestEdge(gv)/std::sqrt(double(GV::dimension))/c0;
  unsigned steps = duration/Delta_t;
  if(Delta_t*steps < duration) {
    //adjust such that duration is hit exactly
    ++steps;
    Delta_t = duration/steps;
  }
  error = electrodynamic
    <GV,FEM,Dune::PDELab::OverlappingConformingDirichletConstraints,ResonatorSolutionFactory<GV,double>,Probe>
    (gv, FEM(gv), quadrature_order,
     ResonatorSolutionFactory<GV,double>(),
     Delta_t, steps, levelprefix.str(), *probe);
  mean_h = std::pow(1/double(gv.size(0)), 1/double(GV::dimension));
  std::cout << "L2 error: " 
            << std::setw(8) << gv.size(0) << " elements, <h>=" 
            << std::scientific << mean_h << ", error="
            << std::scientific << error << std::endl;
  dat << mean_h << "\t" << error << std::endl;
}
template<typename GV, typename LPF>
void testLevel(const GV &gv, unsigned level, const std::string &prefix, std::ostream &dat,
               LPF &lpf)
{
  double error, mean_h;
  testLevel(gv, level, prefix, dat, lpf, error, mean_h);
}


template<typename Grid, typename GPF>
void test(Grid &grid, int &result, GnuplotGraph &graph, GPF &gpf,
          double conv_limit, std::string name = "")
{
  typedef typename Grid::LeafGridView GV;

  if(name == "") name = grid.name();

  typedef typename GPF::template Traits<Grid>::LevelProbeFactory LPF;
  Dune::SmartPointer<LPF> lpf(gpf.levelProbeFactory(grid, name));

  std::cout << std::endl
            << "Testing Electrodynamic problem with EdgeS03D and " << name << std::endl;

  std::string filename = "electrodynamic-" + name;
  graph.addPlot("'" + graph.datname() + "' title '" + name + "' with linespoints");

  graph.dat() << "#<h>\terror" << std::endl;
  graph.dat().precision(8);

  unsigned level = 0;
  double error0, mean_h0;
  testLevel(grid.leafView(), level, filename, graph.dat(), *lpf, error0, mean_h0);

  while(true) {
    ++level;
    grid.globalRefine(1);
    if(unsigned(grid.leafView().size(0)) >= maxelements) break;

    if(do_all_levels)
      testLevel(grid.leafView(), level, filename, graph.dat(), *lpf);
  }
  
  double errorf, mean_hf;
  testLevel(grid.leafView(), level, filename, graph.dat(), *lpf, errorf, mean_hf);

  double total_convergence = std::log(errorf/error0)/std::log(mean_hf/mean_h0);
  std::cout << "electrodynamic total convergence: "
            << std::scientific << total_convergence << std::endl;

  if(result != 1)
    result = 0;

  if(total_convergence < conv_limit) {
    std::cout << "Error: electrodynamic total convergence < " << conv_limit << std::endl;
    result = 1;
  }
}


//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);
  if(mpiHelper.rank() == 0)
    std::cout << "Number of processes: " << mpiHelper.size() << std::endl;
  std::cout << "Rank: " << mpiHelper.rank() << std::endl;
//   if(mpiHelper.rank() == 1) {
//     int i = 0;
//     std::cout << "PID " << getpid() << " ready for attach\n", getpid();
//     while (0 == i)
//       sleep(5);
//   }

  try{
    // 77 is special and means "test was skipped".  Return that if non of the
    // supported grids were available
    int result = 77;

    // Initialize probe_location_fv from probe_location
    for(unsigned i = 0; i < 3; ++i)
      probe_location_fv[i] = probe_location[i];
    typedef Dune::PDELab::GnuplotGridProbeFactory<Dune::FieldVector<double, 3> > PF1;
    Dune::SmartPointer<PF1> pf1 = new PF1("electrodynamic-probe", probe_location_fv);

    typedef ResonatorVTKGridProbeFactory<double> PF2;
    Dune::SmartPointer<PF2> pf2 = new PF2("electrodynamic");

    typedef ResonatorGlobalErrorGridProbeFactory<double> PF3;
    Dune::SmartPointer<PF3> pf3 = new PF3("electrodynamic-globalerror", quadrature_order);

    typedef Dune::PDELab::GridProbeFactoryListTraits<PF1, PF2, PF3>::GPF PFList;
    Dune::SmartPointer<PFList> pflist
      = makeGridProbeFactoryList(pf1, pf2, pf3);
      
    GnuplotGraph graph("testelectrodynamic");
    graph.addCommand("set logscale xy");
    graph.addCommand("set xlabel '<h>'");
    graph.addCommand("set ylabel 'L2 error'");
    graph.addCommand("set key left top reverse Left");
    graph.addCommand("");
    graph.addCommand("set terminal postscript eps color solid");
    graph.addCommand("set output 'electrodynamic.eps'");
    graph.addCommand("");

#ifdef HAVE_ALBERTA
#if (ALBERTA_DIM != 3)
#error ALBERTA_DIM is not set to 3 -- please check the Makefile.am
#endif
//     test(*UnitTetrahedronMaker         <Dune::AlbertaGrid<3, 3>    >::create(),
//          result, graph, conv_limit,    "alberta-tetrahedron");
//     test(*KuhnTriangulatedUnitCubeMaker<Dune::AlbertaGrid<3, 3>    >::create(),
//          result, graph, .7*conv_limit, "alberta-triangulated-cube-6");
//     {
//       Dune::GridPtr<Dune::AlbertaGrid<3, 3> > gridptr("grids/brick.dgf");
//       test(*gridptr,
//            result, graph, conv_limit,    "alu-triangulated-brick-6");
//     }
#endif

#ifdef HAVE_ALUGRID
//     test(*UnitTetrahedronMaker         <Dune::ALUSimplexGrid<3, 3> >::create(),
//          result, graph, conv_limit,    "alu-tetrahedron");
//     test(*KuhnTriangulatedUnitCubeMaker<Dune::ALUSimplexGrid<3, 3> >::create(),
//          result, graph, conv_limit,    "alu-triangulated-cube-6");
#endif // HAVE_ALUGRID

#ifdef HAVE_UG
//     test(*UnitTetrahedronMaker         <Dune::UGGrid<3>            >::create(),
//          result, graph, conv_limit,    "ug-tetrahedron");
    test(*KuhnTriangulatedUnitCubeMaker<Dune::UGGrid<3>            >::create(),
         result, graph, *pflist, conv_limit,    "ug-triangulated-cube-6");
#endif // HAVE_ALBERTA
	return result;
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
