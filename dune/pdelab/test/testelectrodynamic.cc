// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#define DISABLE_LOCAL_INTERFACE

#include <iostream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/smartpointer.hh>

#include <dune/grid/albertagrid/agrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>

#include "../backend/istlmatrixbackend.hh"
#include "../backend/istlsolverbackend.hh"
#include "../backend/istlvectorbackend.hh"
#include "../common/function.hh"
#include "../common/geometrywrapper.hh"
#include "../common/vtkexport.hh"
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
const unsigned maxelements = 1000;

// whether to measure the error after every refinement, or just once at the
// beginning and once at the end
const bool measure_after_every_refinement = false;

//
//  CODE
//

//===============================================================
// Define parameter function mu
//===============================================================

// function for defining the source term
template<typename GV, typename RF>
class Mu
  : public Dune::PDELab::AnalyticGridFunctionBase<
      Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
      Mu<GV,RF>
    >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Mu<GV,RF> > BaseT;

  Mu (const GV& gv) : BaseT(gv) {}

  inline void
  evaluateGlobal (const typename Traits::DomainType& x, 
                  typename Traits::RangeType& y) const
  {
//     typename Traits::DomainType myx(0.5);
//     myx -= x;
//     y = 1+exp(-3.0*myx.two_norm2());
    y = 1;
  }

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
  inline const GV& getGridView ()
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
      y[0] = 1;
  }

private:
  typename Traits::DomainType center;
  typename Traits::DomainFieldType radius2;
};

//! calculate length of smallest edge.  If not edge was found, return +inf
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
template<typename GV, typename FEM, typename CON, int q> 
double electrodynamic (const GV& gv, const FEM& fem, double Delta_t, unsigned steps, const std::string &filename = "")
{
  // constants and types
  typedef typename GV::Grid::ctype DF;
  typedef typename FEM::Traits::LocalFiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS; 
  GFS gfs(gv,fem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<R>::Type C;
  C cg;
  cg.clear();
  typedef B<GV> BType;
  BType b(gv);
  Dune::PDELab::constraints(b,gfs,cg);

  // make coefficent Vector and initialize it from a function
  typedef typename GFS::template VectorContainer<R>::Type V;
  Dune::SmartPointer<V> xprev(new V(gfs));
  *xprev = 0.0;
  Dune::PDELab::interpolateGlobal(Init<GV,double>(gv,1),gfs,*xprev);
  //Dune::PDELab::set_nonconstrained_dofs(cg,1.0,*xprev);
  Dune::PDELab::set_constrained_dofs(cg,0.0,*xprev);
  // we're using dirichlet 0 everywhere, simply leave everything as 0
  //Dune::PDELab::interpolateGlobal(g,gfs,*xprev);
  //Dune::PDELab::set_nonconstrained_dofs(cg,0.0,*xprev);
  Dune::SmartPointer<V> xcur(xprev);
  Dune::SmartPointer<V> xnext(0);

  // make grid function operator
  typedef Mu<GV,R> MuType;
  MuType mu(gv);
  typedef Dune::PDELab::Electrodynamic<MuType,V,q> LOP; 
  LOP lop(mu, Delta_t);
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    LOP,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // represent operator as a matrix
  typedef typename GOS::template MatrixContainer<R>::Type M;
  M m(gos);
  m = 0;
  // just use some values here
  lop.setEprev(*xprev);
  lop.setEcur(*xcur);
  gos.jacobian(*xcur,m);
//   Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);
//   Dune::Richardson<V,V> prec(1.0);
  Dune::SeqILU0<M,V,V> prec(m,1.0);
  Dune::MatrixAdapter<M,V,V> op(m);
  Dune::CGSolver<V> solver(op,prec,1E-10,5000,0);

//   typedef Dune::SuperLU<typename M::BaseT> Solver;
//   Solver solver(m);

  // used to get the rhs from the residual
  V zero(gfs);
  zero = 0.0;

  std::cout << "Number of steps " << steps << std::endl;
  for(unsigned step = 0; step < steps; ++step) {
    if(filename != "") {
      // make discrete function object
      typedef Dune::PDELab::DiscreteGridFunctionGlobal<GFS,V> DGF;
      DGF dgf(gfs,*xcur);

      std::ostringstream s;
      s << filename << "." << step;

      // output grid function with VTKWriter
      Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,0);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
      vtkwriter.write(s.str(),Dune::VTKOptions::ascii);
    }

    std::cout << "Doing step " << step << std::endl;


    lop.setEprev(*xprev);
    lop.setEcur(*xcur);

    // evaluate residual w.r.t initial guess
    V rhs(gfs);
    rhs = 0.0;
    gos.residual(zero,rhs);
    rhs *= -1;

    Dune::InverseOperatorResult stat;

    xnext = new V(gfs);
    *xnext = *xcur; // copy values to have a better start value for the solver
    solver.apply(*xnext,rhs,stat);
    std::cout << "Solver results:" << std::endl
              << "  iterations: " << stat.iterations << std::endl
              << "  reduction:  " << stat.reduction << std::endl
              << "  converged:  " << stat.converged << std::endl
              << "  conv_rate:  " << stat.conv_rate << std::endl
              << "  elapsed:    " << stat.elapsed << std::endl;
    if(!stat.converged)
      DUNE_THROW(Dune::Exception, "Solver did not converge");

    //x += x0;  // set constrained dofs

    xprev = xcur;
    xcur = xnext;
  }

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunctionGlobal<GFS,V> DGF;
  DGF dgf(gfs,*xcur);

  if(filename != "") {
    std::ostringstream s;
    s << filename << "." << steps;

    // output grid function with VTKWriter
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,0);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.write(s.str(),Dune::VTKOptions::ascii);
  }
  return 0; //l2difference(gv,g,dgf,4);
}

template<typename Grid>
void test(Dune::SmartPointer<Grid> grid, int &result, GnuplotGraph &graph, double conv_limit, std::string name = "")
{
  typedef typename Grid::LeafGridView GV;

  if(name == "") name = grid->name();

  std::cout << std::endl
            << "Testing Electrodynamic problem with EdgeS03D and " << name << std::endl;

  std::string filename = "electrodynamic-" + name;
  std::ostringstream plot;
  plot << "'" << filename << ".dat' title '" << name << "' with linespoints";
  graph.addPlot(plot.str());

  std::ofstream dat((filename+".dat").c_str());
  dat << "#h\terror" << std::endl;
  dat.precision(8);

  typedef Dune::PDELab::EdgeS03DLocalFiniteElementMap<typename Grid::LeafGridView, double> FEM;

  std::cout << "electrodynamic level 0" << std::endl;
  // time step
  double Delta_t = smallestEdge(grid->leafView())/std::sqrt(double(Grid::dimension));
  unsigned steps = 1/Delta_t;
  double error0 = electrodynamic
    <GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>
    (grid->leafView(), FEM(grid->leafView()), Delta_t, steps, filename+"-coarse");
  double mean_h0 = std::pow(1/double(grid->leafView().size(0)), 1/double(Grid::dimension));
  std::cout << "L2 error: " 
            << std::setw(8) << grid->leafView().size(0) << " elements, <h>=" 
            << std::scientific << mean_h0 << ", Delta t="
            << std::scientific << Delta_t << ", error="
            << std::scientific << error0 << std::endl;
  dat << mean_h0 << "\t" << error0 << std::endl;

  while(1) {
    grid->globalRefine(1);

    if((unsigned int)(grid->leafView().size(0)) >= maxelements)
      break;

    if(measure_after_every_refinement) {
      std::cout << "electrodynamic level " << grid->maxLevel() << std::endl;
      Delta_t = smallestEdge(grid->leafView())/std::sqrt(double(Grid::dimension));
      steps = 1/Delta_t;
      double error = electrodynamic
        <GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>
        (grid->leafView(), FEM(grid->leafView()), Delta_t, steps);
      double mean_h = std::pow(1/double(grid->leafView().size(0)), 1/double(Grid::dimension));
      std::cout << "L2 error: " 
                << std::setw(8) << grid->leafView().size(0) << " elements, <h>=" 
                << std::scientific << mean_h << ", Delta t="
                << std::scientific << Delta_t << ", error="
                << std::scientific << error << std::endl;
      dat << mean_h << "\t" << error << std::endl;
    }
  }

  std::cout << "electrodynamic level " << grid->maxLevel() << std::endl;
  Delta_t = smallestEdge(grid->leafView())/std::sqrt(double(Grid::dimension));
  steps = 1/Delta_t*10;
  double errorf = electrodynamic
    <GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>
    (grid->leafView(), FEM(grid->leafView()), Delta_t, steps, filename+"-fine");
  double mean_hf = std::pow(1/double(grid->leafView().size(0)), 1/double(Grid::dimension));
  std::cout << "L2 error: " 
            << std::setw(8) << grid->leafView().size(0) << " elements, <h>=" 
            << std::scientific << mean_hf << ", Delta t="
            << std::scientific << Delta_t << ", error="
            << std::scientific << errorf << std::endl;
  dat << mean_hf << "\t" << errorf << std::endl;

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
  try{
    // 77 is special and means "test was skipped".  Return that if non of the
    // supported grids were available
    int result = 77;

    GnuplotGraph graph("testelectrodynamic.gnuplot");
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
//     test(UnitTetrahedronMaker         <Dune::AlbertaGrid<3, 3>    >::create(),
//          result, graph, conv_limit,    "alberta-tetrahedron");
//     test(KuhnTriangulatedUnitCubeMaker<Dune::AlbertaGrid<3, 3>    >::create(),
//          result, graph, .7*conv_limit, "alberta-triangulated-cube-6");
#endif

#ifdef HAVE_ALUGRID
//     test(UnitTetrahedronMaker         <Dune::ALUSimplexGrid<3, 3> >::create(),
//          result, graph, conv_limit,    "alu-tetrahedron");
    test(KuhnTriangulatedUnitCubeMaker<Dune::ALUSimplexGrid<3, 3> >::create(),
         result, graph, conv_limit,    "alu-triangulated-cube-6");
#endif // HAVE_ALUGRID

#ifdef HAVE_UG
//     test(UnitTetrahedronMaker         <Dune::UGGrid<3>            >::create(),
//          result, graph, conv_limit,    "ug-tetrahedron");
//     test(KuhnTriangulatedUnitCubeMaker<Dune::UGGrid<3>            >::create(),
//          result, graph, conv_limit,    "ug-triangulated-cube-6");
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
