// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

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
#include "../localoperator/electrostatic.hh"

#include "gnuplotgraph.hh"
#include "gridexamples.hh"
#include "l2difference.hh"


//===============================================================
//===============================================================
// Solve the Electrostatic "wave" equation
//           rot(1/mu * rot E) = 0 in \Omega, 
//                 n x (n x E) = g on \partial\Omega_D
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
const unsigned maxelements = 10000;

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
    y = 1; // Dirichlet
  }

  //! get a reference to the GridView
  inline const GV& getGridView ()
  {
    return gv;
  }
};

// function for Dirichlet boundary conditions and initialization
template<typename GV, typename RF>
class G
  : public Dune::PDELab::AnalyticGridFunctionBase<
      Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,3>,
      G<GV,RF>
    >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,3> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G<GV,RF> > BaseT;

  G (const GV& gv)
    : BaseT(gv)
  {
    prescribedE[0] = 1;
    prescribedE[1] = 0;
    prescribedE[2] = 0;
  }

  inline void
  evaluateGlobal (const typename Traits::DomainType& x, 
                  typename Traits::RangeType& y) const
  {
    typename Traits::DomainType center(0.5); center[0] = -0.5;
    y = x;
    y -= center;
    y /= std::pow(y.two_norm(),3);

//     y = prescribedE;
  }

private:
  typename Traits::RangeType prescribedE;
};

//===============================================================
// Problem setup and solution 
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename CON, int q> 
double electrostatic (const GV& gv, const FEM& fem, const std::string &filename = "")
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
  V x0(gfs);
  x0 = 0.0;
  typedef G<GV,R> GType;
  GType g(gv);
  Dune::PDELab::interpolateGlobal(g,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);

  // make grid function operator
  typedef Mu<GV,R> MuType;
  MuType mu(gv);
  typedef Dune::PDELab::Electrostatic<MuType,q> LOP; 
  LOP lop(mu);
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    LOP,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // represent operator as a matrix
  typedef typename GOS::template MatrixContainer<R>::Type M;
  M m(gos);
  m = 0;
  gos.jacobian(x0,m);
  //Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // evaluate residual w.r.t initial guess
  V rhs(gfs);
  rhs = 0.0;
  gos.residual(x0,rhs);
  rhs *= -1;

  typedef Dune::SuperLU<typename M::BaseT> SLUSolver;
  SLUSolver sluSolver(m);
  Dune::InverseOperatorResult stat;

  V x(gfs);
  x = 0.0;
  sluSolver.apply(x,rhs,stat);
  std::cout << "SuperLU results:" << std::endl
            << "  iterations: " << stat.iterations << std::endl
            << "  reduction:  " << stat.reduction << std::endl
            << "  converged:  " << stat.converged << std::endl
            << "  conv_rate:  " << stat.conv_rate << std::endl
            << "  elapsed:    " << stat.elapsed << std::endl;
  x += x0;  // set constrained dofs???

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunctionGlobal<GFS,V> DGF;
  DGF dgf(gfs,x);
  
  if(filename != "") {
    // output grid function with VTKWriter
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,0);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.write(filename,Dune::VTKOptions::ascii);
  }

  return l2difference(gv,g,dgf,4);
}

template<typename Grid>
void test(Dune::SmartPointer<Grid> grid, int &result, GnuplotGraph &graph, double conv_limit, std::string name = "")
{
  typedef typename Grid::LeafGridView GV;

  if(name == "") name = grid->name();

  std::cout << std::endl
            << "Testing Electrostatic problem with EdgeS03D and " << name << std::endl;

  std::string filename = "electrostatic-" + name;
  graph.addPlot("'" + graph.datname() + "' title '" + name + "' with linespoints");

  graph.dat() << "#h\terror" << std::endl;
  graph.dat().precision(8);

  typedef Dune::PDELab::EdgeS03DLocalFiniteElementMap<typename Grid::LeafGridView, double> FEM;

  std::cout << "electrostatic level 0" << std::endl;
  double error0 = electrostatic
    <GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>
    (grid->leafView(), FEM(grid->leafView()), filename+"-coarse");
  double mean_h0 = std::pow(1/double(grid->leafView().size(0)), 1/double(Grid::dimension));
  std::cout << "L2 error: " 
            << std::setw(8) << grid->leafView().size(0) << " elements, <h>=" 
            << std::scientific << mean_h0 << ", error="
            << std::scientific << error0 << std::endl;
  graph.dat() << mean_h0 << "\t" << error0 << std::endl;

//   if(Dune::FloatCmp::eq(error0, 0.0)) {
//     std::cerr << "Error: The analytic function was perfectly interpolated." << std::endl
//               << "Error: This makes the interpolation convergence test meaningless." << std::endl
//               << "Error: Please change this test program to use an analyting function which cannot be" << std::endl
//               << "Error: represented exacly by this basis, i.e. something containing exp()" << std::endl;
//     result = 1;
//     return;
//   }

  while(1) {
    grid->globalRefine(1);

    if((unsigned int)(grid->leafView().size(0)) >= maxelements)
      break;

    if(measure_after_every_refinement) {
      std::cout << "electrostatic level " << grid->maxLevel() << std::endl;
      double error = electrostatic
        <GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>
        (grid->leafView(), FEM(grid->leafView()));
      double mean_h = std::pow(1/double(grid->leafView().size(0)), 1/double(Grid::dimension));
      std::cout << "L2 error: " 
                << std::setw(8) << grid->leafView().size(0) << " elements, <h>=" 
                << std::scientific << mean_h << ", error="
                << std::scientific << error << std::endl;
      graph.dat() << mean_h << "\t" << error << std::endl;
    }
  }

  std::cout << "electrostatic level " << grid->maxLevel() << std::endl;
  double errorf = electrostatic
    <GV,FEM,Dune::PDELab::ConformingDirichletConstraints,2>
    (grid->leafView(), FEM(grid->leafView()), filename+"-fine");
  double mean_hf = std::pow(1/double(grid->leafView().size(0)), 1/double(Grid::dimension));
  std::cout << "L2 error: " 
            << std::setw(8) << grid->leafView().size(0) << " elements, <h>=" 
            << std::scientific << mean_hf << ", error="
            << std::scientific << errorf << std::endl;
  graph.dat() << mean_hf << "\t" << errorf << std::endl;

  double total_convergence = std::log(errorf/error0)/std::log(mean_hf/mean_h0);
  std::cout << "electrostatic total convergence: "
            << std::scientific << total_convergence << std::endl;

  if(result != 1)
    result = 0;

  if(total_convergence < conv_limit) {
    std::cout << "Error: electrostatic total convergence < " << conv_limit << std::endl;
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

    GnuplotGraph graph("testelectrostatic");
    graph.addCommand("set logscale xy");
    graph.addCommand("set xlabel '<h>'");
    graph.addCommand("set ylabel 'L2 error'");
    graph.addCommand("set key left top reverse Left");
    graph.addCommand("");
    graph.addCommand("set terminal postscript eps color solid");
    graph.addCommand("set output 'electrostatic.eps'");
    graph.addCommand("");

#ifdef HAVE_ALBERTA
#if (ALBERTA_DIM != 3)
#error ALBERTA_DIM is not set to 3 -- please check the Makefile.am
#endif
//     test(UnitTetrahedronMaker         <Dune::AlbertaGrid<3, 3>    >::create(),
//          result, graph, conv_limit,    "alberta-tetrahedron");
    test(KuhnTriangulatedUnitCubeMaker<Dune::AlbertaGrid<3, 3>    >::create(),
         result, graph, .7*conv_limit, "alberta-triangulated-cube-6");
#endif

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
