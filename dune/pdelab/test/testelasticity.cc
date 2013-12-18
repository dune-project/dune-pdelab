// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <vector>
#include <map>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/linearelasticity.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>

#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

static const double szX = 10.0;

// define some grid functions to interpolate from
template<typename GV, typename RF>
class G
    : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension>,
                                                  G<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G<GV,RF> > BaseT;

  G (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    y = 0.0;
    if (x[0] == szX)
    {
        y[0] = 0.1;
    }
  }
};

// define some boundary grid functions to define boundary conditions
class BCType
  : public Dune::PDELab::DirichletConstraintsParameters
{
public:
  template<typename I>
  bool isDirichlet(const I & intersection,
    const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension>
      xg = intersection.geometry().global( coord );
    return (xg[0] == 0.0); // || xg[0] == szX);  // Dirichlet b.c. on left & right boundary
  }
};

// generate a P1 function and output it
template<class GV>
void testp1 (const GV& gv, double mu, double lambda, double constG)
{
  typedef typename GV::Grid::ctype DF;

  const int dim = GV::dimension;

  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
  FEM fem(gv);

  // make function space
  typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
  typedef Dune::PDELab::ISTLVectorBackend<> ComponentVectorBackend;

  typedef Dune::PDELab::DefaultLeafOrderingTag Mapper;

  typedef Dune::PDELab::LexicographicOrderingTag OrderingTag;
  typedef Dune::PDELab::VectorGridFunctionSpace<
    GV,
    FEM,
    dim,
    Dune::PDELab::ISTLVectorBackend<>,
    ComponentVectorBackend,
    Constraints,
    OrderingTag,
    Mapper
    > GFS;
  GFS gfs(gv,fem);
  gfs.name("displacement");

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<double>::Type C;
  C cg;
  cg.clear();
  // Dune::PDELab::PowerConstraintsParameters<BCType,dim> b;
  BCType b;
  Dune::PDELab::constraints(b,gfs,cg);

  std::cout << gfs.size() << " DOFs\n";
  std::cout << cg.size() << " constraint DOFs\n";

  // make local operator
  typedef Dune::PDELab::LinearElasticity LO;
  LO lo(mu, lambda, constG);

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(9); // 2D Q1

  // make grid operator
  typedef Dune::PDELab::GridOperator<
    GFS,GFS,LO,
    MBE,
    double,double,double,
    C,C> GOS;
  GOS gos(gfs,cg,gfs,cg,lo,mbe);

  // make coefficent Vector and initialize it from a function
  typedef typename GOS::Traits::Domain V;
  V x0(gfs);
  x0 = 0.0;
  typedef G<GV,double> GType;
  GType g(gv);
  Dune::PDELab::interpolate(g,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);

  // represent operator as a matrix
  typedef typename GOS::Traits::Jacobian M;
  M m(gos);
  m = 0.0;
  gos.jacobian(x0,m);
  if (gfs.size() <= 16)
    Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,3);

  // evaluate residual w.r.t initial guess
  V r(gfs);
  r = 0.0;
  gos.residual(x0,r);

  // make ISTL solver
  typedef typename M::BaseT ISTL_M;
  typedef typename V::BaseT ISTL_V;
  Dune::MatrixAdapter<ISTL_M,ISTL_V,ISTL_V> opa(m.base());
  Dune::SeqILU0<ISTL_M,ISTL_V,ISTL_V> ilu0(m.base(),1e-2);

  Dune::CGSolver<ISTL_V> solver(opa,ilu0,1E-20,5000,2);
  Dune::InverseOperatorResult stat;

  // solve the jacobian system
  r *= -1.0; // need -residual
  V x(gfs,0.0);
  Dune::PDELab::set_nonconstrained_dofs(cg,1.0,x);
  solver.apply(x.base(),r.base(),stat);
  x += x0;

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x);
  vtkwriter.write("testelasticity",Dune::VTK::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    double mu = 1.0;
    double lambda = 1.0;
    double g = 1.0;
    int level=5;

    if (argc > 1)
        mu = atof(argv[1]);
    if (argc > 2)
        lambda = atof(argv[2]);
    if (argc > 3)
        g = atof(argv[3]);
    if (argc > 4)
        level = atoi(argv[4]);

    std::cout << "mu     = " << mu << "\n"
              << "lambda = " << lambda << "\n"
              << "g = " << g << std::endl;

    Dune::FieldVector<double,2> L(1); L[0] = szX;
    Dune::array<int,2> N(Dune::fill_array<int,2>(1)); N[0] = szX;
    std::bitset<2> B(false);
    Dune::YaspGrid<2> grid(L,N,B,0);
    grid.globalRefine(level);

    testp1(grid.leafView(), mu, lambda, g);

    // test passed
    return 0;
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
