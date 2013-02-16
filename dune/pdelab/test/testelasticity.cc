// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<vector>
#include<map>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/shared_ptr.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>

#include"../finiteelementmap/q1fem.hh"
#include"../finiteelementmap/conformingconstraints.hh"
#include"../gridfunctionspace/gridfunctionspace.hh"
#include"../gridfunctionspace/gridfunctionspaceutilities.hh"
#include"../gridfunctionspace/interpolate.hh"
#include"../constraints/constraints.hh"
#include"../common/function.hh"
#include"../common/vtkexport.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"../localoperator/linearelasticity.hh"
#include"../backend/istlvectorbackend.hh"
#include"../backend/istlmatrixbackend.hh"
#include"../backend/istlsolverbackend.hh"

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
  typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,double, dim> FEM;
  FEM fem;
  
  // make function space
  typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
  typedef Dune::PDELab::ISTLVectorBackend<1> VectorBackend;
  typedef Dune::PDELab::SimpleGridFunctionStaticSize GFSSize;

  typedef Dune::PDELab::GridFunctionSpace
    <GV, FEM, Constraints, VectorBackend, GFSSize> Q1GFS;
  Q1GFS q1gfs(gv,fem);

  typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper GFMapper;
  typedef Dune::PDELab::PowerGridFunctionSpace<Q1GFS,dim,GFMapper> GFS;
  GFS gfs(q1gfs);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<double>::Type C;
  C cg;
  cg.clear();
  // Dune::PDELab::PowerConstraintsParameters<BCType,dim> b;
  BCType b;
  Dune::PDELab::constraints(b,gfs,cg);

  std::cout << gfs.size() << " DOFs\n";
  std::cout << cg.size() << " constraint DOFs\n";

  // make coefficent Vector and initialize it from a function
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,double>::Type V;
  V x0(gfs);
  x0 = 0.0;
  typedef G<GV,double> GType;
  GType g(gv);
  Dune::PDELab::interpolate(g,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);

  // make grid function operator
  Dune::PDELab::LinearElasticity la(mu, lambda, constG);
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    Dune::PDELab::LinearElasticity,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(gfs,cg,gfs,cg,la);

  // represent operator as a matrix
  typedef typename GOS::template MatrixContainer<double>::Type M;
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
  Dune::MatrixAdapter<M,V,V> opa(m);
  Dune::SeqILU0<M,V,V> ilu0(m,1e-2);

  Dune::CGSolver<V> solver(opa,ilu0,1E-20,5000,2);
  Dune::InverseOperatorResult stat;

  // solve the jacobian system
  r *= -1.0; // need -residual
  V x(gfs,0.0);
  Dune::PDELab::set_nonconstrained_dofs(cg,1.0,x);
  solver.apply(x,r,stat);
  x += x0;

  // make discrete function object
  typedef Dune::PDELab::VectorDiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x);
  
  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"displacement"));
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
    Dune::FieldVector<int,2> N(1); N[0] = szX;
    Dune::FieldVector<bool,2> B(false);
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
